#!/usr/bin/env python3
"""Tissue WGBS: run inputs -> library-aware duplicate marking -> sample outputs.

No software installation or environment activation; commands use argv, never shell=True.
Profiles are documented defaults, not a claim of per-library empirical optimization.
"""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time

HERE=Path(__file__).resolve().parent
METADATA_DIR=HERE.parents[1]/"metadata/public_tissue_wgbs"
RUN_STEPS={'validate','fastq','trim','bsmap'}
REFERENCE_STEPS={'bsmap','mcall','wig','cna'}

def table(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as h:return list(csv.DictReader(h,delimiter='\t'))

def need(path):
    p=Path(path)
    if not p.is_file() or p.stat().st_size==0:raise ValueError('Missing/empty input: '+str(p))
    return p

def sig(path):
    p=need(path);s=p.stat();return [str(p.resolve()),s.st_size,s.st_mtime_ns]

def read_list(path):
    if not path:return []
    return [l.strip() for l in need(path).read_text().splitlines() if l.strip() and not l.lstrip().startswith('#')]

def load_manifest():
    samples=table(METADATA_DIR/'samples.tsv');runs=table(METADATA_DIR/'runs.tsv');profiles=table(METADATA_DIR/'profiles.tsv')
    for data,key in [(samples,'sample_id'),(runs,'run_id'),(profiles,'processing_group')]:
        ids=[x[key] for x in data]
        if len(ids)!=len(set(ids)):raise ValueError('Duplicate '+key)
        if any(not re.fullmatch(r'[A-Za-z0-9_.-]+',v) for v in ids):raise ValueError('Unsafe identifier in '+key)
    sm={s['sample_id']:s for s in samples};pm={p['processing_group']:p for p in profiles}
    for r in runs:
        if r['sample_id'] not in sm or r['processing_group'] not in pm:raise ValueError('Unmatched run '+r['run_id'])
        if r['layout'] not in {'SINGLE','PAIRED'}:raise ValueError('Invalid layout '+r['run_id'])
        if r['umi_status']!='NOT_REPORTED_NO_UMI_PROCESSING':raise ValueError('UMI workflow needs explicit implementation: '+r['run_id'])
    for p in profiles:
        for k in ['clip_r1','clip_r2','tail_r1','tail_r2','quality','min_length','mcall_pe1','mcall_pe2']:
            if int(p[k])<0:raise ValueError('Negative parameter '+k)
        if p['bsmap_n'] not in {'0','1'}:raise ValueError('Invalid BSMAP direction')
    chosen=read_list(os.environ.get('SAMPLE_LIST',''))
    if chosen and set(chosen)-set(sm):raise ValueError('Unknown sample IDs: '+str(set(chosen)-set(sm)))
    available=read_list(os.environ.get('AVAILABLE_RUN_LIST',''))
    if available and set(available)-{r['run_id'] for r in runs}:raise ValueError('Unknown available Run ID')
    selected=[s for s in samples if s['enabled']=='1' and (not chosen or s['sample_id'] in chosen)]
    allowed_cell_lines={'SAMN04889842','SAMN04889845'}
    if any(s['material']!='tissue' and not (s['material']=='cell_line' and s['biosample'] in allowed_cell_lines) for s in selected):
        raise ValueError('Only tissue and the two explicitly included cell lines are allowed')
    if available:
        selected=[s for s in selected if {r['run_id'] for r in runs if r['sample_id']==s['sample_id']}<=set(available)]
    ids={s['sample_id'] for s in selected};rr=[r for r in runs if r['sample_id'] in ids]
    for s in selected:
        sr=[r for r in rr if r['sample_id']==s['sample_id']]
        if not sr:raise ValueError('No runs for '+s['sample_id'])
        if len({r['layout'] for r in sr})!=1 or len({r['processing_group'] for r in sr})!=1:
            raise ValueError('Mixed layout/profile requires explicit split libraries: '+s['sample_id'])
    if not selected:raise ValueError('No complete enabled samples selected')
    return selected,rr,pm

def cutadapt_args(profile,paired,inputs,outputs,report,threads):
    p=profile
    cmd=['cutadapt','-j',str(threads),'-q',str(p['quality']),'-m',str(p['min_length']),'--json',str(report)]
    for flag,front,tail in [('-u','clip_r1','tail_r1')]+([('-U','clip_r2','tail_r2')] if paired else []):
        if int(p[front]):cmd += [flag,str(p[front])]
        if int(p[tail]):cmd += [flag,str(-int(p[tail]))]
    if paired:cmd+=['--pair-filter=any']
    cmd+=['-o',str(outputs[0])]
    if paired:cmd+=['-p',str(outputs[1])]
    return cmd+list(map(str,inputs))

def dump_outputs(folder,layout):
    """Resolve a fresh fasterq output directory without assuming its filename stem.

    Different tool versions can retain .sra in the stem. Never accept an
    ambiguous pair or silently treat a single/interleaved file as paired data.
    """
    folder=Path(folder)
    files=sorted(p for p in folder.iterdir() if p.is_file() and re.search(r'\.(?:fastq|fq)(?:\.gz)?$',p.name))
    if layout=='SINGLE':
        if len(files)!=1:raise ValueError('Expected one SINGLE FASTQ; found '+str([p.name for p in files]))
        return [need(files[0])],[]
    groups={};single=[]
    for p in files:
        match=re.fullmatch(r'(.+)_([12])\.(fastq|fq)(\.gz)?',p.name)
        if match:groups.setdefault(match[1],{}).setdefault(match[2],[]).append(p)
        else:single.append(p)
    pairs=[v for v in groups.values() if len(v.get('1',[]))==1 and len(v.get('2',[]))==1]
    if len(pairs)!=1 or len(groups)!=1:
        raise ValueError('Cannot uniquely identify R1/R2 FASTQ; found '+str([p.name for p in files]))
    pair=[need(pairs[0]['1'][0]),need(pairs[0]['2'][0])]
    return pair,[p for p in single if p.stat().st_size]

def merge_cna_files(inputs,out):
    # Preserve the user's working pandas merge, including original column names.
    import pandas as pd
    dfs=[]
    for sample,path in inputs:
        df=pd.read_csv(path,sep='\t')
        target_cols=[col for col in df.columns if 'Corrected_Copy_Number' in col]
        target_col=target_cols[0]
        df_sub=df[['chr','start','end',target_col]].copy()
        df_sub.set_index(['chr','start','end'],inplace=True)
        dfs.append(df_sub)
    merged_df=pd.concat(dfs,axis=1).reset_index()
    merged_df['chr']=pd.Categorical(merged_df['chr'],categories=[f'chr{i}' for i in range(1,23)],ordered=True)
    merged_df=merged_df.sort_values(by=['chr','start','end'])
    merged_df=merged_df.dropna(subset=['chr']).reset_index(drop=True)
    out=Path(out)
    if out.exists():raise FileExistsError(out)
    out.parent.mkdir(parents=True,exist_ok=True)
    merged_df.to_csv(out,index=False)

class Pipeline:
    def __init__(self):
        self.samples,self.runs,self.profiles=load_manifest()
        self.root=Path(os.environ.get('PROJECT_DIR','results/19_public_tissue_wgbs')).expanduser().resolve()
        self.threads=max(1,int(os.environ.get('THREADS','8')))
        self.bin=int(os.environ.get('BIN_SIZE','5000000'))
        self.ref=Path(os.environ.get('REF','references/hg38.fa')).resolve()
        self.cna_input=Path(os.environ.get('CNA_INPUT_DIR',str(self.root/'results/03_nanopore_cna/3.1_wig')))
        self.cna_output=Path(os.environ.get('CNA_OUTPUT_DIR',str(self.root/'results/03_nanopore_cna/3.3_ichorCNA')))
        self.cna_merge=Path(os.environ.get('CNA_MERGE_DIR',str(self.root/'results/03_nanopore_cna')))
        self.errors=[]
        self.code_hash=hashlib.sha256(b''.join((HERE/n).read_bytes() for n in ['_pipeline.py','_dmr_counts.sh'])).hexdigest()

    def run(self,args,log,stdout=None,cwd=None):
        log.write('$ '+shlex.join(list(map(str,args)))+'\n');log.flush()
        subprocess.run(list(map(str,args)),stdout=stdout or log,stderr=log,cwd=cwd,check=True)

    def pipe(self,left,right,log):
        log.write('$ '+shlex.join(list(map(str,left)))+' | '+shlex.join(list(map(str,right)))+'\n');log.flush()
        process=subprocess.Popen(list(map(str,left)),stdout=subprocess.PIPE,stderr=log)
        try:
            result=subprocess.run(list(map(str,right)),stdin=process.stdout,stdout=log,stderr=log)
        finally:
            process.stdout.close()
            code=process.wait()
        if result.returncode or code:raise RuntimeError(f'Pipeline failed: upstream={code}, downstream={result.returncode}')

    def raw(self,r):
        if r['sra_path']:return need(r['sra_path']).resolve()
        a=r['run_id'];base=Path(os.environ.get('SRA_DIR',''))
        for p in [base/(a+'.sra'),base/a/(a+'.sra')]:
            if p.is_file():return need(p).resolve()
        raise ValueError('SRA missing for '+a+'; set SRA_DIR or runs.tsv sra_path')

    def fq(self,r,trim=False):
        a=r['run_id']
        if not trim and r['fastq_r1']:
            files=[Path(r['fastq_r1'])]
            if r['layout']=='PAIRED':files.append(Path(r['fastq_r2']))
            return files
        base=self.root/('trimmed' if trim else 'fastq')/a
        return [base/(a+('.fastq.gz' if r['layout']=='SINGLE' else '_1.fastq.gz'))]+([] if r['layout']=='SINGLE' else [base/(a+'_2.fastq.gz')])

    def sample_bam(self,s):return self.root/'bam_processed'/(s['sample_id']+'.clean.q20.bam')

    def cna_file(self,s):
        return need(self.cna_output/(s['sample_id']+'.cna.seg'))

    def plan(self):
        print(f'{len(self.samples)} entries ({sum(s["material"]=="tissue" for s in self.samples)} tissue + {sum(s["material"]=="cell_line" for s in self.samples)} cell lines); {len(self.runs)} runs; {len({r["processing_group"] for r in self.runs})} processing groups')
        print('SINGLE runs:',sum(r['layout']=='SINGLE' for r in self.runs))
        print('PROJECT_DIR:',self.root)
        print('Library assumption: same sample+SRX -> same LB; NOT independently confirmed physical libraries.')
        print('No UMI extraction: status NOT_REPORTED, not proven absent. EM disabled.')
        for s in self.samples:
            rr=[r for r in self.runs if r['sample_id']==s['sample_id']]
            print(s['sample_id'],s['biosample'],s['group'],rr[0]['layout'],rr[0]['processing_group'],','.join(r['run_id'] for r in rr),sep='\t')

    def expected(self,step,item):
        a=item.get('run_id',item.get('sample_id'));base=self.root
        if step=='validate':return [base/'qc'/(a+'.validate.txt')]
        if step=='fastq':return self.fq(item)
        if step=='trim':return self.fq(item,True)+[base/'qc'/(a+'.fastp.json'),base/'qc'/(a+'.cutadapt.json')]
        if step=='bsmap':return [base/'aligned'/a/(a+'.rg.bam')]
        if step=='bam':return [self.sample_bam(item),base/'bam_processed'/(a+'.duplicate_metrics.txt'),Path(str(self.sample_bam(item))+'.bai')]
        if step=='mcall':return [Path(str(self.sample_bam(item))+'.G.bed')]
        if step=='dmr':return [base/'bam_processed/methylation/dmr_stats'/(a+'.tsv')]
        if step=='wig':return [self.cna_input/(a+'.wig')]
        if step=='cna':return [self.cna_output/(a+'.cna.seg')]
        raise ValueError(step)

    def inputs(self,step,item):
        a=item.get('run_id',item.get('sample_id'))
        if step=='validate':return self.fq(item) if item['fastq_r1'] else [self.raw(item)]
        if step=='fastq':return self.fq(item) if item['fastq_r1'] else [self.raw(item)]
        if step=='trim':return self.fq(item)
        if step=='bsmap':return self.fq(item,True)+[self.ref]
        if step=='bam':return [self.root/'aligned'/r['run_id']/(r['run_id']+'.rg.bam') for r in self.runs if r['sample_id']==a]
        if step=='mcall':return [self.sample_bam(item),Path(str(self.sample_bam(item))+'.bai'),self.ref]
        if step=='dmr':return [Path(str(self.sample_bam(item))+'.G.bed'),Path(os.environ['REGIONS'])]
        if step=='wig':return [self.sample_bam(item),Path(os.environ['GC_WIG']),Path(os.environ['MAP_WIG'])]
        if step=='cna':
            return [self.cna_input/(a+'.wig')]+[Path(os.environ[k]) for k in ['GC_WIG','MAP_WIG','CENTROMERE','ICHOR_SCRIPT']]
        return []

    def fingerprint(self,step,item):
        rr=[item] if 'run_id' in item else [r for r in self.runs if r['sample_id']==item['sample_id']]
        config={k:os.environ.get(k,'') for k in ['THREADS','BIN_SIZE','REGIONS','REF','GC_WIG','MAP_WIG','CENTROMERE','ICHOR_SCRIPT','CNA_INPUT_DIR','CNA_OUTPUT_DIR','CNA_MERGE_DIR','NORMAL_PANEL','ICHOR_NORMAL','ICHOR_PLOIDY','ICHOR_MAX_CN','SITE_INTERVAL_MODE','THRESHOLD']}
        payload=[self.code_hash,step,item,rr,[self.profiles[r['processing_group']] for r in rr],config,[sig(p) for p in self.inputs(step,item)]]
        return hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()

    def execute(self,step,item):
        a=item.get('run_id',item.get('sample_id'));outputs=self.expected(step,item)
        fingerprint=self.fingerprint(step,item);mark=self.root/'logs'/(a+'.'+step+'.ok.json')
        if mark.exists():
            record=json.loads(mark.read_text())
            if record['fingerprint']!=fingerprint:raise ValueError('Configuration/input changed. Use a NEW PROJECT_DIR: '+str(mark))
            if record['outputs']!=[sig(p) for p in outputs]:raise ValueError('Completed output changed: '+str(mark))
            print('SKIP',step,a,flush=True);return
        # External FASTQ validation does not overwrite user inputs.
        if step!='fastq' or not item.get('fastq_r1'):
            if any(p.exists() for p in outputs):raise ValueError('Unmarked output exists; use a new output directory: '+str(outputs[0]))
        with (self.root/'logs'/(a+'.'+step+'.log')).open('a') as log:
            tmp=Path(tempfile.mkdtemp(prefix=a+'.'+step+'.',dir=self.root/'tmp'))
            print('RUN',step,a,flush=True)
            try:self.body(step,item,tmp,log)
            except Exception:
                log.write('FAILED: work files retained at '+str(tmp)+'\n');log.flush()
                print('PRESERVED_WORKDIR',tmp,file=sys.stderr,flush=True)
                raise
            else:shutil.rmtree(tmp)
        for p in outputs:need(p)
        mark.write_text(json.dumps(dict(fingerprint=fingerprint,outputs=[sig(p) for p in outputs],finished=time.strftime('%Y-%m-%dT%H:%M:%S')),indent=2)+'\n')

    def publish(self,src,dst):
        dst=Path(dst);dst.parent.mkdir(parents=True,exist_ok=True)
        if dst.exists():raise FileExistsError(dst)
        shutil.move(str(src),str(dst))

    def body(self,step,item,tmp,log):
        a=item.get('run_id',item.get('sample_id'));outs=self.expected(step,item)
        if step=='validate':
            if item['fastq_r1']:
                self.check_fastq(self.fq(item),item['layout']);(tmp/'validated.txt').write_text('External FASTQ fully parsed; paired identifiers checked\n')
            else:
                with (tmp/'validated.txt').open('w') as h:self.run([Path(os.environ['SRA_BIN'])/'vdb-validate',self.raw(item)],log,h)
                # Some vdb-validate builds send their entire report to stderr.
                with (tmp/'validated.txt').open('a') as h:h.write('vdb-validate exit code 0\n')
            self.publish(tmp/'validated.txt',outs[0])
        elif step=='fastq':
            if item['fastq_r1']:
                self.check_fastq(self.fq(item),item['layout']);return
            need(self.root/'logs'/(a+'.validate.ok.json'))
            source=self.raw(item)
            # Older fasterq-dump builds reject accession directories. Pass the
            # local .sra file and an explicit output filename to avoid deriving
            # output paths from an absolute input accession.
            dump=tmp/'dump';dump.mkdir()
            self.run([(Path(os.environ['SRA_BIN'])/'fasterq-dump').resolve(),source,'--split-3','--threads',self.threads,'--temp',tmp/'scratch','--outfile',dump/(a+'.fastq')],log,cwd=dump)
            inventory=[{'name':p.name,'bytes':p.stat().st_size} for p in dump.iterdir() if p.is_file()]
            (self.root/'qc'/(a+'.fasterq_files.json')).write_text(json.dumps(inventory,indent=2)+'\n')
            raw,orphans=dump_outputs(dump,item['layout'])
            self.check_fastq(raw,item['layout'])
            ready=tmp/'ready';ready.mkdir()
            for src,dst in zip(raw,outs):
                target=ready/dst.name
                if src.name.endswith('.gz'):shutil.move(str(src),str(target))
                else:
                    plain=ready/dst.name[:-3];shutil.move(str(src),str(plain))
                    self.run(['pigz','-p',self.threads,plain],log)
            self.check_fastq([ready/p.name for p in outs],item['layout'])
            for i,orphan in enumerate(orphans,1):
                if not orphan.name.endswith('.gz'):
                    self.run(['pigz','-p',self.threads,orphan],log);orphan=Path(str(orphan)+'.gz')
                self.publish(orphan,self.root/'orphans'/(a+'.singleton'+str(i)+'.fastq.gz'))
            if orphans:(self.root/'qc'/(a+'.orphans.txt')).write_text('Singletons retained separately and excluded from PE analysis\n')
            # Publish the complete pair as one directory, not one mate at a time.
            self.publish(ready,outs[0].parent)
        elif step=='trim':
            p=self.profiles[item['processing_group']];paired=item['layout']=='PAIRED';inputs=self.fq(item)
            adapter=[tmp/('adapter_'+str(i)+'.fq.gz') for i in range(len(inputs))]
            cmd=['fastp','--in1',inputs[0],'--out1',adapter[0],'--thread',min(self.threads,int(os.environ.get('FASTP_THREADS','4'))),
                 '--disable_quality_filtering','--disable_length_filtering','--disable_trim_poly_g','--dont_eval_duplication',
                 '--json',tmp/'fastp.json','--html',tmp/'fastp.html']
            if paired:cmd+=['--in2',inputs[1],'--out2',adapter[1],'--detect_adapter_for_pe']
            if p['adapter_r1']!='AUTO':cmd+=['--adapter_sequence',p['adapter_r1']]
            if paired and p['adapter_r2']!='AUTO':cmd+=['--adapter_sequence_r2',p['adapter_r2']]
            self.run(cmd,log)
            trimmed=[tmp/('final_'+str(i)+'.fq.gz') for i in range(len(inputs))]
            self.run(cutadapt_args(p,paired,adapter,trimmed,tmp/'cutadapt.json',self.threads),log)
            self.check_fastq(trimmed,item['layout'])
            for src,dst in zip(trimmed,self.fq(item,True)):self.publish(src,dst)
            for name in ['fastp.json','fastp.html','cutadapt.json']:self.publish(tmp/name,self.root/'qc'/(a+'.'+name))
        elif step=='bsmap':
            p=self.profiles[item['processing_group']];ff=self.fq(item,True)
            cmd=['bsmap','-a',ff[0]]
            if len(ff)==2:cmd+=['-b',ff[1]]
            cmd+=['-d',self.ref,'-z','33','-S','123','-n',p['bsmap_n'],'-r','0','-u','-p',self.threads]
            self.pipe(cmd,['samtools','sort','-@',max(0,self.threads-1),'-m','512M','-o',tmp/'sorted.bam','-'],log)
            # Preserve one RG per run, shared LB across assumed same sample+SRX.
            self.run(['picard','-Xmx'+os.environ.get('PICARD_HEAP','16g'),'AddOrReplaceReadGroups',
                'I='+str(tmp/'sorted.bam'),'O='+str(tmp/'rg.bam'),'RGID='+a,'RGLB='+item['library_id'],
                'RGSM='+item['sample_id'],'RGPL='+item['platform'],'RGPU='+a,'SORT_ORDER=coordinate','VALIDATION_STRINGENCY=LENIENT'],log)
            self.run(['samtools','quickcheck','-v',tmp/'rg.bam'],log);self.publish(tmp/'rg.bam',outs[0])
        elif step=='bam':
            inputs=self.inputs(step,item);merged=tmp/'merged.bam'
            if len(inputs)>1:self.run(['samtools','merge','-@',max(0,self.threads-1),'-o',merged,*inputs],log)
            else:merged=inputs[0]
            self.run(['picard','-Xmx'+os.environ.get('PICARD_HEAP','16g'),'MarkDuplicates','I='+str(merged),
                'O='+str(tmp/'marked.bam'),'M='+str(tmp/'metrics.txt'),'REMOVE_DUPLICATES=false','READ_NAME_REGEX=null',
                'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=512','TMP_DIR='+str(tmp)],log)
            self.run(['samtools','view','-b','-F','3844','-q','20','-o',tmp/'clean.bam',tmp/'marked.bam'],log)
            self.run(['samtools','quickcheck','-v',tmp/'clean.bam'],log)
            self.run(['samtools','index',tmp/'clean.bam'],log)
            with (tmp/'flagstat.txt').open('w') as h:self.run(['samtools','flagstat',tmp/'clean.bam'],log,h)
            self.publish(tmp/'clean.bam',outs[0]);self.publish(tmp/'clean.bam.bai',Path(str(outs[0])+'.bai'))
            self.publish(tmp/'metrics.txt',outs[1]);self.publish(tmp/'flagstat.txt',self.root/'qc'/(a+'.flagstat.txt'))
        elif step=='mcall':
            rr=next(r for r in self.runs if r['sample_id']==a);p=self.profiles[rr['processing_group']]
            if not Path(str(self.ref)+'.fai').exists():raise ValueError('Create REF.fai once using samtools faidx before running')
            # Private BAM path prevents mcall builds that ignore outputDir writing beside production BAM.
            local=tmp/self.sample_bam(item).name;local.symlink_to(self.sample_bam(item).resolve())
            Path(str(local)+'.bai').symlink_to(need(str(self.sample_bam(item))+'.bai').resolve())
            cmd=['mcall','-m',local,'-r',self.ref,'-g','hg38','--outputDir',tmp,
                 '--qualityScoreBase','33','--cytosineMinScore','20','--nextBaseMinScore','3','--reportCpX','G',
                 '--trimRRBSEndRepairSeq','0','--trimWGBSEndRepairPE1Seq',p['mcall_pe1'] if rr['layout']=='PAIRED' else '0',
                 '--trimWGBSEndRepairPE2Seq',p['mcall_pe2'] if rr['layout']=='PAIRED' else '0',
                 '--processPEOverlapSeq','1','--fullMode','0','-p',os.environ.get('MCALL_THREADS',str(self.threads))]
            self.run(cmd,log)
            source=need(Path(str(local)+'.G.bed'));self.publish(source,outs[0])
            for f in tmp.glob('*stat*'):
                if f.is_file():self.publish(f,self.root/'qc'/(a+'.'+f.name))
        elif step=='dmr':
            env=os.environ.copy();env['SAMPLE_ONLY']=a
            subprocess.run(['bash',str(HERE/'_dmr_counts.sh')],stdout=log,stderr=log,env=env,check=True)
        elif step=='wig':
            self.check_wigs()
            with (tmp/'counts.wig').open('w') as h:self.run(['readCounter','--window',self.bin,'--quality','20','--chromosome',','.join('chr'+str(i) for i in range(1,23)),self.sample_bam(item)],log,h)
            if 'fixedStep' not in (tmp/'counts.wig').read_text():raise ValueError('Invalid readCounter WIG')
            self.publish(tmp/'counts.wig',outs[0])
        elif step=='cna':
            self.check_wigs()
            cmd=['Rscript',os.environ['ICHOR_SCRIPT'],'--id',a,'--WIG',self.cna_input/(a+'.wig'),
                 '--gcWig',os.environ['GC_WIG'],'--mapWig',os.environ['MAP_WIG'],'--centromere',os.environ['CENTROMERE'],
                 '--ploidy',os.environ.get('ICHOR_PLOIDY','c(2)'),'--normal',os.environ.get('ICHOR_NORMAL','c(0.4,0.5,0.6,0.7,0.8,0.9)'),'--maxCN',os.environ.get('ICHOR_MAX_CN','3'),
                 '--includeHOMD','FALSE','--chrs','c('+','.join(map(str,range(1,23)))+')','--chrTrain','c('+','.join(map(str,range(1,23)))+')','--genomeStyle','UCSC',
                 '--estimateNormal','TRUE','--estimatePloidy','TRUE','--estimateScPrevalence','FALSE','--altFracThreshold','0.01',
                 '--rmCentromereFlankLength','1000000','--outDir',tmp/'result']
            self.run(cmd,log);candidates=list((tmp/'result').rglob(a+'.cna.seg'))
            if len(candidates)!=1:raise ValueError('Missing/ambiguous CNA result')
            if candidates[0].parent!=tmp/'result':shutil.copy(candidates[0],tmp/'result'/(a+'.cna.seg'))
            self.cna_output.mkdir(parents=True,exist_ok=True)
            for child in (tmp/'result').iterdir():self.publish(child,self.cna_output/child.name)

    @staticmethod
    def check_fastq(paths,layout):
        def records(path):
            opener=gzip.open if str(path).endswith('.gz') else open
            with opener(need(path),'rt') as h:
                while True:
                    name=h.readline()
                    if not name:break
                    seq=h.readline().rstrip();plus=h.readline();quality=h.readline().rstrip()
                    if not name.startswith('@') or not plus.startswith('+') or len(seq)!=len(quality):raise ValueError('Invalid FASTQ '+str(path))
                    yield re.sub(r'/[12]$','',name.split()[0])
        count=0
        if layout=='SINGLE':
            for _ in records(paths[0]):count+=1
        else:
            from itertools import zip_longest
            for r1,r2 in zip_longest(records(paths[0]),records(paths[1])):
                if r1 is None or r2 is None or r1!=r2:raise ValueError('R1/R2 identifiers or counts differ')
                count+=1
        if count==0:raise ValueError('No reads retained')

    def check_wigs(self):
        reference_headers=[]
        for k in ['GC_WIG','MAP_WIG']:
            headers=[]
            with need(os.environ[k]).open() as h:
                for line in h:
                    if line.startswith('fixedStep'):
                        attrs=dict(x.split('=',1) for x in line.split()[1:])
                        if int(attrs.get('step','0'))!=self.bin or int(attrs.get('span',str(self.bin)))!=self.bin:raise ValueError(k+' does not match BIN_SIZE')
                        if not attrs.get('chrom','').startswith('chr'):raise ValueError(k+' must use UCSC chr names')
                        headers.append(line.strip())
            if not headers:raise ValueError('No fixedStep headers in '+k)
            reference_headers.append(headers)
        if reference_headers[0]!=reference_headers[1]:raise ValueError('GC/map WIG chromosome grids differ')

    def merge_guard(self,step,inputs,dest,action):
        signature=hashlib.sha256(json.dumps([self.code_hash,step,[sig(p) for p in inputs],
            [s['sample_id'] for s in self.samples],os.environ.get('THRESHOLD','1')],sort_keys=True).encode()).hexdigest()
        mark=self.root/'logs'/(step+'.ok.json')
        if mark.exists():
            record=json.loads(mark.read_text())
            if record['fingerprint']!=signature or record['output']!=sig(dest):raise ValueError('Merged output/input changed; use a new PROJECT_DIR')
            print('SKIP',step);return
        action();mark.write_text(json.dumps(dict(fingerprint=signature,output=sig(dest)))+'\n')

    def merge_dmr(self):
        records=[];base=self.root/'bam_processed/methylation/dmr_stats';coords=None
        for s in self.samples:
            rows=table(need(base/(s['sample_id']+'.tsv')));values={}
            for r in rows:
                key=(r['chrom'],int(r['start']),int(r['end']))
                if key in values:raise ValueError('Duplicate DMR coordinate')
                values[key]=r['percent_m'] if int(r['count_valid_m'])>=int(os.environ.get('THRESHOLD','1')) else 'NA'
            if coords is not None and set(values)!=coords:raise ValueError('Inconsistent DMR regions')
            coords=set(values);records.append((s['sample_id'],values))
        dest=base/'merge'/('threshold_'+os.environ.get('THRESHOLD','1')+'_wgbs_methylationfile.tsv.gz')
        if dest.exists():raise FileExistsError(dest)
        dest.parent.mkdir(parents=True,exist_ok=True)
        with gzip.open(dest,'wt',newline='') as h:
            w=csv.writer(h,delimiter='\t');w.writerow(['chr','start','end']+[s for s,v in records])
            for k in sorted(coords):w.writerow([*k]+[v[k] for s,v in records])

    def qc(self):
        out=[]
        for r in self.runs:
            a=r['run_id'];p=self.profiles[r['processing_group']];warning=[];fraction='NA'
            f=self.root/'qc'/(a+'.fastp.json')
            if f.exists():
                j=json.loads(f.read_text());before=j.get('summary',{}).get('before_filtering',{}).get('total_reads',0)
                c=self.root/'qc'/(a+'.cutadapt.json')
                if c.exists():
                    counts=json.loads(c.read_text()).get('read_counts',{});n=counts.get('input',0);kept=counts.get('output',0)
                    if n:fraction=round(kept/n,4)
                    if n and kept/n<0.5:warning.append('LOW_RETENTION_REVIEW')
                ad=j.get('adapter_cutting',{})
                if not ad.get('adapter_trimmed_reads',0):warning.append('NO_ADAPTER_TRIMMING_REPORTED_NOT_PROOF_OF_ABSENCE')
            else:warning.append('TRIM_NOT_COMPLETED')
            if p['bsmap_n']=='1' and r['processing_group']!='PG13':warning.append('DIRECTION_UNRESOLVED_FOUR_STRAND_SEARCH')
            warning+=['LIBRARY_ID_ASSUMED_FROM_SRX','UMI_NOT_REPORTED','FIXED_CLIPPING_NOT_MBIAS_VALIDATED']
            out.append(dict(run_id=a,sample_id=r['sample_id'],layout=r['layout'],processing_group=r['processing_group'],retained_fraction=fraction,warnings=';'.join(warning)))
        with (self.root/'qc/run_qc_summary.tsv').open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
        summary=[]
        for s in self.samples:
            a=s['sample_id'];f=self.root/'bam_processed/methylation/dmr_stats'/(a+'.tsv')
            total=missing='NA';frac='NA'
            if f.exists():
                rows=table(f);total=len(rows);missing=sum(r['percent_m']=='NA' for r in rows)
                frac=round(missing/total,4) if total else 'NA'
            duplicate='NA';metrics=self.root/'bam_processed'/(a+'.duplicate_metrics.txt')
            if metrics.exists():
                lines=metrics.read_text().splitlines()
                for i,line in enumerate(lines):
                    if line.startswith('LIBRARY\t') and 'PERCENT_DUPLICATION' in line:
                        col=line.split('\t').index('PERCENT_DUPLICATION');values=[]
                        for row in lines[i+1:]:
                            if not row or row.startswith('#'):break
                            try:values.append(float(row.split('\t')[col]))
                            except (ValueError,IndexError):break
                        if values:duplicate=max(values)
                        break
            summary.append(dict(sample_id=a,group=s['group'],material=s['material'],pooling=s['pooling'],region_count=total,
                no_coverage_regions=missing,no_coverage_fraction=frac,max_library_duplicate_fraction=duplicate))
        with (self.root/'qc/sample_qc_summary.tsv').open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=list(summary[0]),delimiter='\t');w.writeheader();w.writerows(summary)

    def main(self,step):
        if step=='plan':self.plan();return
        if '/path/to/' in str(self.root):raise ValueError('Set PROJECT_DIR in config.sh first')
        self.root.mkdir(parents=True,exist_ok=True)
        for d in ['logs','qc','tmp']:(self.root/d).mkdir(exist_ok=True)
        if step=='check':
            for tool in ['fastp','cutadapt','bsmap','samtools','picard','mcall','pigz','bedtools']:
                if not shutil.which(tool):raise ValueError('Required executable missing: '+tool)
            need(self.ref);need(str(self.ref)+'.fai');need(os.environ['REGIONS'])
            contigs={l.split('\t')[0]:int(l.split('\t')[1]) for l in Path(str(self.ref)+'.fai').read_text().splitlines()}
            if not all('chr'+str(i) in contigs for i in range(1,23)):raise ValueError('Reference must contain UCSC chr1..chr22')
            count=0
            for l in Path(os.environ['REGIONS']).read_text().splitlines():
                if not l.strip() or l.startswith(('#','track','browser')):continue
                f=l.split()
                if len(f)<3 or f[0] not in contigs or not (0<=int(f[1])<int(f[2])<=contigs[f[0]]):raise ValueError('BED/reference coordinate mismatch: '+l)
                count+=1
            if not count:raise ValueError('Empty region BED')
            for r in self.runs:
                if r['fastq_r1']:
                    for f in self.fq(r):need(f)
                else:self.raw(r)
            if any(not r['fastq_r1'] for r in self.runs):
                for tool in ['vdb-validate','fasterq-dump']:need(Path(os.environ['SRA_BIN'])/tool)
            versions={}
            for tool,args in [('fastp',['--version']),('cutadapt',['--version']),('samtools',['--version']),('picard',['MarkDuplicates','--version']),('mcall',['--help']),('bsmap',['-h'])]:
                try:
                    result=subprocess.run([tool,*args],capture_output=True,text=True,timeout=20)
                    versions[tool]=(result.stdout+'\n'+result.stderr)[:30000]
                except Exception as e:versions[tool]='Version/help capture failed: '+str(e)
            (self.root/'qc/tool_versions_and_help.json').write_text(json.dumps(versions,indent=2)+'\n')
            for name in ['samples.tsv','runs.tsv','profiles.tsv']:
                shutil.copy(METADATA_DIR/name,self.root/'qc'/('snapshot_'+name))
            print('Paths/executables present. Versions and biological QC still require actual execution.');return
        if step=='qc':self.qc();return
        if step=='em':raise ValueError('EM disabled for tissue WGBS')
        if step=='merge_dmr':
            base=self.root/'bam_processed/methylation/dmr_stats'
            dest=base/'merge'/('threshold_'+os.environ.get('THRESHOLD','1')+'_wgbs_methylationfile.tsv.gz')
            self.merge_guard(step,[base/(s['sample_id']+'.tsv') for s in self.samples],dest,self.merge_dmr);return
        if step=='merge_cna':
            inputs=[(s['sample_id'],self.cna_file(s)) for s in self.samples];dest=self.cna_merge/'discovery_Corrected_Copy_Number.csv'
            self.merge_guard(step,[p for s,p in inputs],dest,lambda:merge_cna_files(inputs,dest));return
        if step not in RUN_STEPS|{'bam','mcall','dmr','wig','cna'}:raise ValueError('Unknown step '+step)
        for item in self.runs if step in RUN_STEPS else self.samples:
            try:self.execute(step,item)
            except Exception as e:
                a=item.get('run_id',item.get('sample_id'));self.errors.append((step,a,str(e)))
                print('FAILED',step,a,str(e),file=sys.stderr,flush=True)
                if os.environ.get('FAIL_FAST')=='1':break
        with (self.root/'logs'/(step+'.failures.tsv')).open('w',newline='') as h:
            w=csv.writer(h,delimiter='\t');w.writerow(['step','id','error']);w.writerows(self.errors)
        if self.errors:raise RuntimeError(f'{len(self.errors)} failed items; other independent items were attempted. See logs/{step}.failures.tsv')

if __name__=='__main__':
    try:Pipeline().main(sys.argv[1] if len(sys.argv)>1 else 'plan')
    except Exception as e:print('ERROR:',e,file=sys.stderr);sys.exit(1)
