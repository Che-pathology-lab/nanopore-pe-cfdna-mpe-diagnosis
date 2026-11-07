#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/19.00_settings.sh"
INPUT_DIR="$PROJECT_DIR/bam_processed"
OUTPUT_DIR="$INPUT_DIR/methylation/dmr_stats"
EXPECTED_REGIONS="${EXPECTED_REGIONS:-3995}"


command -v bedtools >/dev/null || {
    echo "错误：没有找到bedtools，请激活包含bedtools的环境" >&2
    exit 1
}

[[ -s "${REGIONS}" ]] || {
    echo "错误：找不到DMR文件：${REGIONS}" >&2
    exit 1
}

[[ -d "${INPUT_DIR}" ]] || {
    echo "错误：输入目录不存在：${INPUT_DIR}" >&2
    exit 1
}

case "${SITE_INTERVAL_MODE}" in
    original|anchor) ;;
    *)
        echo "错误：SITE_INTERVAL_MODE必须是original或anchor" >&2
        exit 1
        ;;
esac

mkdir -p "${OUTPUT_DIR}"

# 本脚本由run_dmr.sh按列表逐样本调用。
BED_FILES=("$INPUT_DIR/$SAMPLE_ONLY.clean.q20.bam.G.bed")
[[ -s "$INPUT_DIR/$SAMPLE_ONLY.clean.q20.bam.G.bed" ]] || {
    echo "缺少mcall汇总：$SAMPLE_ONLY" >&2
    exit 1
}
export LC_ALL=C

WORK_DIR=$(mktemp -d "${OUTPUT_DIR}/.mcall_dmr_tmp.XXXXXX")
trap 'rm -r -- "${WORK_DIR}"' EXIT

SORTED_REGIONS="${WORK_DIR}/regions.sorted.bed"

echo "整理DMR区域……"

# 去除Windows换行符，转换为BED4并排序
sed 's/\r$//' "${REGIONS}" |
awk 'BEGIN {
    OFS="\t"
}
/^[[:space:]]*$/ {next}
/^[[:space:]]*#/ {next}
/^track[[:space:]]/ {next}
/^browser[[:space:]]/ {next}

NF<3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ {
    print "错误：非BED记录：" $0 > "/dev/stderr"
    exit 1
}

$3<=$2 {
    print "错误：区域坐标异常：" $0 > "/dev/stderr"
    exit 1
}

{
    name=(NF>=4 && $4!="" ? $4 : $1 ":" $2 "-" $3)
    print $1,$2,$3,name
}' |
bedtools sort -i - > "${SORTED_REGIONS}"

REGION_COUNT=$(wc -l < "${SORTED_REGIONS}")

echo "有效DMR数量：${REGION_COUNT}"

if [[ ${REGION_COUNT} -eq 0 ]]; then
    echo "错误：没有得到有效DMR区域" >&2
    exit 1
fi

if [[ ${REGION_COUNT} -ne ${EXPECTED_REGIONS} ]]; then
    echo "警告：预期${EXPECTED_REGIONS}个DMR，实际${REGION_COUNT}个" >&2
fi

echo "发现 ${#BED_FILES[@]} 个mcall汇总文件"
echo "CpG区间模式：${SITE_INTERVAL_MODE}"

for BED_FILE in "${BED_FILES[@]}"; do
    SAMPLE=$(basename "${BED_FILE}" .clean.q20.bam.G.bed)

    SITE_BED="${WORK_DIR}/${SAMPLE}.sites.bed"
    MAPPED="${WORK_DIR}/${SAMPLE}.mapped.bed"
    TEMP_OUTPUT="${WORK_DIR}/${SAMPLE}.tsv"
    OUTPUT="${OUTPUT_DIR}/${SAMPLE}.tsv"

    echo "正在处理：${SAMPLE}"

    if [[ -e "${OUTPUT}" ]]; then
        echo "错误：结果已存在，为避免覆盖停止：${OUTPUT}" >&2
        exit 1
    fi

    # mcall输入：
    # chrom start end ratio totalC methC strand next ...
    #
    # 转换为：
    # chrom start end methC totalC
    sed 's/\r$//' "${BED_FILE}" |
    awk -v mode="${SITE_INTERVAL_MODE}" '
    BEGIN {
        OFS="\t"
    }
    /^[[:space:]]*$/ {next}
    /^[[:space:]]*#/ {next}

    NF<8 {
        print "错误：mcall记录列数不足：" $0 > "/dev/stderr"
        exit 1
    }

    $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3<=$2 {
        print "错误：mcall坐标异常：" $0 > "/dev/stderr"
        exit 1
    }

    $5 !~ /^[0-9]+$/ || $6 !~ /^[0-9]+$/ || $6>$5 {
        print "错误：mcall计数异常：" $0 > "/dev/stderr"
        exit 1
    }

    # 仅使用CpG，忽略无有效覆盖记录
    $8!="G" || $5==0 {next}

    {
        end=(mode=="anchor" ? $2+1 : $3)
        print $1,$2,end,$6,$5
    }' |
    bedtools sort -i - > "${SITE_BED}"

    if [[ ! -s "${SITE_BED}" ]]; then
        echo "错误：没有得到有效CpG记录：${BED_FILE}" >&2
        exit 1
    fi

    # 每个DMR汇总甲基化计数和有效覆盖计数
    bedtools map \
        -a "${SORTED_REGIONS}" \
        -b "${SITE_BED}" \
        -c 4,5 \
        -o sum,sum \
        -null 0 \
        > "${MAPPED}"

    # 保持与之前脚本相同的8列输出结构
    awk 'BEGIN {
        OFS="\t"
        print "chrom","start","end","name","strand", \
              "count_m","count_valid_m","percent_m"
    }
    {
        count_m=$5+0
        count_valid=$6+0

        percent=(count_valid>0 \
            ? sprintf("%.6f",100*count_m/count_valid) \
            : "NA")

        print $1,$2,$3,$4,".", \
              sprintf("%.0f",count_m), \
              sprintf("%.0f",count_valid),percent
    }' "${MAPPED}" > "${TEMP_OUTPUT}"

    OUTPUT_COUNT=$(awk 'END {print NR-1}' "${TEMP_OUTPUT}")

    if [[ ${OUTPUT_COUNT} -ne ${REGION_COUNT} ]]; then
        echo "错误：输出区域数不正确：${SAMPLE}" >&2
        exit 1
    fi

    mv "${TEMP_OUTPUT}" "${OUTPUT}"

    # 保存所用CpG区间定义
    printf "SITE_INTERVAL_MODE=%s\n" "${SITE_INTERVAL_MODE}" \
        > "${OUTPUT}.meta.txt"

    echo "完成：${OUTPUT}"
done

echo "所有样本处理完成"
echo "结果目录：${OUTPUT_DIR}"
