"""Run copied ROC/bootstrap scripts on invented data, never manuscript records."""
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    output = ROOT / "examples/output"
    output.mkdir(exist_ok=True)
    metadata = str(ROOT / "examples/synthetic_metadata.csv")
    predictions = str(ROOT / "examples/synthetic_predictions.csv")
    roc = load("demo_roc", "scripts/11_roc_plot/11.01_roc_plot.py")
    roc.GROUP_FILE = metadata
    roc.FRAG_CSV = roc.CNA_CSV = roc.METH_CSV = roc.MULTI_CSV = predictions
    roc.OUTPUT_ROC = str(output / "synthetic_roc.pdf")
    roc.mpl.rcParams["font.family"] = "DejaVu Sans"
    roc.main()
    roc.plt.gcf().suptitle("Synthetic demonstration — not participant data")
    roc.plt.savefig(roc.OUTPUT_ROC)
    ci = load("demo_ci", "scripts/15_bootstrap_ci/15.01_bootstrap_auc_ci.py")
    ci.GROUP_FILE = metadata
    for info in ci.MODEL_FILES.values():
        info["file"] = predictions
    ci.OUTPUT_DIR = str(output)
    ci.OUTPUT_CSV = str(output / "synthetic_metrics.csv")
    ci.N_BOOTSTRAP = 200  # Demo speed only; original script remains unchanged.
    ci.main()
    print("Synthetic example completed; outputs are not manuscript findings.")


if __name__ == "__main__":
    main()
