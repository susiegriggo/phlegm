#!/usr/bin/env bash
# Set up PHLEGM inside the currently active conda environment.
#
# First create and activate an environment with the dependencies, e.g.:
#   conda env create -f environment.yml        # add '-n NAME' or '-p PATH' to choose
#   conda activate phlegm
# Then:
#   bash install.sh          # Linux: GPU (CUDA 12) JAX, which also works without a GPU
#                            # macOS: CPU JAX
#   bash install.sh --cpu    # Linux: skip GPU JAX (saves ~3 GB if you will never fold)
#
# AlphaFold weights (~4 GB) are downloaded on the first 'phlegm fold', or ahead
# of time with 'phlegm download-weights' (do this on HPC login nodes).
#
# Afterwards:
#   phlegm run -i proteins.fasta -o phlegm_out
set -euo pipefail

GPU_JAX=1
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() { sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cpu) GPU_JAX=0; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

# ---- check we are inside a suitable environment ------------------------------
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "Error: no conda environment is active. Create and activate one first:" >&2
    echo "  conda env create -f environment.yml && conda activate phlegm" >&2
    exit 1
fi
if [[ "${CONDA_DEFAULT_ENV:-}" == "base" ]]; then
    echo "Error: refusing to install into the conda 'base' environment." >&2
    echo "  conda env create -f environment.yml && conda activate phlegm" >&2
    exit 1
fi
if ! command -v colabfold_batch >/dev/null 2>&1; then
    echo "Error: ColabFold was not found in the active environment ($CONDA_PREFIX)." >&2
    echo "Add PHLEGM's dependencies to it with:" >&2
    echo "  conda env update -f $SCRIPT_DIR/environment.yml" >&2
    exit 1
fi
echo "==> Installing PHLEGM into $CONDA_PREFIX"

PY="$CONDA_PREFIX/bin/python"

echo "==> Installing the 'phlegm' command"
"$PY" -m pip install --no-deps --no-build-isolation -e "$SCRIPT_DIR"

# ---- GPU support -------------------------------------------------------------
# ColabFold's conda package ships CPU-only JAX, which is far too slow for
# multimers. On Linux we always install CUDA JAX: it falls back to the CPU when
# no GPU is present, so one environment works on login, CPU and GPU nodes.
# Whether a GPU is actually available is checked by 'phlegm fold' at run time.
if [[ "$(uname -s)" == "Linux" && "$GPU_JAX" -eq 1 ]]; then
    echo "==> Installing CUDA-enabled JAX"
    # The bioconda colabfold package cannot be solved together with conda's
    # CUDA jaxlib builds, so (as LocalColabFold does) install CUDA JAX from pip.
    # cuDNN is pinned to 8.9 to match the jaxlib 0.4.25 cudnn89 build.
    "$PY" -m pip install --upgrade \
        "jax[cuda12_pip]==0.4.25" "numpy<2" "nvidia-cudnn-cu12>=8.9,<9" \
        -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
    if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
        "$PY" -c "import jax; d = jax.devices(); print('JAX devices:', d); assert d[0].platform == 'gpu', 'JAX cannot see the GPU'"
    fi
elif [[ "$(uname -s)" == "Linux" ]]; then
    echo "==> Skipping GPU JAX (--cpu). 'phlegm fold' will need --allow-cpu and will be slow."
else
    echo "==> $(uname -s): using CPU JAX (GPU JAX requires Linux with an NVIDIA GPU)."
    echo "    Fine for 'phlegm msa' and 'phlegm predict'; 'phlegm fold' needs --allow-cpu and is slow."
fi

echo "==> Checking installation"
colabfold_batch --help >/dev/null
"$CONDA_PREFIX/bin/phlegm" --help >/dev/null

cat <<EOF

PHLEGM installed into $CONDA_PREFIX.

  phlegm run -i proteins.fasta -o phlegm_out

See 'phlegm --help' for the separate msa / fold / predict commands.
AlphaFold weights (~4 GB) download automatically on the first 'phlegm fold'.
On HPC, fetch them now on the login node: phlegm download-weights
By default MSAs are generated with the public ColabFold MMseqs2 server.
To use a local ColabFold database instead, add: --db /path/to/colabfold_db
EOF
