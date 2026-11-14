#!/usr/bin/env bash
# One-step installer for PHLEGM.
#
#   bash install.sh                 # creates the conda env "phlegm"
#   bash install.sh -n myenv        # custom environment name
#   bash install.sh --gpu           # install GPU (CUDA 12) JAX even if no GPU is
#                                   #   visible, e.g. on an HPC login node
#   bash install.sh --cpu           # do not install GPU JAX
#   bash install.sh --no-weights    # skip downloading AlphaFold weights
#
# By default GPU JAX is installed only if an NVIDIA GPU is visible here.
#
# Afterwards:
#   conda activate phlegm
#   phlegm run -i proteins.fasta -o phlegm_out
set -euo pipefail

ENV_NAME="phlegm"
GPU_MODE=auto
DOWNLOAD_WEIGHTS=1
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() { sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--name) ENV_NAME="$2"; shift 2 ;;
        --gpu) GPU_MODE=gpu; shift ;;
        --cpu) GPU_MODE=cpu; shift ;;
        --no-weights) DOWNLOAD_WEIGHTS=0; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

# Prefer mamba (faster solver) if available.
if command -v mamba >/dev/null 2>&1; then
    CONDA=mamba
elif command -v conda >/dev/null 2>&1; then
    CONDA=conda
else
    echo "Error: conda or mamba is required. Install Miniforge: https://github.com/conda-forge/miniforge" >&2
    exit 1
fi

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "Error: a conda environment named '$ENV_NAME' already exists." >&2
    echo "Remove it with 'conda env remove -n $ENV_NAME' or choose another name with -n." >&2
    exit 1
fi

echo "==> Creating conda environment '$ENV_NAME' with $CONDA"
"$CONDA" env create -n "$ENV_NAME" -f "$SCRIPT_DIR/environment.yml"

echo "==> Installing the 'phlegm' command"
conda run -n "$ENV_NAME" python -m pip install --no-deps --no-build-isolation -e "$SCRIPT_DIR"

# GPU support: ColabFold's conda package ships CPU-only JAX, which is far too
# slow for multimers. 'phlegm fold' needs GPU JAX; 'phlegm msa' does not.
HAS_GPU=0
if [[ "$(uname -s)" == "Linux" ]] && command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
    HAS_GPU=1
fi
if [[ "$GPU_MODE" == "gpu" || ( "$GPU_MODE" == "auto" && "$HAS_GPU" -eq 1 ) ]]; then
    if [[ "$(uname -s)" != "Linux" ]]; then
        echo "Error: GPU JAX is only available on Linux." >&2
        exit 1
    fi
    echo "==> Installing CUDA-enabled JAX"
    # The bioconda colabfold package cannot be solved together with conda's
    # CUDA jaxlib builds, so (as LocalColabFold does) install CUDA JAX from pip.
    # cuDNN is pinned to 8.9 to match the jaxlib 0.4.25 cudnn89 build.
    conda run -n "$ENV_NAME" python -m pip install --upgrade \
        "jax[cuda12_pip]==0.4.25" "numpy<2" "nvidia-cudnn-cu12>=8.9,<9" \
        -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
    if [[ "$HAS_GPU" -eq 1 ]]; then
        conda run -n "$ENV_NAME" python -c \
            "import jax; d = jax.devices(); print('JAX devices:', d); assert d[0].platform == 'gpu', 'JAX cannot see the GPU'"
    else
        echo "    No GPU visible here; 'phlegm fold' will check for the GPU when it runs on a GPU node."
    fi
else
    echo "==> Using CPU-only JAX (no NVIDIA GPU visible, or --cpu given)."
    echo "    That is fine for 'phlegm msa' and 'phlegm predict'. For 'phlegm fold' on a"
    echo "    cluster, re-install with --gpu (GPU nodes usually aren't visible from login nodes)."
fi

if [[ "$DOWNLOAD_WEIGHTS" -eq 1 ]]; then
    echo "==> Downloading AlphaFold-Multimer v3 weights (cached for later runs)"
    conda run -n "$ENV_NAME" python -m colabfold.download alphafold2_multimer_v3
fi

echo "==> Checking installation"
conda run -n "$ENV_NAME" colabfold_batch --help >/dev/null
conda run -n "$ENV_NAME" phlegm --help >/dev/null

cat <<EOF

PHLEGM installed successfully.

  conda activate $ENV_NAME
  phlegm run -i proteins.fasta -o phlegm_out

See 'phlegm --help' for the separate msa / fold / predict commands.

By default MSAs are generated with the public ColabFold MMseqs2 server.
To use a local ColabFold database instead, add: --db /path/to/colabfold_db
EOF
