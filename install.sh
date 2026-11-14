#!/usr/bin/env bash
# One-step installer for PHLEGM.
#
#   bash install.sh                 # creates the conda env "phlegm"
#   bash install.sh -n myenv        # custom environment name
#   bash install.sh --cpu           # do not install GPU (CUDA) JAX
#   bash install.sh --no-weights    # skip downloading AlphaFold weights
#
# Afterwards:
#   conda activate phlegm
#   python phlegm.py -i proteins.fasta -o phlegm_out
set -euo pipefail

ENV_NAME="phlegm"
FORCE_CPU=0
DOWNLOAD_WEIGHTS=1
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() { sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--name) ENV_NAME="$2"; shift 2 ;;
        --cpu) FORCE_CPU=1; shift ;;
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

# GPU support: ColabFold runs on CPU by default, which is very slow for
# multimers. If an NVIDIA GPU is present, swap in the CUDA build of jaxlib.
if [[ "$FORCE_CPU" -eq 0 && "$(uname -s)" == "Linux" ]] && command -v nvidia-smi >/dev/null 2>&1; then
    echo "==> NVIDIA GPU detected; installing CUDA-enabled JAX"
    # The bioconda colabfold package cannot be solved together with conda's
    # CUDA jaxlib builds, so (as LocalColabFold does) install CUDA JAX from pip.
    # cuDNN is pinned to 8.9 to match the jaxlib 0.4.25 cudnn89 build.
    conda run -n "$ENV_NAME" python -m pip install --upgrade \
        "jax[cuda12_pip]==0.4.25" "numpy<2" "nvidia-cudnn-cu12>=8.9,<9" \
        -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
    conda run -n "$ENV_NAME" python -c \
        "import jax; d = jax.devices(); print('JAX devices:', d); assert d[0].platform == 'gpu', 'JAX cannot see the GPU'"
else
    echo "==> No NVIDIA GPU detected (or --cpu given); using CPU JAX."
    echo "    Structure prediction will be slow. Re-run on a GPU machine for practical runtimes."
fi

if [[ "$DOWNLOAD_WEIGHTS" -eq 1 ]]; then
    echo "==> Downloading AlphaFold-Multimer v3 weights (cached for later runs)"
    conda run -n "$ENV_NAME" python -m colabfold.download alphafold2_multimer_v3
fi

echo "==> Checking installation"
conda run -n "$ENV_NAME" colabfold_batch --help >/dev/null
conda run -n "$ENV_NAME" python "$SCRIPT_DIR/phlegm.py" --help >/dev/null

cat <<EOF

PHLEGM installed successfully.

  conda activate $ENV_NAME
  python $SCRIPT_DIR/phlegm.py -i proteins.fasta -o phlegm_out

By default MSAs are generated with the public ColabFold MMseqs2 server.
To use a local ColabFold database instead, add: --db /path/to/colabfold_db
EOF
