#!/usr/bin/env bash
# Submit the PHLEGM pipeline to a SLURM cluster as three chained jobs:
#   1. phlegm msa      one CPU job (local ColabFold database, or the public server)
#   2. phlegm fold     a GPU job array, one task per shard of proteins
#   3. phlegm predict  a short CPU job that runs once all fold tasks finish
#
# Usage:
#   bash examples/submit_slurm.sh proteins.fasta phlegm_out /path/to/colabfold_db [N_SHARDS]
#   bash examples/submit_slurm.sh proteins.fasta phlegm_out server [N_SHARDS]
#
# 'server' uses the public ColabFold MMseqs2 server instead of a local
# database. The MSA job then needs internet access from a compute node; if
# your compute nodes have none, run 'phlegm msa' on the login node first
# (it's light: it mostly waits on the server). The job then finds the MSAs
# already there and does nothing.
#
# The AlphaFold weights are downloaded here first (on the login node), because
# GPU nodes often have no internet access.
#
# Edit the resource settings below (account, partitions, memory, time) for your
# cluster. Install PHLEGM with 'bash install.sh' on Linux, so the same conda
# environment works on CPU and GPU nodes.
set -euo pipefail

if [[ $# -lt 3 ]]; then
    sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi
FASTA=$(realpath "$1")
OUT=$(realpath -m "$2")
if [[ "$3" == "server" ]]; then
    DB_OPT=""
else
    DB_OPT="--db '$(realpath "$3")'"
fi
SHARDS=${4:-10}

# ---- edit for your cluster -------------------------------------------------
ENV_NAME=${PHLEGM_ENV:-phlegm}
ACCOUNT_OPT=""                      # e.g. "--account=myproject"
MSA_OPTS="--cpus-per-task=32 --mem=128G --time=24:00:00"
FOLD_OPTS="--gres=gpu:1 --cpus-per-task=8 --mem=64G --time=24:00:00"   # add --partition=gpu if needed
PREDICT_OPTS="--cpus-per-task=1 --mem=8G --time=01:00:00"
WEIGHTS_DIR=""                      # optional shared AlphaFold weights directory
# ----------------------------------------------------------------------------

RUN="conda run --no-capture-output -n $ENV_NAME phlegm"
WOPT=${WEIGHTS_DIR:+--weights-dir '$WEIGHTS_DIR'}
mkdir -p "$OUT/logs"

echo "Checking AlphaFold weights (downloads ~4 GB the first time)..."
eval "$RUN download-weights $WOPT"

msa=$(sbatch --parsable $ACCOUNT_OPT $MSA_OPTS --job-name=phlegm_msa \
    --output="$OUT/logs/slurm_msa_%j.out" \
    --wrap="$RUN msa -i '$FASTA' -o '$OUT' $DB_OPT --threads \$SLURM_CPUS_PER_TASK")

fold=$(sbatch --parsable $ACCOUNT_OPT $FOLD_OPTS --job-name=phlegm_fold \
    --dependency=afterok:$msa --array=1-$SHARDS \
    --output="$OUT/logs/slurm_fold_%A_%a.out" \
    --wrap="$RUN fold -o '$OUT' --shard \$SLURM_ARRAY_TASK_ID/$SHARDS --keep-going $WOPT")

# afterany: still summarise the proteins that folded if some tasks failed.
predict=$(sbatch --parsable $ACCOUNT_OPT $PREDICT_OPTS --job-name=phlegm_predict \
    --dependency=afterany:$fold \
    --output="$OUT/logs/slurm_predict_%j.out" \
    --wrap="$RUN predict -o '$OUT'")

echo "Submitted: msa=$msa fold=$fold (array 1-$SHARDS) predict=$predict"
echo "Results will be in $OUT/phlegm_predictions.tsv"
