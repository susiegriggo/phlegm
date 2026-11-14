#!/usr/bin/env bash
# Submit the PHLEGM pipeline to a PBS Pro cluster (e.g. NCI Gadi) as chained jobs:
#   1. phlegm msa      one CPU job (local ColabFold database, or the public server)
#   2. phlegm fold     N GPU jobs, one per shard of proteins
#   3. phlegm predict  a short job that runs once all fold jobs finish
#
# Separate fold jobs are used instead of a job array because some PBS sites,
# including Gadi, do not support job arrays.
#
# Usage:
#   bash examples/submit_pbs.sh proteins.fasta phlegm_out /path/to/colabfold_db [N_SHARDS]
#   bash examples/submit_pbs.sh proteins.fasta phlegm_out server [N_SHARDS]
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
# Edit the resource settings below (project, queues, storage) for your cluster.
# Install PHLEGM with 'bash install.sh' on Linux, so the same conda environment
# works on CPU and GPU nodes.
set -euo pipefail

if [[ $# -lt 3 ]]; then
    sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'
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
PROJECT=your_project                                     # NCI/PBS project code
COMMON="-P $PROJECT -l storage=gdata/$PROJECT -l wd"     # project and storage (Gadi)
MSA_RES="-q normal -l ncpus=48,mem=190GB,walltime=24:00:00"
FOLD_RES="-q gpuvolta -l ncpus=12,ngpus=1,mem=96GB,walltime=24:00:00"
PREDICT_RES="-q normal -l ncpus=1,mem=8GB,walltime=01:00:00"
WEIGHTS_DIR=""                      # optional shared AlphaFold weights directory
# Command that makes 'conda' available inside a batch job:
CONDA_SETUP='eval "$(conda shell.bash hook)"'
# ----------------------------------------------------------------------------

WOPT=${WEIGHTS_DIR:+--weights-dir '$WEIGHTS_DIR'}
mkdir -p "$OUT/logs"

echo "Checking AlphaFold weights (downloads ~4 GB the first time)..."
eval "conda run --no-capture-output -n $ENV_NAME phlegm download-weights $WOPT"

submit() {  # submit <name> <resources> <depend-option> <command>
    echo "$CONDA_SETUP; conda activate $ENV_NAME; $4" | \
        qsub -N "$1" $COMMON $2 $3 -j oe -o "$OUT/logs/pbs_$1.out"
}

msa=$(submit phlegm_msa "$MSA_RES" "" \
    "phlegm msa -i '$FASTA' -o '$OUT' $DB_OPT --threads \$PBS_NCPUS")

folds=()
for i in $(seq 1 "$SHARDS"); do
    folds+=("$(submit "phlegm_fold_$i" "$FOLD_RES" "-W depend=afterok:$msa" \
        "phlegm fold -o '$OUT' --shard $i/$SHARDS --keep-going $WOPT")")
done

# afterany: still summarise the proteins that folded if some jobs failed.
deps=$(IFS=:; echo "${folds[*]}")
predict=$(submit phlegm_predict "$PREDICT_RES" "-W depend=afterany:$deps" \
    "phlegm predict -o '$OUT'")

echo "Submitted: msa=$msa, $SHARDS fold jobs, predict=$predict"
echo "Results will be in $OUT/phlegm_predictions.tsv"
