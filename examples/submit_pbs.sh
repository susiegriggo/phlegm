#!/usr/bin/env bash
# Submit the PHLEGM pipeline to a PBS Pro cluster (e.g. NCI Gadi) as chained jobs:
#   1. phlegm msa      one CPU job with access to the ColabFold database
#   2. phlegm fold     N GPU jobs, one per shard of proteins
#   3. phlegm predict  a short job that runs once all fold jobs finish
#
# Separate fold jobs are used instead of a job array because some PBS sites,
# including Gadi, do not support job arrays.
#
# Usage:
#   bash examples/submit_pbs.sh proteins.fasta phlegm_out /path/to/colabfold_db [N_SHARDS]
#
# Edit the resource settings below (project, queues, storage) for your cluster.
# Install PHLEGM with 'bash install.sh --gpu' so the same conda environment
# works on CPU and GPU nodes.
set -euo pipefail

if [[ $# -lt 3 ]]; then
    sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi
FASTA=$(realpath "$1")
OUT=$(realpath -m "$2")
DB=$(realpath "$3")
SHARDS=${4:-10}

# ---- edit for your cluster -------------------------------------------------
ENV_NAME=${PHLEGM_ENV:-phlegm}
PROJECT=your_project                                     # NCI/PBS project code
COMMON="-P $PROJECT -l storage=gdata/$PROJECT -l wd"     # project and storage (Gadi)
MSA_RES="-q normal -l ncpus=48,mem=190GB,walltime=24:00:00"
FOLD_RES="-q gpuvolta -l ncpus=12,ngpus=1,mem=96GB,walltime=24:00:00"
PREDICT_RES="-q normal -l ncpus=1,mem=8GB,walltime=01:00:00"
# Command that makes 'conda' available inside a batch job:
CONDA_SETUP='eval "$(conda shell.bash hook)"'
# ----------------------------------------------------------------------------

mkdir -p "$OUT/logs"
submit() {  # submit <name> <resources> <depend-option> <command>
    echo "$CONDA_SETUP; conda activate $ENV_NAME; $4" | \
        qsub -N "$1" $COMMON $2 $3 -j oe -o "$OUT/logs/pbs_$1.out"
}

msa=$(submit phlegm_msa "$MSA_RES" "" \
    "phlegm msa -i '$FASTA' -o '$OUT' --db '$DB' --threads \$PBS_NCPUS")

folds=()
for i in $(seq 1 "$SHARDS"); do
    folds+=("$(submit "phlegm_fold_$i" "$FOLD_RES" "-W depend=afterok:$msa" \
        "phlegm fold -o '$OUT' --shard $i/$SHARDS --keep-going")")
done

# afterany: still summarise the proteins that folded if some jobs failed.
deps=$(IFS=:; echo "${folds[*]}")
predict=$(submit phlegm_predict "$PREDICT_RES" "-W depend=afterany:$deps" \
    "phlegm predict -o '$OUT'")

echo "Submitted: msa=$msa, $SHARDS fold jobs, predict=$predict"
echo "Results will be in $OUT/phlegm_predictions.tsv"
