#!/bin/bash
#SBATCH --job-name=pre_bb_cv
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=3:00:00           # Max runtime (HH:MM:SS)
#SBATCH --mem=60G                  # Memory per node
#SBATCH --array=0-49                 # should be up to nmarkers*nfolds
#SBATCH --cpus-per-task=20         # Number of CPU cores per task

#SBATCH --mail-type=BEGIN,END,FAIL  # Send email when job begins, ends, or fails
#SBATCH --mail-user=lucy.harrison@wwarn.org  # Replace with your email address

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate greta-env

pwd

cd ../

markers=(
#   mdr86
#   mdr184
#   mdr1246
#   crt76
#   k13_marcse
  k13snp_R561H
  k13snp_R622I
  k13snp_C469Y
  k13snp_P441L
  k13snp_A675V
)

nfolds=10

idx=$SLURM_ARRAY_TASK_ID

#fold=$(( idx % ${#markers[@]}  + 1))
fold=$(( idx % nfolds  + 1))

# Select marker for this array task
#marker=${markers[idx / ${#markers[@]}]}
marker=${markers[idx / nfolds]}

echo $idx $marker $fold

make_preds() {
    local MARKER=$1
    local MODEL=$2
    local MIN_YEAR=$3
    local MAX_YEAR=$4
    local FOLD=$5
    
    echo "Running task for marker: $MARKER, model: $MODEL, years: $MIN_YEAR to $MAX_YEAR, fold: $FOLD"

    # Loop over the years from MIN_YEAR to MAX_YEAR
    for year in $(seq $MIN_YEAR $MAX_YEAR); do
        echo "Running for year: $year"
    
        # Run the R script with the parameters
        Rscript code/pred_slurm_cv.R $MARKER $MODEL $year $FOLD
    done
}

make_preds "$marker" "bb_gne" 2000 2028 "$fold"


echo "Job done! ${SLURM_ARRAY_TASK_ID}"