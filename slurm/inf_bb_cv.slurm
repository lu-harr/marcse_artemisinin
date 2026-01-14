#!/bin/bash
#SBATCH --job-name=inf_cv
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=24:00:00           # Max runtime (HH:MM:SS)
#SBATCH --mem=50000                  # Memory per node
#SBATCH --array=0-4
#SBATCH --cpus-per-task=30         # Number of CPU cores per task

#SBATCH --mail-type=BEGIN,END,FAIL  # Send email when job begins, ends, or fails
#SBATCH --mail-user=lucy.harrison@wwarn.org  # Replace with your email address

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate greta-env

cd ../

pwd

DATASETS=(
  mdr86
  mdr184
  mdr1246
  crt76
  k13_marcse
)

# Select dataset for this array task
DATASET=${DATASETS[$SLURM_ARRAY_TASK_ID]}

# Fixed seed ...
SEED=124

# Run R script
Rscript code/fit_betabinom_cv.R "$DATASET" "$SEED"

# provide marker and seed ...
# Rscript code/fit_betabinom_cv.R mdr86 124 &
# Rscript code/fit_betabinom_cv.R mdr184 124 &
# Rscript code/fit_betabinom_cv.R mdr1246 124 &
# Rscript code/fit_betabinom_cv.R crt76 124 &
# Rscript code/fit_betabinom_cv.R k13_marcse 124 &
# wait

echo "Job done!"
