#!/bin/bash
#SBATCH --job-name=pred_bb
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=8:00:00           # Max runtime (HH:MM:SS)
#SBATCH --mem=60G                  # Memory per node
#SBATCH --array=0-5
#SBATCH --cpus-per-task=30         # Number of CPU cores per task

#SBATCH --mail-type=BEGIN,END,FAIL  # Send email when job begins, ends, or fails
#SBATCH --mail-user=lucy.harrison@wwarn.org  # Replace with your email address

# takes about 6.5h each at the minute

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate greta-env

pwd

cd ../

# Note: if you change what's happening here, be sure to edit array numbers
DATASETS=(
    # mdr86
    # mdr184
    # mdr1246
    # crt76
    k13_marcse
    k13snp_C469Y  
    k13snp_R561H  
    k13snp_A675V  
    k13snp_P441L  
    k13snp_R622I
)

# Select dataset for this array task
DATASET=${DATASETS[$SLURM_ARRAY_TASK_ID]}

# echo "Watch out: filtering unsuccessful runs below"

make_preds() {
    local MARKER=$1
    local MODEL=$2
    local MIN_YEAR=$3
    local MAX_YEAR=$4
    
    echo "Running task for marker: $MARKER, model: $MODEL, years: $MIN_YEAR to $MAX_YEAR"

    # Loop over the years from MIN_YEAR to MAX_YEAR
    for year in $(seq $MIN_YEAR $MAX_YEAR); do

        # sneaky check to hurry things up - remove when I'm ready to rerun
        # if [ ! -f "output/${MARKER}/${MODEL}/${year}_preds.grd" ]; then
        #     echo "Running for year: $year"
    
        # Run the R script with the parameters
        Rscript code/pred_slurm.R $MARKER $MODEL $year
        # fi
        
    done
}

make_preds "$DATASET" "bb_gne" 2000 2028

echo "Job done!"
