#!/bin/bash
#SBATCH --job-name=pdps
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --time=2:00:00           # Max runtime (HH:MM:SS)
#SBATCH --mem=40G                  # Memory per node
#SBATCH --ntasks=1                # Number of tasks

#SBATCH --mail-type=BEGIN,END,FAIL  # Send email when job begins, ends, or fails
#SBATCH --mail-user=lucy.harrison@wwarn.org  # Replace with your email address

# should take under 1 h

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate greta-env

cd ../

pwd

# provide marker and model ...

Rscript code/coverages.R



echo "Job done!"
