#!/bin/bash
#SBATCH --job-name=mat52
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --time=02:00:00           # Max runtime (HH:MM:SS)
#SBATCH --mem=5G                  # Memory per node
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --cpus-per-task=1         # Number of CPU cores per task

#SBATCH --mail-type=BEGIN,END,FAIL  # Send email when job begins, ends, or fails
#SBATCH --mail-user=lucy.harrison@wwarn.org  # Replace with your email address


# Load R module (adjust version as needed)
module load R/4.2.0

cd ../
echo pwd

Rscript code/install_all_packages.R
Rscript code/setup.R
Rscript code/mat52.R

echo "Job done!"