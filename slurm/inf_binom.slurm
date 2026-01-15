#!/bin/bash
#SBATCH --job-name=inf_all
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --time=1:15:00           # Max runtime (HH:MM:SS)
#SBATCH --mem=50000                  # Memory per node
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --cpus-per-task=30         # Number of CPU cores per task

#SBATCH --mail-type=BEGIN,END,FAIL  # Send email when job begins, ends, or fails
#SBATCH --mail-user=lucy.harrison@wwarn.org  # Replace with your email address

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate greta-env

cd ../

pwd

# provide marker and seed ...
Rscript code/fit_binom.R mdr86 124 &
Rscript code/fit_binom.R mdr184 124 &
Rscript code/fit_binom.R mdr1246 124 &
Rscript code/fit_binom.R crt76 124 &
Rscript code/fit_binom.R k13_marcse 124 &
wait

echo "Job done!"
