#!/bin/bash
#SBATCH --job-name=inf_all
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err
#SBATCH --time=5:00:00           # Max runtime (HH:MM:SS)
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
# Rscript code/fit_betabinom.R mdr86 124 &
# Rscript code/fit_betabinom.R mdr184 124 &
# Rscript code/fit_betabinom.R mdr1246 124 &
# Rscript code/fit_betabinom.R crt76 124 &
# Rscript code/fit_betabinom.R k13_marcse 124 &
Rscript code/fit_betabinom.R k13snp_R561H 124 &
Rscript code/fit_betabinom.R k13snp_R622I 124 &
Rscript code/fit_betabinom.R k13snp_A675V 124 &
Rscript code/fit_betabinom.R k13snp_P441L 124 &
Rscript code/fit_betabinom.R k13snp_C469Y 124 &
wait

echo "Job done!"
