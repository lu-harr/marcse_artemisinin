#!/bin/bash
#SBATCH --job-name=pdps
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --time=2:00:00           # Max runtime (HH:MM:SS)
#SBATCH --mem=20G                  # Memory per node
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

# Rscript code/partial_dependence_plots.R k13_marcse
# Rscript code/partial_dependence_plots.R mdr86
# Rscript code/partial_dependence_plots.R mdr184
# Rscript code/partial_dependence_plots.R mdr1246
# Rscript code/partial_dependence_plots.R crt76
# Rscript code/partial_dependence_plots.R k13snp_A675V 
# Rscript code/partial_dependence_plots.R k13snp_C469Y 
# Rscript code/partial_dependence_plots.R k13snp_P441L 
# Rscript code/partial_dependence_plots.R k13snp_R561H
# Rscript code/partial_dependence_plots.R k13snp_R622I

Rscript code/pdps_makeup.R



echo "Job done!"
