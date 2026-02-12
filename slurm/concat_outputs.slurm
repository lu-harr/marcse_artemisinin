#!/bin/bash
#SBATCH --job-name=concat
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --time=2:00:00           # Max runtime (HH:MM:SS)
#SBATCH --mem=5G                  # Memory per node
#SBATCH --ntasks=1                # Number of tasks
#SBATCH --cpus-per-task=1         # Number of CPU cores per task

#SBATCH --mail-type=BEGIN,END,FAIL  # Send email when job begins, ends, or fails
#SBATCH --mail-user=lucy.harrison@wwarn.org  # Replace with your email address

# should take under 1 h

module load Miniforge3/24.1.2-0
eval "$(conda shell.bash hook)"
conda activate greta-env

cd ../

pwd

# provide marker and model ...

Rscript code/concat_outputs_slurm.R k13_marcse bb_gne
Rscript code/concat_outputs_slurm.R mdr86 bb_gne
Rscript code/concat_outputs_slurm.R mdr184 bb_gne
Rscript code/concat_outputs_slurm.R mdr1246 bb_gne
Rscript code/concat_outputs_slurm.R crt76 bb_gne
Rscript code/concat_outputs_slurm.R k13snp_A675V bb_gne
Rscript code/concat_outputs_slurm.R k13snp_C469Y bb_gne
Rscript code/concat_outputs_slurm.R k13snp_P441L bb_gne
Rscript code/concat_outputs_slurm.R k13snp_R561H bb_gne
Rscript code/concat_outputs_slurm.R k13snp_R622I bb_gne



echo "Job done!"
