#!/bin/bash
#SBATCH --job-name=pred_bin
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log
#SBATCH --time=30:00:00           # Max runtime (HH:MM:SS)
#SBATCH --mem=50G                  # Memory per node
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
(for year in $(seq 2000 2027);
do
    echo $year
    Rscript code/pred_slurm.R mdr86 gneiting_sparse $year
done;
Rscript code/concat_outputs_slurm.R mdr86 gneiting_sparse) &
(for year in $(seq 2000 2027);
do
    echo $year
    Rscript code/pred_slurm.R mdr184 gneiting_sparse $year
done;
Rscript code/concat_outputs_slurm.R mdr184 gneiting_sparse) &
(for year in $(seq 2000 2027);
do
    echo $year
    Rscript code/pred_slurm.R mdr1246 gneiting_sparse $year
done;
Rscript code/concat_outputs_slurm.R mdr1246 gneiting_sparse) &
(for year in $(seq 2000 2027);
do
    echo $year
    Rscript code/pred_slurm.R crt76 gneiting_sparse $year
done;
Rscript code/concat_outputs_slurm.R crt76 gneiting_sparse) &
(for year in $(seq 2000 2027);
do
    echo $year
    Rscript code/pred_slurm.R k13_marcse gneiting_sparse $year
done;
Rscript code/concat_outputs_slurm.R k13_marcse gneiting_sparse) &

wait

echo "Job done!"
