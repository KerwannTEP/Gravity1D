#!/bin/bash
#SBATCH --job-name=dos_energy
#SBATCH --partition=pscomp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH --mem=64G
#SBATCH --array=0-11                # 12 jobs total
#SBATCH --output=logs/dos_energy_%A_%a.out
#SBATCH --error=logs/dos_energy_%A_%a.err
#SBATCH --exclusive

module load julia/1.10   # adjust if needed

# Build the list of seeds dynamically
SEEDS=( $(ls ../data/output/ | grep seed_ | sed 's/seed_//' | sort -n) )
NSEEDS=${#SEEDS[@]}

# Distribute seeds across the 12 array tasks
for ((i=${SLURM_ARRAY_TASK_ID}; i<$NSEEDS; i+=12)); do
    SEED=${SEEDS[$i]}
    echo "[$SLURM_ARRAY_TASK_ID] Processing seed $SEED"
    julia analyze_seed.jl --output output --tmax 2000 --N 1000 --seed $SEED
done
