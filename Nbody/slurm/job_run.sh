#!/bin/bash
#SBATCH --job-name=N_1D
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --exclusive
#SBATCH --array=0-127
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=pscomp
#SBATCH --time=02:00:00
#SBATCH --mail-type=TIME_LIMIT,END,FAIL
#SBATCH --mail-user=<your email>

module purge
module load julia/1.11.3

N=100
TMAX=5000000
SAVE_FREQ=10000

SEED_INIT=0
MODEL=plummer
OUTPUT=${MODEL}_${N}

SAVE_FINAL_STATE=true

PATH_TO_RUN=/path/to/run/
RUN=Main.jl

SEED=$((SEED_INIT + SLURM_ARRAY_TASK_ID))
RESTART=restart_data_${OUTPUT}_seed_${SEED}.jld2

echo "Running task ${SLURM_ARRAY_TASK_ID} on $(hostname) with seed = $SEED"

cd ${PATH_TO_RUN}

julia ${RUN} --N ${N} --tmax ${TMAX} --save_freq ${SAVE_FREQ} \
            --model ${MODEL} --seed ${SEED} --output ${OUTPUT} \
            --save_final_state ${SAVE_FINAL_STATE} \
            # --restart ${RESTART}

