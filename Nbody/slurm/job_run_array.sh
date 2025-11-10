#!/bin/bash
#SBATCH --job-name=N_1D
#SBATCH --output=logs/job_1d_%A_%a.out
#SBATCH --error=logs/job_1d_%A_%a.err
#SBATCH --partition=<your partition>
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4-04:00:00
#SBATCH --mem=2G
#SBATCH --mail-type=TIME_LIMIT,END
#SBATCH --mail-user=<your email>
#SBATCH --array=0-63  # 64 array jobs

module purge
module load julia/1.11.3

# Prevent Julia from spawning additional threads
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Simulation parameters
N=1000
TMAX=2000
SAVE_FREQ=1
MODEL=plummer
OUTPUT=${MODEL}_N_${N}
VERBOSITY=false
SAVE_FINAL_STATE=true

PATH_TO_RUN=/path/to/run
RUN=Main.jl

# Seeds to process
SEED_INIT=0
SEED_FINAL=1000
TOTAL_SEEDS=$((SEED_FINAL - SEED_INIT + 1))
NUM_ARRAY_TASKS=64

# Distribute seeds as evenly as possible
BASE_SEEDS_PER_TASK=$((TOTAL_SEEDS / NUM_ARRAY_TASKS))
REMAINDER=$((TOTAL_SEEDS % NUM_ARRAY_TASKS))

# Compute how many seeds this task gets
if [ $SLURM_ARRAY_TASK_ID -lt $REMAINDER ]; then
    SEEDS_PER_TASK=$((BASE_SEEDS_PER_TASK + 1))
    START_SEED=$((SEED_INIT + SLURM_ARRAY_TASK_ID * SEEDS_PER_TASK))
else
    SEEDS_PER_TASK=$((BASE_SEEDS_PER_TASK))
    START_SEED=$((SEED_INIT + REMAINDER * (BASE_SEEDS_PER_TASK + 1) + \
                  (SLURM_ARRAY_TASK_ID - REMAINDER) * BASE_SEEDS_PER_TASK))
fi

END_SEED=$((START_SEED + SEEDS_PER_TASK - 1))

cd ${PATH_TO_RUN}

# Run Julia sequentially over the subset of seeds for this task
for SEED in $(seq $START_SEED $END_SEED); do
    julia ${RUN} --N ${N} --tmax ${TMAX} --save_freq ${SAVE_FREQ} \
                 --model ${MODEL} --seed ${SEED} --output ${OUTPUT} \
                 --verbose ${VERBOSITY} --save_final_state ${SAVE_FINAL_STATE}
done