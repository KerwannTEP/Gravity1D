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
export JULIA_PRECOMPILE_THREADS=1

# Disable automatic precompilation in workers
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NO_PRECOMPILE=1
export JULIA_LOAD_PRECOMPILED_CODE=1

# ------------------------
# Simulation parameters
# ------------------------
N=1000
TMAX=2000
SAVE_FREQ=1
MODEL=plummer
OUTPUT=${MODEL}_N_${N}
VERBOSITY=false
SAVE_FINAL_STATE=true

# ------------------------
# Paths
# ------------------------
PATH_TO_RUN=/path/to/run
RUN=Main.jl

# ------------------------
# Seeds
# ------------------------
SEED_INIT=0
SEED_FINAL=1000
TOTAL_SEEDS=$((SEED_FINAL - SEED_INIT + 1))
NUM_ARRAY_TASKS=64

# Distribute seeds evenly across array tasks
BASE_SEEDS_PER_TASK=$((TOTAL_SEEDS / NUM_ARRAY_TASKS))
REMAINDER=$((TOTAL_SEEDS % NUM_ARRAY_TASKS))

if [ $SLURM_ARRAY_TASK_ID -lt $REMAINDER ]; then
    SEEDS_PER_TASK=$((BASE_SEEDS_PER_TASK + 1))
    START_SEED=$((SEED_INIT + SLURM_ARRAY_TASK_ID * SEEDS_PER_TASK))
else
    SEEDS_PER_TASK=$((BASE_SEEDS_PER_TASK))
    START_SEED=$((SEED_INIT + REMAINDER * (BASE_SEEDS_PER_TASK + 1) + \
                  (SLURM_ARRAY_TASK_ID - REMAINDER) * BASE_SEEDS_PER_TASK))
fi

END_SEED=$((START_SEED + SEEDS_PER_TASK - 1))


# Create a separate environment per N
BASE_ENV="/path/to/LocalEnv"
JULIA_ENV="${BASE_ENV}_N${N}"
JULIA_DEPOT_PATH="/path/to/data/julia_depot_N${N}"  # Job-specific depot

if [ ! -d "$JULIA_ENV" ]; then
    echo ">>> Creating separate environment for N=${N}"
    cp -rn "$BASE_ENV" "$JULIA_ENV"
fi

# Create depot directory if it doesn't exist
mkdir -p "$JULIA_DEPOT_PATH"

# Precompile marker and lockfile
PRECOMPILE_MARKER="$JULIA_ENV/.precompiled_ok"
LOCKFILE="$JULIA_ENV/.precompile_lock"

# ------------------------
# Conditional precompilation (safe for array jobs)
# ------------------------
if [ ! -f "$PRECOMPILE_MARKER" ]; then
    if ( set -o noclobber; echo "$$" > "$LOCKFILE" ) 2> /dev/null; then
        trap 'rm -f "$LOCKFILE"; exit $?' INT TERM EXIT

        echo ">>> SLURM task $SLURM_ARRAY_TASK_ID is precompiling environment for N=${N}..."
        julia --project="$JULIA_ENV" -e '
            using Pkg
            Pkg.instantiate()
            Pkg.precompile()
        '

	touch "$PRECOMPILE_MARKER"
        rm -f "$LOCKFILE"
        trap - INT TERM EXIT
    else
        echo ">>> Waiting for precompilation by another task for N=${N}..."
        while [ ! -f "$PRECOMPILE_MARKER" ]; do
            sleep 2
        done
	echo ">>> Precompilation finished; continuing for N=${N}."
    fi
fi

# ------------------------
# Run simulation over assigned seeds
# ------------------------
cd ${PATH_TO_RUN}

for SEED in $(seq $START_SEED $END_SEED); do
    julia --project=${JULIA_ENV} \
          --compiled-modules=yes --pkgimages=no \
          ${RUN} \
          --N ${N} --tmax ${TMAX} --save_freq ${SAVE_FREQ} \
          --model ${MODEL} --seed ${SEED} --output ${OUTPUT} \
          --verbose ${VERBOSITY} --save_final_state ${SAVE_FINAL_STATE} #\
#          --restart restart_data_${OUTPUT}_seed_${SEED}.jld2
done








