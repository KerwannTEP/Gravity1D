#!/bin/bash

N=100
TMAX=500
SAVE_FREQ=1

SEED_INIT=0
SEED_FINAL=1000
MODEL=plummer
OUTPUT=${MODEL}_N_${N}

VERBOSITY=false
SAVE_FINAL_STATE=true

PATH_TO_RUN=../src
RUN=Main.jl

# Number of cores to use
NCORES=12

cd ${PATH_TO_RUN}

seq ${SEED_INIT} ${SEED_FINAL} | xargs -n1 -P${NCORES} -I{} \
    julia ${RUN} --N ${N} --tmax ${TMAX} --save_freq ${SAVE_FREQ} \
                 --model ${MODEL} --seed {} --output ${OUTPUT} \
                 --verbose ${VERBOSITY} --save_final_state ${SAVE_FINAL_STATE}