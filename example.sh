#!/bin/bash

SAMPLE_SIZE=5
POPULATION_SIZE=32
BASE_SIZE=100
SEED=1
ITERATIONS=10

echo ./optimized/run_experiments --n=$SAMPLE_SIZE --N=$POPULATION_SIZE --k=$BASE_SIZE --seed=$SEED --i=$ITERATIONS
./optimized/run_experiments --n=$SAMPLE_SIZE --N=$POPULATION_SIZE --k=$BASE_SIZE --seed=$SEED --i=$ITERATIONS
