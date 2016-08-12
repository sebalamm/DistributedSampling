#!/bin/bash

# Init parameters
PE=(32 64 128)
SAMPLE=(16 18 20)
POPULATION=30
BASE=10
ITERATIONS=30

python generate_jobs.py -P ${PE[*]} -n ${SAMPLE[*]} -N $POPULATION -k $BASE -i $ITERATIONS 
