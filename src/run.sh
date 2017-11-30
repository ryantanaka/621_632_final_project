#!/bin/bash

# Parameters
HOSTFILE="./hostfile_64.txt"
PLATFORM="./cluster_crossbar_64.xml"
NP=64

# Reduce Implementation
S="smpi_binomial_reduce"
B="binomial_reduce"
G="greedy_reduce"

VAR=$1
REDUCE_TYPE=${!VAR}
NUM_INTS=$2

# Compile and Run
smpicc -O3 -DNUM_INTS=$NUM_INTS reduce_skeleton.c greedy_reduce.c binomial_reduce.c -o reduce -lm \
&& \
smpirun --cfg=smpi/reduce:binomial \
  -np $NP \
  -hostfile $HOSTFILE \
  -platform $PLATFORM \
  ./reduce $REDUCE_TYPE
