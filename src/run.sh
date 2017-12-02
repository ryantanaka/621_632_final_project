#!/bin/bash

# Parameters
HOSTFILE="./hostfile_64.txt"
PLATFORM="./cluster_crossbar_64.xml"
NP=64

# Reduce Implementation
S="smpi_binomial_reduce"
B="binomial_reduce"
P="pipeline_reduce"
G="greedy_reduce"

VAR=$1
REDUCE_TYPE=${!VAR}

# ********* 4 BYTE INTS *********
#  NUM_INTS    MSG SIZE in BYTES
#  --------    -----------------
#  100         400
#  1000        4000 ~ 4 KB
#  10000       40000 ~ 40 KB
#  100000      400000 ~ 400 KB
#  1000000     4000000 ~ 4 MB
#  10000000    40000000 ~ 40 MB
#  100000000   400000000 ~ 400 MB
#
# *******************************
NUM_INTS=$2

# Compile and Run
smpicc -O3 -DNUM_INTS=$NUM_INTS reduce_skeleton.c greedy_reduce.c my_reduce_implementations.c -o reduce -lm \
&& \
smpirun --cfg=smpi/reduce:binomial \
  --cfg=smpi/host-speed:10000000000 \
  -np $NP \
  -hostfile $HOSTFILE \
  -platform $PLATFORM \
  ./reduce $REDUCE_TYPE \
  2> /dev/null
