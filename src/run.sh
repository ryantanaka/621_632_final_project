#!/bin/bash

# Parameters
HOSTFILE="./hostfile_64.txt"
PLATFORM="./cluster_crossbar_64.xml"
NP=64

# Compile and Run
smpicc -O3 reduce_skeleton.c greedy_reduce.c binomial_reduce.c -o reduce -lm \
&& \
smpirun --cfg=smpi/reduce:binomial \
  --cfg=smpi/privatization:yes \
  -np $NP \
  -hostfile $HOSTFILE \
  -platform $PLATFORM \
  ./reduce
