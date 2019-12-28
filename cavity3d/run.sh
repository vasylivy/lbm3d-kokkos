#!/bin/bash

release=../lbm3d/release
[ ! -d "output" ] && mkdir output

nx=128; 
ny=128;
nz=128;
re=400; # Reynolds number
tol=0.001; # steady state tolerance 
max_steps=1000000;
output_rate=5000; 

export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=32

$release/lbm.host $nx $ny $nz $re $tol $max_steps $output_rate
