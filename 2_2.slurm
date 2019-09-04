#!/bin/bash 
#SBATCH --ntasks=10
#SBATCH --time=24:00:00
#SBATCH --job-name=2_2
cd /auto/rcf-proj2/ma2/azizim/GTSP
source /usr/usc/julia/1.1.1/setup.sh
julia AdMST.jl 2 2> 2_2.txt 
