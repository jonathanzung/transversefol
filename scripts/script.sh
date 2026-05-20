#!/bin/bash

#SBATCH -n 64 #Request 4 tasks (cores)
#SBATCH -N 1 #Request 1 node
#SBATCH -t 0-00:30 #Request runtime of 30 minutes
#SBATCH -p sched_engaging_default #Run on sched_engaging_default partition
#SBATCH --mem-per-cpu=4000 #Request 4G of memory per CPU
#SBATCH -o output_%j.txt #redirect output to output_JOBID.txt
#SBATCH -e error_%j.txt #redirect errors to error_JOBID.txt
#SBATCH --mail-type=BEGIN,END #Mail when job starts and ends
#SBATCH --mail-user=jzung@mit.edu #email recipient

module add julia

JULIA_NUM_THREADS=64 julia search.jl
