#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -l mem=4.5G
#$ -pe mpi 40   # it can be 40,80,120,160 and os on
#$ -N my_simulation
#$ -A KCL_CHANGE_FOR_YOUR_GROUP_NAME
#$ -P Gold
#$ -cwd 

module purge
module load compilers/gnu/4.9.2
module load orca/5.0.4-sbindist
module load gerun

export OMP_NUM_THREADS=1
ORCA_EXEC=$(which orca)
$ORCA_EXEC input.inp > output.out
