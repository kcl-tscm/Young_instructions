#!/bin/bash -l
#$ -l h_rt=1:00:00
#$ -l mem=4G
#$ -pe mpi 40
#$ -N orca
#$ -A KCL_pilot
#$ -P Free
#$ -wd /home/mmm0666/orca


module unload default-modules/2018
module unload compilers
module load gerun
module load gcc-libs
module load compilers/gnu/4.9.2
module load mpi/openmpi/3.1.4/gnu-4.9.2
module load orca/4.2.1-bindist/gnu-4.9.2

export OMP_NUM_THREADS=1
ORCA_EXEC=$(which orca)
$ORCA_EXEC example1.inp > aug.out
