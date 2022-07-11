#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -l mem=4.5G
#$ -pe mpi 80
#$ -N cp2k_blyp_optimization
#$ -A YOUR_GROUP_BUDGET
#$ -P Gold
#$ -cwd 

module purge
module load gerun
module load gcc-libs
module load compilers/gnu/4.9.2
module load mpi/openmpi/3.1.4/gnu-4.9.2
module load openblas/0.3.7-openmp/gnu-4.9.2
module load cp2k/7.1/ompi/gnu-4.9.2

export OMP_NUM_THREADS=1
gerun cp2k.popt -inp geo.inp > output.out 
