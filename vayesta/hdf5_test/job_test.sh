#!/bin/bash -l
#$ -l h_rt=5:00:00
#$ -l mem=4.5G
#$ -pe mpi 80
#$ -N test_mpi_h5py
#$ -A KCL_YOUR_PROJECT
#$ -P Gold
#$ -cwd 

module purge
module load beta-modules
module load gcc-libs/10.2.0
module load python/3.9.6
module load openblas/0.3.13-openmp/gnu-10.2.0
module load compilers/gnu/10.2.0
module load mpi/openmpi/4.0.5/gnu-10.2.0
module load gerun

source $HOME/vayesta_installation/work/bin/activate

export HDF5_INCLUDE_PATH=$HOME/vayesta_installation/libs/include
export HDF5_LIBRARY_PATH=$HOME/vayesta_installation/libs/lib

export OMP_NUM_THREADS=1

gerun python li.py 


