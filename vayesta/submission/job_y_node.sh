#!/bin/bash -l                                                                                                                                         
#$ -l h_rt=48:00:00                                                                                                                                    
#$ -l mem=35G                                                                                                                                          
#$ -pe mpi 40                                                                                                                                          
#$ -N mpi_vayesta
#$ -A KCL_YOUR_GROUP                                                                                                                                    
#$ -ac allow=Y                                                                                                                                         
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

export HDF5_INCLUDE_PATH=$HOME/work/ext_libs/include
export HDF5_LIBRARY_PATH=$HOME/work/ext_libs/lib

source $HOME/work/pol/bin/activate

export PATH=$HOME/work/pol/bin:$PATH
export PYTHONPATH=$HOME/work/pol/lib/python3.9/site-packages:$PYTHONPATH

export PYTHONPATH=$PYTHONPATH:$HOME/work/pyscf
export PYTHONPATH=$PYTHONPATH:$HOME/work/Vayesta

OMP_NUM_THREADS=$(ppn)
gerun python 90-mpi.py > output.out
