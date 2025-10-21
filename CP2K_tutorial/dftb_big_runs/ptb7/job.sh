#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -l mem=4.5G
#$ -pe mpi 40
#$ -N cp2k_optimization
#$ -A KCL_Lorenz
#$ -P Gold
#$ -cwd 

module purge
module load python/miniconda3/24.3.0-0
source $UCL_CONDA_PATH/etc/profile.d/conda.sh

conda activate cp2k_env

export OMP_NUM_THREADS=1
mpirun -np 40 cp2k.psmp -inp md.inp > output.out  
