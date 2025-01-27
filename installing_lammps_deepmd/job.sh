#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -l mem=4.5G
#$ -pe mpi 40
#$ -N mpi_6
#$ -A KCL_YOUR_GROUP
#$ -P Gold
#$ -cwd

module purge
module load beta-modules
module load gcc-libs/10.2.0
module load python/3.9.6-gnu-10.2.0
module load openblas/0.3.13-openmp/gnu-10.2.0
module load compilers/gnu/10.2.0
module load mpi/openmpi/4.0.5/gnu-10.2.0
module load fftw/3.3.9/
module load gerun

source $HOME/new_deepmdkit/envdeep/bin/activate

export PATH=$HOME/new_deepmdkit/envdeep/bin:$PATH
export PYTHONPATH=$HOME/new_deepmdkit/envdeep/lib/python3.9/site-packages:$PYTHONPATH

gerun $HOME/new_deepmdkit/mylammps/src/lmp_mpi < many_particles.lam
