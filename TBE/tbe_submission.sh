#!/bin/bash -l
#$ -N test_w120_6h_20x2
#$ -P Gold
#$ -A KCL_Paxton
#$ -l h_rt=12:00:00
#$ -l mem=4.5G
#$ -pe mpi 40 
#$ cwd
                                                                        
# In Young, each node has 40 physical cores. To exploit 
# openmp+mpi parallelization, one needs to define how many 
# cores are used to be bound by openmp. This is defined by
# OMP_NUM_THREADS. Depending on the problem it can be, 1 or 20.
# rearely something in between will have gains. 

export OMP_NUM_THREADS=20 #$(ppn)

# gerun is an alias for mpirun -np .... with fine-tuned options of the command
gerun ./run_nvt
