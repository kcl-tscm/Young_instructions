#!/bin/bash -l                                                                                                                                   
#$ -l h_rt=48:00:00                                                                                                                              
#$ -l mem=4G                                                                                                                                     
#$ -pe smp 40                                                                                                                                    
#$ -N negf_test                                                                                                                                  
#$ -A KCL_                                                                                                                             
#$ -P Gold                                                                                                                                       
#$ -cwd                                                                                                                                          

module purge
module load beta-modules
module load gcc-libs/10.2.0
module load python/3.9.10
module load gerun

source $HOME/matt_negf/negf/bin/activate
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/matt_negf/lib_negf/lib

time OMP_NUM_THREADS=40 ./negf3t
