#!/bin/bash -l
#$ -l h_rt=1:00:00
#$ -l threads=40
#$ -l mem=4G
#$ -pe wss 1
#$ -N test_nompi
#$ -A KCL_YOURPILASTNAME
#$ -P Gold
#$ -wd ./result_1.out

# Modules

module purge
module load beta-modules
module load gcc-libs/9.2.0
module load compilers/gnu/9.2.0
module load openblas/0.3.7-openmp/gnu-9.2.0
module load mpi/openmpi/3.1.5/gnu-9.2.0

# Write file

touch test_nompi.c
echo "#include<stdlib.h>                                                                                                                               
#include<stdint.h>                                                                                                                   
#include<stdio.h>                                                                                                                              
#include<math.h>                                                                                                                                      
#include<omp.h>                                                                                                                               
#define MAX 100000000                                                                                                                          
int main()                                                                                                 
{                                                                                                                                                                                     
    float sum = 0;                                                                                                                                                                    
                      
#pragma omp parallel for reduction (+:sum)                                                                                                                                            
                      
    for (int i = 0; i < MAX; i++)                                                                                                                                                     
                      
    {                                                                                                                                                                                 
                      
        float x = (i + 0.5) / MAX;                                                                                               
        sum += sqrt(1 - x * x);                                                                             
    }                                                                                                                                                     
}                                                                                                                                        
" >> test_nompi.c

# Compile                                                                                                                                                    
gcc  -Wall -Wextra -std=c11 -O3 -ffast-math  -fPIC -fopenmp -lm  test_nompi.c -o test_nompi

# Run                                                                                                     
time OMP_NUM_THREADS=10 ./test_nompi

# Clean
rm test_nompi.c
rm test_nompi


