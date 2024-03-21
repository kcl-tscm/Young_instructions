# Filling in parallel a matrix with random numbers

This code is compiled using the following libraries in Young:

```
module purge
module load beta-modules
module load gcc-libs/10.2.0
module load compilers/gnu/10.2.0
module load mpi/openmpi/4.0.5/gnu-10.2.0
module load emacs/28.1
```

The compilation command is:

```
mpif90 -fbounds-check -o fill.out .f90
```
