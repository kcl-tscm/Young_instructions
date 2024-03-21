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
mpif90 -fbounds-check -o fill.out filling_matrix_ran.f90
```

The code can be run using:

```
mpirun -np #number_of_desired_cores ./fill.out
```

## Example:

with a 4x4 matrix (using two cores), we got this result:

```
4,4
 Filled matrix:
  0.117871046       4.25382257E-02   1.44552302       1.97893417    

  0.332209587      0.111207962       1.35407186       1.07353926    

  0.122348487       3.90575528E-02   1.59612536       1.47379994    

  0.250616789      0.150681615       1.81927860       1.29736590   
```
