Installing HDF5 enabling mpi 
=====================================

This library is very important for compiling tbe and unfortunately the version in Young is serial and a mpi version is required. HDF5 needs also 
the **./configure** command and therefore the previous steps are needed for using that command. 

1. Download the library:

In this webpage the library can be downloaded.

.. code-block:: bash

   https://www.hdfgroup.org/downloads/hdf5/
   
A registration is needed to download the library. After this, one can transfer the library from your PC to YOung using scp. 


2. Compiling the library

This can be done via this command:

.. code-block:: bash

   CC=mpicc FC=mpif90 ./configure --prefix=$INSTALL_PREFIX --enable-parallel --enable-shared --enable-fortran  
   make 
   make check-s 
   make check-p RUNPARALLEL='mpirun -n 4'

where the **$INSTALL_PREFIX** has been defined as:

.. code-block:: bash

   $HOME/young_compilation/opt/hdf/gcc/9.2.0/

The commands **make check-s** checks the single-core version and **make check-p RUNPARALLEL='mpirun -n 4 '** checks the parallel tests. The library 
should be installed after these steps.

3. Define the corresponding environmental variables

This can be done using the following command:

.. code-block:: bash

   export PATH="$HOME/young_compilation/opt/hdf/gcc/9.2.0/bin:$PATH"
   export INCLUDE="$HOME/young_compilation/opt/hdf/gcc/9.2.0/include:$INCLUDE"
   export LD_LIBRARY_PATH="$HOME/young_compilation/opt/hdf/gcc/9.2.0/lib:$LD_LIBRARY_PATH"
   export HDF5_DIR=$HOME/young_compilation/opt/hdf

This will provide the required variables to be used by either **cmake** or **./configure**.
