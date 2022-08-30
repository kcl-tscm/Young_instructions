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

   CC=mpicc FC=mpif90 ./configure --prefix=$INSTALL_PREFIX --enable-parallel --enable-shared --enable-fortran  --enable-fortran2003




