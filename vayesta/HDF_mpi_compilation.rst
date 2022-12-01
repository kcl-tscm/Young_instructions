How to install HDF5 and h5py supporting MPI on OS X
=======================================================

Works on hdf5-1.10.5 and modules used for vayesta compilation. This should be done **BEFORE** installing **PySCF** and **Vayesta**.

## HDF5
0. Load all the modules used to install Vayesta such as openmpi, beta-modules, etc. 

1. Download source tarball from [https://www.hdfgroup.org/downloads/hdf5/]. Untar the folder in $HOME. Create a new directory inside 
   $HOME/work directory that can be named ext_libs.
   
.. code-block:: bash

   cd $HOME/work
   mkdir ext_libs
   

2. Configure Makefile using the following command:  

.. code-block:: bash
    
    $ CC=mpiCC ./configure --prefix=$HOME/work/ext_libs --enable-shared --enable-parallel

This will configure a shared library HDF5 installation that uses MPI. Note that this is provided by loading all the requested modules.

3. Compile and install HDF5.  
.. code-block:: bash
  $ make
  $ OMP_NUM_THREADS=4 make check
  $ make install

4. Check HDF5 installation using `h5pcc --showconfig` on the command line.

## h5py
0. Download and install numpy. Use pip or conda.
1. Declare environmental variables used by h5py  

.. code-block:: bash
  
  export CC=mpicc
  export HDF5_MPI="ON"

2. Download h5py from GitHub

.. code-block:: bash

   git clone 

3. Configure using setup.py and install.  
```bash
$ python setup.py configure --mpi
$ python setup.py install
```
3. Check the installation using the demo code on this [page](http://docs.h5py.org/en/latest/mpi.html)
