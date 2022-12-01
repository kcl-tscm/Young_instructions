How to install HDF5 and h5py supporting MPI on OS X
=======================================================

Works on hdf5-1.10.5 and modules used for Vayesta compilation. This should be done **BEFORE** installing **PySCF** and **Vayesta**.
Activate also the **Python virtual environment**. 

## HDF5
0. Load all the modules used to install Vayesta such as openmpi, beta-modules, etc. 

1. Download source tarball from [https://www.hdfgroup.org/downloads/hdf5/]. Untar the folder in $HOME. Create a new directory inside 
   $HOME/work directory that can be named ext_libs.
   
.. code-block:: bash

   cd $HOME/work
   mkdir ext_libs
   

2. Configure Makefile using the following command:  

.. code-block:: bash
    
    CC=mpiCC ./configure --prefix=$HOME/work/ext_libs --enable-shared --enable-parallel

This will configure a shared library HDF5 installation that uses MPI. Note that this is provided by loading all the requested modules.

3. Compile and install HDF5.  

.. code-block:: bash

  $ make
  $ OMP_NUM_THREADS=4 make check
  $ make install

4. Check HDF5 installation using `h5pcc --showconfig` on the command line.

## h5py
0. Download and install numpy. Activate the virtual environment as done in any Vayesta calculation.

1. Declare environmental variables used by h5py  

.. code-block:: bash
  
  export CC=mpicc
  export HDF5_MPI="ON"

2. Download h5py from GitHub

.. code-block:: bash

   git clone https://github.com/h5py/h5py.git

3. Configure using setup.py and install.  

.. code-block:: bash
   
   CC="mpicc" HDF5_MPI="ON" HDF5_DIR=$HOME/vayesta_installation/libs/ pip install --no-cache-dir --no-binary=h5py h5py


3. Check the installation using the demo codes stored in hdf5_test.
