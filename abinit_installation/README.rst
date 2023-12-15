*Installation of abinit :*
==========================

These are the instructions to install abinit (latest version) in Young. 

First of all, we need to unload all modules which are uploaded when we initiate our session. This is done by the command:

.. code-block:: bash

   module purge
   module load emacs

Abinit requires a text editor (in this case is emacs).

*Preparing the environement for installing the code :*
=======================================================

We need to download *abinit*. This is done by typping the following commands:

.. code-block:: bash

    mkdir abinit_latest
    cd abinit_latest
    wget https://www.abinit.org/sites/default/files/packages/abinit-9.10.3.tar.gz
    tar -xzvf abinit-9.10.3.tar.gz

once this is done, we have a folder called *abinit-9.10.3*. We need to access this folder by typping *cd abinit-9.10.3*.

*Loading the required modules :*
=====================================

The following list of modules need to be loaded for installing abinit:

.. code-block:: bash

   module load beta-modules
   module load gcc-libs/10.20.0
   module load compilers/gnu/10.2.0
   module load openblas/0.3.13-openmp/gnu-10.2.0
   module load mpi/openmpi/4.0.5/gnu-10.2.0
   module load fftw/3.3.9/gnu-10.2.0 
   module load cmake
   

*Compiling abinit:*
==========================

Special attention is required in this section, since a preliminary attempt for building abinit is carried out. This attempt will fail, but 
abinit enables the use of a folder called **fallbacks**. This folder contains the missing libraries and it will be automatically compiled.
This is required, since it is more complicated to compile 4 libraries (i.e. hdf5, netcdf, netcdf_fortran, and LibXC) which are not recognised
by abinit from the loaded modules. 

The **failed attempt** can be performed in this way:

.. code-block:: bash

   cd abinit-9.10.3
   ./configure --prefix="$PATH_TO_DIR_WHERE_TO_INSTALL_ABINIT" FC=mpif90 CC=mpicc CXX=mpicxx 

This will start to check the paths and libraries. The flag *$PATH_TO_DIR_WHERE_TO_INSTALL_ABINIT* is defined by the user, namely you need to define where to store the executables of Abinit. You can use the first folder abinit_folder ($HOME/abinit_folder).
At some point is going to break due to not gfortran command found. Once this happens, we need  to type the following commands:

.. code-block:: bash

   cd fallbacks
   ./build-abinit-fallbacks.sh

This will install many libraries, it can take up to 1-hour to complete. All the 4 libraries will be installed. Once this process finishes, we need
to copy the following messages provided by this script:

.. code-block:: rst

   The fallbacks are ready to use.
   You can link these fallbacks with Abinit by copying the following options:

   with_libxc=$PATH_GIVEN_BY_THE_SCRIPT
   with_hdf5=$PATH_GIVEN_BY_THE_SCRIPT
   with_netcdf=$PATH_GIVEN_BY_THE_SCRIPT
   with_netcdf_fortran=$PATH_GIVEN_BY_THE_SCRIPT


Based on this information, we need to proceed in the following manner:

.. code-block:: bash

   cd ..
   ./configure --prefix="$PATH_TO_DIR_WHERE_TO_INSTALL_ABINIT"  --with-libxc="$PATH_GIVEN_IN_THE_FALLBACKS" --with-hdf5="$PATH_GIVEN_IN_THE_FALLBACKS" --with-netcdf="$PATH_GIVEN_IN_THE_FALLBACKS" --with-netcdf_fortran="$PATH_GIVEN_IN_THE_FALLBACKS" FC=mpif90 CC=mpicc CXX=mpicxx 
   
This command should contruct a valid **MakeFile** enabling the compilation of abinit. It is important to notice that there are differences with the path reported by the fallback script (**with_libxc**) and the one declared for the **configure** command (**--with-libxc=**). Also
it is important to notice that for the path in the configure command, one needs to declare it in between quotation marks **"PATH_HERE"**. 

finally, once this is finished. One can compile the code by typping:

.. code-block:: bash

   make -j8

and then one can test the installation using the following command:

.. code-block:: bash

   make test_fast

This is recommended to test that all executables have been correctly built. If this passes, we can perform the most extensive testing of the executable by typping the following 
commands:

.. code-block:: bash

   cd tests
   python runtests.py -j8


This will execute and run the whole test_suite (around 6 different sets of tests). It can take a while (around 2 hours). Finally, the library can be fully installed 
by typping this command:

.. code-block:: bash

   make install -j8

And this will build inside the directory declared in the **prefix** flag the executables and stored in folders named **bin**, **lib**, and **shared**.






