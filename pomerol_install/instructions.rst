How to install pomerol in Young
=================================
This document explains step-by-step how to install the library pomarol_ in the HPC facility Young. There are three main steps:

* Installing libcommand_ library
  
* Installing eigen3_ library

* Installing boost_ library 

Importing the required modules to compile the libraries libcommand_, eigen3_ and boost_:

.. code-block:: bash

   module purge		
   module load beta-modules
   module load gcc-libs/9.2.0
   module load compilers/gnu/9.2.0
   module load mpi/openmpi/3.1.5/gnu-9.2.0
   module load cmake/3.21.1
   module load emacs/28.1
   module load openblas/0.3.7-openmp/gnu-9.2.0

Exporting the right compilers and declaring them correctly to be used in the compilation of the 4 libraries:

.. code-block:: bash

   export CC=mpicc
   export CXX=mpicxx


Preparations to install the libraries correctly:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To install pomerol_, one needs to store the libraries in a provided directory. To do so,
one can create a directory called **libs-pomarol**. This should preferebly be done in the
`$HOME` directory. This is done by typing the following command:

.. code-block:: bash

   mkdir $HOME/libs_pomarol
   cd libs_pomarol
   mkdir eigen
   mkdir boost
   mkdir libcommand
   mkdir pomerol

The total path `$HOME/libs_pomarol` will be used as a `PREFIX_PATH` for eigen3_, boost_, and libcommand_, with the corresponding
folders **eigen**, **boost**, and **libcommand**. 


Eigen3 compilation:
^^^^^^^^^^^^^^^^^^^^

The eigen3_ library can be installed by firstly downloading the eigen3 library from the following link_.
After the library has been dowloaded, one needs to send the `.tar.gz` file to Young by using *scp -r eigen-3.4.0.tar.gz mmmXXX@young....*.
Then, one can unpack the library by typing and access the folder by typing:

.. code-block:: bash

   tar -xzvf eigen-3.4.0.tar.gz
   cd eigen-3.4.0

To compile it, one needs to create a new folder and create a new cmake file as shown in these lines of code:

.. code-block:: bash

   mkdir build_dir
   cd build_dir
   cmake ../ -DCMAKE_INSTALL_PREFIX=$HOME/pomerol_libs/eigen
   make install

This suffies for compiling eigen3_.
   
Installing Boost:
^^^^^^^^^^^^^^^^^^

The library boost is more complex to install, and therefore more steps are needed to perform this process. As starting point, one needs to
download the library in link1_ (I have used the version 1.79.0), transfer it to Young (see the scp instructions in eigen3), unpack it, and create
a valid library.  

.. code-block:: bash

   tar -xzvf boost_1_79_0.tar.gz
   cd boost_1_79_0
   ./bootstrap.sh --prefix=$HOME/pomerol_libs/boost/ --with-libraries=all

To use the **MPI**  implementation, one needs to add a line in the script **project-config.jam** which is generated after invoking the command **./bootstrap.sh**.  
This is performed in the following way:

.. code-block:: bash

   emacs project-config.jam
   (after line 12 which states : using gcc ; one needs to include)
   using mpi ;

save the changes and one can proceed with the installation as it follows:

.. code-block:: bash

   ./b2 install --prefix=$HOME/pomerol_libs/boost --target=shared,static
   

Installing libcommand:
^^^^^^^^^^^^^^^^^^^^^^

To install libcommand_, we need to create a new folder called libcom in this tutorial. This is done in `$HOME`.
This can be done in the following way:

.. code-block:: bash

   mkdir $HOME/libcom
   cd libcom
   git clone https://github.com/krivenko/libcommute.git libcommute.src 
   mkdir libcommute.build
   cd libcommute.build
   cmake ../libcommute.src -DCMAKE_INSTALL_PREFIX=$HOME/pomerlo_libs/libcommand -DEXAMPLES=ON -DTESTS=ON 

   
Then, we can create the binary files as shown in the following lines od code:

.. code-block:: bash

   make
   make test
   make install

This creates the library.
   
Installing pomarol:
^^^^^^^^^^^^^^^^^^^

Finally, this is the way to install pomerol:

* We need to download the library from GitHub within a folder created in `$HOME`. This can be done as:

.. code-block:: bash
		
   mkdir $HOME/pomerol_try
   cd pomerol_try
   git clone https://github.com/aeantipov/pomerol.git
   cd pomerol
   mkdir build
   cd build

This downloads and creates the `build` folder where the **CMake** file will be created. Afterwards, we need to link correctly the previously compiled librries. To do so, we need to
follow these steps:

* This is the command to export the eigen3 library correctly:

.. code-block:: bash

   export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$HOME/pomerol_libs/eigen"

Once, we can create the required **CMake** in the following way:

.. code-block:: bash
		
   cmake .. -DCMAKE_BUILD_TYPE=Release -Dlibcommute_DIR=$HOME/pomerol_libs/libcommand/lib/cmake/libcommute/ -DEigen3_DIR=$HOME/pomerol_libs/eigen/share/eigen3/cmake -DBOOST_ROOT=$HOME/pomerol_libs/boost -DCMAKE_INSTALL_PREFIX=$HOME/pomerol_libs/pomerol -DDocumentation=OFF

To create the library one needs to follow these steps:

.. code-block:: bash

   make
   make test
   make install

This creates the library `libpomerol.so`, which can be used to create the Hamiltonians and solve the Green's function.



.. _pomarol: https://github.com/aeantipov/pomerol
.. _boost: https://kratos-wiki.cimne.upc.edu/index.php/How_to_compile_the_Boost_if_you_want_to_use_MPI
.. _libcommand: https://github.com/krivenko/libcommute
.. _eigen3: http://eigen.tuxfamily.org/index.php?title=Main_Page#Compiler_support
.. _link: http://eigen.tuxfamily.org/index.php?title=Main_Page#Download
.. _link1: https://www.boost.org/users/download/

