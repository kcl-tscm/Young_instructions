How to install xtb-6.5.1 in Young
=================================

This is a guide to install the latest version of the semi-empirical code **xtb**. Since this can be used in conjunction with Python scripts, we will
also create a virtual environment. 


Modules required:
------------------
Before starting, please deactivate all modules which are automaticall loaded, by typing:

.. code-block:: bash

   [~] module purge

Subsequently, the following modules to install and submit jobs using **xtb**:

 1. emacs/28.1   
 2. beta-modules   
 3. gcc-libs/10.2.0   
 4. compilers/gnu/10.2.0   
 5. python/3.9.6-gnu-10.2.0   
 6. cmake/3.21.1(default)   
 7. openblas/0.3.13-openmp/gnu-10.2.0
 8. gerun

this is done by using the command:

.. code-block::  bash

   [~] module load MODULE
   
where MODULE is replaced by the aforementined list.

Building XTB-6.5.1 within a virtual environment:
-----------------------------------------------------

The main advantages to build **xtb** within a virtual environment is that we can use the **bin/** folder created by the virtual environment to store the
**xtb** executable. All the previous modules must be loaded at this point. The following steps creates a virtual environment within a folder called 
**work**:

.. code-block:: bash

   [~] mkdir work
   [~] cd work
   [~] python -m venv base
   [~] source base/bin/activate
   [~] export CC=gcc
   [~] export FC=gfortran
   [~] git clone https://github.com/grimme-lab/xtb.git
   [~] cd xtb
   [~] mkdir build
   [~] cmake .. -DCMAKE_BUILD_TYPE=Release
   [~] make -j4
   [~] make test
   [~] chmod +x xtb
   [~] cp -r xtb ../base/bin/
   
If everything went fine, the **make test** command will test the **xtb** executable. All tests must said **passed**. The executable can be found in **/work/base/** and can be used as:

.. code-block:: bash

  [~] ../base/xtb mol.xyz --opt 
  
to run a geometrical optimization, for instance. To submit jobs using the queue system, please send me an email, and I will explain to you the procedure. 

