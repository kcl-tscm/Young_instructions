Instructions for installing Vayesta in Young
==============================================

This notes are inteded to create a seamless workflow to install Vayesta_ with the mpi4py support in Young. Along these notes, a set of files is supplied
to perform the installation steps and a submission file for performing calculations using the cluster. 

Files:
----------

1. **Set of modules:** A condensed set of instructions to install Vayesta_ in Young as text file.

2. **job_normal.sh :** An example of a submission file employing normal nodes (4.5 GB and 40 cores).

3. **job_y_node.sh :** An example of a submission file employing high-memory RAM nodes (20GB and 40 cores).

4. **job_z_nodes.sh:** An example of a submission file employing very high-memory RAM nodes (60GB and 36 cores).


How to install Vayesta with mpi support:
==========================================


1. Create a folder in your **$HOME**, it can be called work for instance.

2. After creating the previous folder, the two programs need to be download using git. This is done in the following way:

.. code-block:: bash
   
   [~] cd work
   [~] git clone https://github.com/pyscf/pyscf.git
   [~] git clone https://github.com/BoothGroup/Vayesta.git

After this, two folders will be shown, **pyscf** and **Vayesta**.

3. A set of selected modules are needed to perform the installation:

.. code-block:: bash

   [~] module purge
   [~] module load beta-modules
   [~] module load gcc-libs/10.2.0
   [~] module load python/3.9.6
   [~] module load emacs
   [~] module load openblas/0.3.13-openmp/gnu-10.2.0
   [~] module load compilers/gnu/10.2.0
   [~] module load mpi/openmpi/4.0.5/gnu-10.2.0
   [~] module load cmake/3.21.1


Note:
=======

There are couple of important things to notice from these instructions. First the command **module purge** unload ALL modules including any text      
editor such as *emacs*. If you use **vim**, this is the default text editor and therefore is not unloaded. 
   
The module **python/3.9.6** is the only python installation that does not have any **openblas/numpy serial** linked module. This is important 
since we want to install **numpy** with **openblas**. 
   
   
4. A virtual environment is needed for installing the corresponding libraries for PySCF_ and Vayesta_. This is done in the following way:

.. code-block:: bash
 
   [~] python3 -m venv pol
   [~] source pol/bin/activate

The name **pol** is optional and can be changed. 
   
5. Change the **path** where the Python path where libraries will be searched. This is presented in the following code snippet. 

.. code-block:: bash
 
   [~] export PATH=$HOME/work/pol/bin:$PATH
   [~] export PYTHONPATH=$HOME/work/pol/lib/python3.9/site-packages:$PYTHONPATH

Note:
========

This step is a crucial step to ensure that libraries are installed locally and can be used for further calculations.

6. Installing **mpi4py** with the provided python version. This done as indicated below:

.. code-block:: bash

   [~] python3 -m pip install --upgrade pip
   [~] env MPICC=mpicc python -m pip install --force mpi4py

7. Installing PySCF_ as indicated in the following commands:

.. code-block:: bash

   [~] cd pyscf/
   [~] cd pyscf/lib/
   [~] mkdir build
   [~] cd build
   [~] cmake ..
   [~] make -j8

Note:
=======

In the last command **make -j8**, the option **-j** indicates the number of cores used for the installation. I suggest to use 8 cores
since PYSCF_ builds very heavy libraries such as **libxc**. 


8. Declare the installation path of PySCF_ 

.. code-block:: bash

   [~] PYTHONPATH=$PYTHONPATH:$HOME/work/pyscf

9. Installing Vayesta_ :

.. code-block:: bash

   [~] cd Vayesta
   [~] cd vayesta/libs
   [~] mkdir build
   [~] cd build
   [~] cmake ..

10. Declare the installation path of Vayesta_

.. code-block:: bash

   [~] PYTHONPATH=$PYTHONPATH:$HOME/work/Vayesta/


After these steps have been performed, Vayesta_ and PySCF_ have been installed inside the **bin** folder in the virtual environment created 
within the **$HOME/work** folder.





.. _PySCF: https://pyscf.org/
.. _Vayesta: https://github.com/BoothGroup/Vayesta

.. role:: python(code)
   :language: python

.. role:: console(code)
   :language: console   
   
   
   


