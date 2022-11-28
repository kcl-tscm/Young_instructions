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

.. warning::

   There are couple of important things to notice from these instructions. First the command **module purge** unload ALL modules including any text      
   editor such as *emacs*. If you use **vim**, this is the default text editor and therefore is not unloaded. 
   
   The module **python/3.9.6** is the only python installation that does not have any **openblas/numpy serial** linked module. This is important 
   since we want to install **numpy** with **openblas**. 
   
   
   
   
   
   

   
   
   
   
   
   


