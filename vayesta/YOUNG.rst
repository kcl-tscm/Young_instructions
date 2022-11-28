How to use Vayesta in Young:
==============================

Vayesta_ can be used in Young via the queue (qsub) system. In the folder **submission** the following files can be found:


1. **job_normal.sh :** An example of a submission file employing 1 normal nodes (4.5 GB and 40 cores).

2. **job_normal2.sh :** An example of a submission file employing 2 normal nodes (4.5 GB and 80 cores).

3. **job_y_node.sh :** An example of a submission file employing high-memory RAM nodes (20GB and 40 cores).

4. **job_z_node.sh:** An example of a submission file employing very high-memory RAM nodes (60GB and 36 cores).

which display different possibilities of calculations using Vayesta_ in the Young submission file system with **mpi** parallelism.

Disecting the submission file:
=================================

1. The submission file can be split into three parts:

.. code-block:: bash

  #!/bin/bash -l                                                                                                                                     
  #$ -l h_rt=48:00:00                                                                                                                                 
  #$ -l mem=30G        # This must be changed accordingly to the used node Y=30GB, Z=60GB and normal= 4.5 GB                                             
  #$ -pe mpi 40        # Number of processes used. It must be multiple of 40, except Z nodes which is 36 cores instead.                                   
  #$ -N mpi_54         # Name provided to this job.                                                                                                       
  #$ -A KCL_YOUR_GROUP                                                                                                      
  #$ -ac allow=Y       # It can be Y or Z for high-memory jobs. If erased, uses normal nodes.                                                    
  #$ -P Gold                                                       
  #$ -cwd                                                                                                                                                
  
Note:
=======

In this block, the important lines are **-ac** which defines the type of node you are using, the option **-A** for setting your
own KCL project budget 

2. Second part of the submission file:

.. code-block:: bash

  module purge
  module load beta-modules
  module load gcc-libs/10.2.0
  module load python/3.9.6
  module load openblas/0.3.13-openmp/gnu-10.2.0
  module load compilers/gnu/10.2.0
  module load mpi/openmpi/4.0.5/gnu-10.2.0
  module load gerun

Note:
=======

In this block, the required modules are loaded. There is nothing to be changed in this block.

3. Third part of the submission file:

.. code-block:: bash

  # Activate the Python virtual environment
  source $HOME/work/pol/bin/activate

  # Declare the PATHS where PySCF and VAYESTA can be found within the virtual environment.
  export PATH=$HOME/work/pol/bin:$PATH
  export PYTHONPATH=$HOME/vayesta_installation/pol/lib/python3.9/site-packages:$PYTHONPATH

  # Declare the library paths to the python path environment variable
  export PYTHONPATH=$PYTHONPATH:$HOME/vayesta_installation/pyscf
  export PYTHONPATH=$PYTHONPATH:$HOME/vayesta_installation/Vayesta

  gerun python vayesta-mpi.py > output.out

Note:
=======

This is the most important part of the submission file. You need to be aware that **path declaration** is very important, since Vayesta_ is installed
**locally**. These three commands are enabling the local environment to be found by the python interpreter once it is launched by gerun. 

gerun is a wrapper of the command **mpirun -np**.  

.. _Vayesta: https://github.com/BoothGroup/Vayesta

