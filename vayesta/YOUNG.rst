How to use Vayesta in Young:
==============================

Vayesta_ can be used in Young via the queue (qsub) system. In the folder **submission** the following files can be found:


1. **job_normal.sh :** An example of a submission file employing 1 normal nodes (4.5 GB and 40 cores).

2. **job_normal2.sh :** An example of a submission file employing 2 normal nodes (4.5 GB and 80 cores).

3. **job_y_node.sh :** An example of a submission file employing high-memory RAM nodes (20GB and 40 cores).

4. **job_z_nodes.sh:** An example of a submission file employing very high-memory RAM nodes (60GB and 36 cores).


Disecting the submission file:
=================================

The submission file can be split into three parts:

.. code-block:: bash

  #!/bin/bash -l                                                                                                                                     
  #$ -l h_rt=48:00:00                                                                                                                                 
  #$ -l mem=30G                                                                                                                                           
  #$ -pe mpi 40                                                                                                                                          
  #$ -N mpi_54                                                                                                                                         
  #$ -A KCL_YOUR_GROUP                                                                                                      
  #$ -ac allow=Y       # It can be Y or Z for high-memory jobs. If erased, uses normal nodes.                                                    
  #$ -P Gold                                                       
  #$ -cwd                                                                                                                                                
  
Note:
=======

In this block, the important 


  module purge
  module load beta-modules
module load gcc-libs/10.2.0
module load python/3.9.6
module load openblas/0.3.13-openmp/gnu-10.2.0
module load compilers/gnu/10.2.0
module load mpi/openmpi/4.0.5/gnu-10.2.0

source $HOME/vayesta_installation/pol/bin/activate

export PATH=$HOME/vayesta_installation/pol/bin:$PATH
export PYTHONPATH=$HOME/vayesta_installation/pol/lib/python3.9/site-packages:$PYTHONPATH

export PYTHONPATH=$PYTHONPATH:$HOME/vayesta_installation/pyscf
export PYTHONPATH=$PYTHONPATH:$HOME/vayesta_installation/Vayesta

gerun python 90-mpi.py > output.out

