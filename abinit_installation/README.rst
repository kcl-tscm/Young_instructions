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
   module load gcc-libs/gnu





