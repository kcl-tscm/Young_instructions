Installing LibXC
=====================

In this section, we will present the procedure to install **LibXC-4.3.4** which is a compatible version for tbe. To use the **./configure**
command, please follow the instructions in the main folder indicating the modules needed and the way to enable that command. I have created 
a folder called  **young_compilation** in **$HOME**, in which the compiled library and the folders are contained. These are the steps to install the
library. 

1. Download the library from this website:

.. code-block:: bash

    https://tddft.org/programs/libxc/download/previous/


Download the version 4.3.4, otherwise some incompatibilities can happen. Use the **scp** command to transfer the folder from your PC to Young. If you
are not sure about this step, please contact me. 

2. Untar the file and enter the folder. After doing this, please type the following commands:

.. code-block:: bash

  ./configure --prefix=PATH/TO/LIBXC
  make
  make check
  make install

**Note:** I have used the PATH/TO/LIBXC as: 

.. code-block:: bash

   $HOME/young_compilation/opt/libxc/4.3.4/gcc/9.2.0/

please change it, accordingly.

3. Create the environment variables to be used for a subsequent program that needs to be installed:

.. code-block:: bash

   export PATH="$HOME/young_compilation/opt/libxc/4.3.4/gcc/9.2.0/bin:$PATH"
   export INCLUDE="$HOME/young_compilation/opt/libxc/4.3.4/gcc/9.2.0/include:$INCLUDE"
   export LD_LIBRARY_PATH="$HOME/young_compilation/opt/libxc/4.3.4/gcc/9.2.0/lib:$LD_LIBRARY_PATH"
   
This last step is very important to link the library to another installation using either **cmake** or **./configure**.
