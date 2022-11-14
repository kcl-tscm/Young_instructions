Instructions:
===================

Steps to install fftw_ in Young.

1. To install the fftw_ library, first one needs to dowload it from:

.. code-block:: bash

   http://www.fftw.org/download.html

the current version is **3.3.10**, but this can change. If you dont know how to submit a tar file to Young, contact me for a tutorial
of using the scp command. 


2. After following the instructions in the previous folder to enable **./configure**, you can type the following 
commands:

.. code-block:: bash

  ./configure --prefix=$HOME/PATH_TO_FOLDER --enable-mpi --enable-shared=yes --disable-doc
   make
   make install


It is important to mention that **PATH_TO_FOLDER** must be declared using the **$HOME** variable, and the path to a user-created folder where the 
compiled files and the folders **bin, include, share and libs** will be placed. I strongly recommend to create a folder in **$HOME** and call it 
**my_conf** meaning you are storing those files there and you can rapidly link the libraries when you are compiling another code that needs them.




.. _fftw: http://www.fftw.org/download.html








