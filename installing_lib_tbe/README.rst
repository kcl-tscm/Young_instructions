
Installing libraries using ./configure command
================================================

Some of these libraries use the ./configure command. To do so, one needs to follow these steps to create the command and install 
the library using this manner.

Modules needed:
----------------

These are the modules needed for the ./configure command to be used:

1. beta-modules     
2. compilers/gnu/9.2.0   
3. cmake/3.21.1     
4. guile/2.0.11/gnu-4.9.2     
5. autoconf/2.69         
6. libtool/2.4.6             
7. autogen/5.18.12                        
8. gcc-libs/9.2.0   
9. automake/1.16.1       
10. libbdwgc/7.4.2/  


Steps to follow:
------------------

Within the folder of the downloaded library, one needs to used the following commands **in this order:**

.. code-block:: bash
  
  libtoolize --force
  aclocal
  autoheader 
  automake --force-missing --add-missing
  autoconf 

Once this is done the **./configure** command can be used in conjuction with make to install the library. 

