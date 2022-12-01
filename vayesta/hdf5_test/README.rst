H5PY test to confirm mpi compiled H5PY installation
===========================================================


There are two files in this repository. 

1. **threads.py:** This is a very simple example to produce a h5py data set. 

2. **job_test.sh:** Submission file for using the queue system in Young.


The ouptup can be read by using the command:

.. code-block:: bash

   h5dump parallel_test.hdf5 
   
   
This command is found in **$HOME/work/ext_lib/libs/bin/h5dump**. The obtained result should be this:

.. code-block:: bash

HDF5 "parallel_test.hdf5" {
GROUP "/" {
   DATASET "test" {
      DATATYPE  H5T_STD_I32LE
      DATASPACE  SIMPLE { ( 80 ) / ( 80 ) }
      DATA {
      (0): 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
      (19): 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
      (35): 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
      (51): 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
      (67): 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79
      }
   }
}
}
  
   
   
