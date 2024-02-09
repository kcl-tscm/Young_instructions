Instructions
===================

This repository includes the input file (input.inp) to perform cluster optimization for Mo clusters using CP2K. 
Also contain an initial structure (mo.xyz) to test and a submission file (job.sh) for submitting jobs.

Changing the structure:
=================================

CP2K is able to read files with different names, to change the name you need to go to this senction in the input file (input.inp):


```properties
    &TOPOLOGY
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME mo.xyz
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
    &END TOPOLOGY
``` 
