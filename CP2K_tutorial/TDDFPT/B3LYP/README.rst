CP2K Excited states calculation via TDDFPT
======================================================


In this section, there is an input template to perform time-dependent density function perturbation theory as implemented in CP2K. 
Issue to be considered:

1. The default calculation is **SINGLETS** (multiplicity 1).  
2. To compute **TRIPLETS** (multiplicity 3), one needs to activate the following keyword:

.. code-block:: bash

  &PROPERTIES
    &TDDFPT
       NSTATES      5
       MAX_ITER    100
       CONVERGENCE [eV] 1.0e-6
       **RKS_TRIPLETS .TRUE.**
    &END TDDFPT
  &END PROPERTIES

3. To compute both **SINGLETS** and **TRIPLETS** excitations you need to carry out separated jobs with the keyword **RKS_TRIPLETS ** set to 
**FALSE** for singlets and **TRUE** for triplets.


Validation of the results
==========================

A complementary calculation has been carried out using **ORCA 4.2.1** employing the same Exchange-Correlation potential (B3LYP). There are substantial
differences between the jobs, since **ORCA** is an all electron code, whereas **CP2K** uses a pseudo-potential (GTH in this case) approach. The following
table reports the results between the two programs:

.. list-table:: Title
   :widths: 25 25 50
   :header-rows: 1

   * - Heading row 1, column 1
     - Heading row 1, column 2
     - Heading row 1, column 3
   * - Row 1, column 1
     -
     - Row 1, column 3
   * - Row 2, column 1
     - Row 2, column 2
     - Row 2, column 3



    


