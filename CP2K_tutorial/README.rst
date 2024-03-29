

CP2K Tutorial for Young users
===============================

In this document, I will introduce the use of **CP2K** a versatile ab-initio modular software to be used for the computation of electronic structure 
calculation.The program is able to perform the following calculations:

1. *Ab-initio* Molecular Dynamics (AIMD) for different ranges of system in both  Periodic Boundary Condition (PBC) and gas-phase systems.

2. Excited state calculations using **TD-DFPT**
  
3. PBC and gas-phase geometry optimizations. The calculations includes Cell optimization, Phonon calculations 
(using the supercell method for non Gamma-points)
  
4. XAS (X-RAy absorption spectraum) for gas-phase and extended (PBC) systems.

5. Plumed Calculations


Among other capabilities.

Structure of this Repo:
^^^^^^^^^^^^^^^^^^^^^^^^^

The different folders contain different kind of working input files. The user can copy and paste the inputs or transfer them to Young
via the git clone command. The calculatations are divided in the following manner:

1. **Geo_Opt:** This folder contains the geometry optimization inputs of isolated (gas-phase) at differen levels of theory such as PM6, BLYP (GGA-DFT), 
and B3LYP (via LibXC )

2. **SP:** Single point calculation using the B3LYP excahnge-corrrelation via LibXC library.

3. **MD:** AIMD calculation of a gas-phase molecule. This empoys the NVT thermostat (Nose-Hoover) at a temperature of 300K.

4. **TD-DFT:** where a TDDFPT calcilation is performed.

All folders contain the file **job.sh**, which includes a template of a submission file to be used with Young.









