

CP2K Tutorial for Young users
===============================

In this document, I will introduce the use of **CP2K** a versatile ab-initio modular software to be used for the computation of electronic structure 
calculation.The program is able to perform the following calculations:

1. *Ab-initio* Molecular Dynamics (AIMD) for different ranges of system in both  Periodic Boundary Condition (PBC) and gas-phase systems.

2. Excited state calculations, in different flavours, namely:
  * TD-DFT
  * TD-DFPT
  * LR-TDDFT
  
3. PBC and gas-phase geometry optimizations. The calculations includes Cell optimization, Phonon calculations 
(using the supercell method for non Gamma-points)
  
4. XAS (X-RAy absorption spectraum) for gas-phase and extended (PBC) systems.

5. Plumed Calculations


Among other capabilities.

Structure of this Repo:
^^^^^^^^^^^^^^^^^^^^^^^^^

The different folders contain different kind of working input files. The user can copy and paste the inputs or transfer them to Young
via the git clone command. The calculatations are divided in the following manner:

* Geo_Opt: This folder contains the geometry optimization inputs of isolated (gas-phase) at different levels of theory, namely, PM6, BLYP (GGA-DFT) and 
B3LYP using libxc library.

* SP: Single point calculation using the B3LYP excahnge-corrrelation via LibXC library.

* MD: AIMD calculation of a gas-phase molecule. This empoys the NVT thermostat (Nose-Hoover) at a temperature of 300K.

* TD-DFT: Includes 3 folders, where calculations using TD-DFT, TD-DFPT and LR-TDDFT are used.

All folders contain the file **job.sh**, which includes a template of a submission file to be used with Young.









