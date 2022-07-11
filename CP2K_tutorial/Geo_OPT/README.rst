
Geometry Optimization for molecular systems
====================================================


1. **PM6:** Semi-empirical method PM6 used to quickly optimize molecules using a MT, non-periodic type of calculation. The input declares a variable 
in the beggining of the file calling for the basename of the file. This is used as a FLAG for the &TOPOLOGY section in which the extension of the file 
(Geometry format) is added.


2. **BLYP:** GGA exchange-correlation potential is used to relax the mol.xyz molucular system. 

3. **B3LYP:** Hybrid exchange-correlation potential is used to relax the mol.xyz molucular system. The library LibXC is used to perform this calculation.
