
Geometry Optimization for molecular systems
====================================================


1. **PM6:** Semi-empirical method PM6 used to quickly optimize molecules using a MT, non-periodic type of calculation. The input declares a variable 
in the beggining of the file calling for the basename of the file. This is used as a FLAG for the &TOPOLOGY section in which the extension of the file 
(Geometry format) is added.


2. **BLYP:** GGA exchange-correlation potential is used to relax the mol.xyz molucular system. 

3. **B3LYP:** Hybrid exchange-correlation potential is used to relax the mol.xyz molucular system. The library LibXC is used to perform this calculation.


Tips for performing these calculations:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


1. For large molecules, it is always good idea to start from the cheapest computational model (i.e PM6 in this case). The relaxed structure
from this method is used as initial guess from a higher-level method (such as BLYP or B3LYP).

2. To perform Hybrid calculations, it is recommened to use an already converged GGA calculation hopefuly same family of XC. For instance, if you
want to perform a B3LYP optimization, then you can use a BLYP GGA wavefunction as starting point. This helps with convergence and avoid errors of
Hartree-Fock not bein occupied 100%.


