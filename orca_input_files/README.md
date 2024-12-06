Contents
============

This repository contains the following files:

1. conformers_ff.py 
2. molecule.xyz
3. orca *.inp files

## How to use conformers_ff ?

For using this code, one needs to install the following python libraries:

``` bash
pip install openbabel-wheel
```
``` bash
pip install rdkit
```
and then you can use the code as python script as:

``` python

python conformers_ff.py
```
## ORCA input files

Many orca input files are added for performing different kind of calculations (single-point (sp), geometry optimisation, excited states, dlpno-ccsd, and freq).
Main blocks are defined for fine tuned parameters in experienced users.

1. orca_sp.inp
2. orca_geo_opt.inp
3. orca_tddft.inp
4. orca_dlpno.inp: Concatenated job. The keyword RICOSX can be also added to speed up calculations. It needs corresponding basis set definition:
   !DLPNO-CCSD(T) def2-pVTZ def2-pVTZ/C RIJCOSX NormalSCF NormalPNO PModel
6. orca_hess.inp
