# setupFEP  

This collection of python command line functions is designed with the
aim to faciliate a robust and fast setup of FEP calculations for the
software package Q. It includes modules which are described in further
detail below:

- setupFEP.py: module to generate ligand FEP calculations using a
dual topology approach, see Jespers et al. (doi:/MS_SUBMITTED).
- resFEP.py: module to setup protein FEP calcultions, based on 
Boukharta et al. (https://doi.org/10.1371/journal.pcbi.1003585) and 
Jespers et al. (see doi://MS_SUBMITTED).
- protprep.py: a module to prepare proteins for simulations under
spherical boundary conditions.
- opls2Q and charmm2Q.py: tools to translate from OPLS (ffld_server) 
and CHARMM (cgenff) to Q readable file formats.

A few toplevel scripts are included in the scripts folder to faciliate
high throughput setup. Additionally, a tutorials folder is included
with a detailed description of the setup procedure as published in
Jespers et al. (resFEP/ligFEP). This tutorial includes the generation
of ligand parameters using OPLS, how to prepare a protein system and
how to run ligand and protein FEP calculations. These examples are 
based on ligand binding of CDk2 inhibitors, and protein stability
calculations on T4-lysozyme.

## Requirements  
- ffld_server
- cgenff
- Python2.7.XX
- Q  


