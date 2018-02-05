# setupFEP  

For an example go to the TEST folder and run:  

    python ../generate_prms.py -l TEST -FF OPLS2005 -o Q

After the .prm files have been generated you can run:  

    $QDIR/qprep < Q_TEST.inp > Q_TEST.out

to generate a water solvated system of the TEST.pdb ligand.  

Currently OPLS2005 and OPLS2015 protein parameters and OPLS2005 ligand parameters
are included. For now only transformation to Q parameters is included, converting
to GROMACS .itp files and other forcefield will be added later.  

NOTE: the charge group assignment is currently being developed, therefore the user
is highly encouraged to check the assigned charge groups for ligands in the .lib
file!  


## Requirements  

- Schrodinger's ffld_server.  
- Q  


