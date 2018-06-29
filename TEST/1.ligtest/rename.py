## TEMPORARY RENAME SCRIPT, needs to be added to the generate_parameters script

import glob
import re

include = ['ATOM  ', 'HETATM']
libnames = {}

for mol in glob.glob('*.pdb'):
    name = mol.split('.')[0]
    with open(name + '.lib') as infile:
        for line in infile:
            line = line.split()
            if len(line) > 3:
                try:
                    at_id = int(line[0])
                    at_name = '{:4}'.format(line[1])
                    libnames[at_id] = at_name
                except:
                    continue
    
    with open(name + '.pdb') as infile, open(name + '_out.pdb', 'w') as outfile:
        for line in infile:
            if line[0:6] in include:
                at_id = int(line[7:11])
                atomname = libnames[at_id]
                outfile.write('ATOM  ' + line[6:13] + atomname + 'LIG' + line[20:23] + '  1' + line[26:])
