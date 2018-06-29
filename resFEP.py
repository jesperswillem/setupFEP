import argparse
import re
import os
import sys

import functions as f
import settings as s
import IO

class Run(object):
    """
    Setup residue FEPs using either a single or dual topology approach.
    """
    def __init__(self, cofactor, mutation, include, forcefield, *args, **kwargs):
        self.cofactor = [cofactor]
        # Check whether all required files are there:
        required = ['protein.pdb', 'water.pdb', 'protPREP.log']
        if self.cofactor != None:
            extension = ['.pdb', '.lib', '.prm']
            for filename in self.cofactor:
                for line in extension:
                    required.append(filename + line)
        for filename in required:
            if os.path.exists(filename) == False:
                required = ' '.join(required)
                print "ERROR: {} are required files, exiting now.".format(required)
                sys.exit()
        
        self.mutation = re.split('(\d+)', mutation)
        self.include = include
        self.forcefield = forcefield
        
        
        self.sphere = None
        self.CYX = []
        self.PDB = {}
        self.PDB2Q = {}
        
    def read_input(self):
        block = 0
        with open('protPREP.log') as infile:
            for line in infile:
                line = line.split()
                if len(line) > 1:
                    if line[0] == 'Sphere':
                        self.sphere = line[2:]
                        
                    if line[0] == 'Q_CYS1':
                        block = 1
                        
                    if line[0] == 'pdbfile':
                        block = 2
                        
                    if line[0][0] == '-':
                        block = 0
                        
                    if block == 1:
                        if line[0].isdigit():
                            self.CYX.append([line[0], line[1]])
                        
                    if block == 2:
                        self.PDB2Q[line[1]] = line[0]
    
    def readpdb(self):
        with open('protein.pdb') as infile:
            for line in infile:
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
                    try:
                        self.PDB[line[6]].append(line)
                        
                    except:
                        self.PDB[line[6]] = [line]
                       
    def read_prm(self):
        prmfiles = [s.FF_DIR + '/' + self.forcefield + '.prm']
        if self.cofactor != None:
            for filename in self.cofactor:
                prmfiles.append(filename + '.prm')
            
        print IO.read_prm(prmfiles)
        
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='protPREP',
        version='1.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Generate FEP files for resFEP, takes output from protPREP.py == ')

    
    parser.add_argument('-c', '--cofactor',
                        dest = "cofactor",
                        required = False,
                        help = "PDB files of cofactors (e.g. ligand) with *.prm and *.lib files")
    
    parser.add_argument('-m', '--mutation',
                        dest = "mutation",
                        required = True,
                        help = "The desired mutation given as $WT$RESN$MUT,e.g. Y197A, resn is taken from original .pdb")
    
    parser.add_argument('-f', '--forcefield',
                        dest = "forcefield",
                        required = True,
                        choices = ['OPLS2015', 'OPLS2005'],
                        help = "Forcefield to use.")    
    
    args = parser.parse_args()
    run = Run(mutation = args.mutation,
              cofactor = args.cofactor,
              forcefield = args.forcefield,
              include = ('ATOM','HETATM')
             )
    
    run.readpdb()
    run.read_input()
    run.read_prm()