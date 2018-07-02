import argparse
import re
import os
import sys
import glob

import functions as f
import settings as s
import IO

class Run(object):
    """
    Setup residue FEPs using either a single or dual topology approach.
    """
    def __init__(self, cofactor, mutation, include, forcefield, windows,
                 sampling, system, cluster, temperature, replicates, *args, **kwargs):
        
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
        
        if windows == None:
            self.windows = s.WINDOWS
            
        else:
            self.windows = int(windows)
            
        if sampling == None:
            self.sampling = s.SAMPLING
            
        else:
            self.sampling = sampling
            
        if temperature == None:
            self.temperature = s.TEMPERATURE
            
        else:
            self.temperature = temperature
            
        if replicates == None:
            self.replicates = s.REPLICATES
            
        else:
            self.replicates = replicates
        
        self.mutation = re.split('(\d+)', mutation)
        self.include = include
        self.forcefield = forcefield
        self.lambdas = []
        self.CYX = []
        self.PDB = {}
        self.PDB2Q = {}
        self.systemsize = 0
        self.system = system
        self.cluster = cluster
        self.FEPlist = []
        
        self.replacements = {'EQ_LAMBDA': '1.000 0.000',
                     'ATOM_START_LIG1':'1',
                     'ATOM_END': '{}'.format(self.systemsize),
                     'WATER_RESTRAINT':'',
                     'TEMP_VAR':self.temperature,
                     'RUN_VAR':self.replicates,
                     'RUNFILE':'run' + self.cluster + '.sh'  
                    }
    
    def checkFEP(self):
        AA_from = IO.AA(self.mutation[0])
        AA_to = IO.AA(self.mutation[2])
        mutation = '{}{}{}'.format(*self.mutation)
        FEPdir = s.FF_DIR + '/.FEP/' + AA_from + '_' + AA_to
        if self.forcefield == 'OPLS2005':
            FEPdir = FEPdir + '-aa'
            
        if not os.path.exists(FEPdir):
            print 'FATAL: no FEP files found for the {} mutation in {} exiting now.'.format(mutation,
                                                                                            FEPdir)
            sys.exit()
        
        else:
            for line in sorted(glob.glob(FEPdir + '/FEP*.fep')):
                line = line.split('/')
                self.FEPlist.append(line[-1])
            self.FEPdir = FEPdir
        
    def create_environment(self):
        self.directory = 'FEP_{}{}{}'.format(*self.mutation)
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
            os.makedirs(self.directory + '/inputfiles')
    
    def read_input(self):
        block = 0
        with open('protPREP.log') as infile:
            for line in infile:
                line = line.split()
                if len(line) > 1:
                    if line[1] == 'center:':
                        self.sphere = line[2:]
                        
                    if line[1] == 'radius:':
                        self.radius = line[2]
                        self.replacements['SPHERE'] = self.radius
                        
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
                    
                    self.systemsize += 1
                       
    def merge_prm(self):
        headers =['[options]', 
                  '[atom_types]',
                  '[bonds]',
                  '[angles]',
                  '[torsions]',
                  '[impropers]'
                 ]
        
        prmfiles = [s.FF_DIR + '/' + self.forcefield + '.prm']
        if self.cofactor != None:
            for filename in self.cofactor:
                prmfiles.append(filename + '.prm')
            
        prms = IO.read_prm(prmfiles)
        
        with open (self.forcefield + '_merged.prm', 'w') as outfile:
            for key in headers:
                outfile.write(key + '\n')
                for line in prms[key]:
                    outfile.write(line)
                    
    def get_lambdas(self):
        self.lambdas = IO.get_lambdas(self.windows, self.sampling)
        
    def write_EQ(self):
        for line in self.PDB[1]:
            if line[2] == 'CA' and self.system == 'water' or self.system == 'vacuum':
                self.replacements['WATER_RESTRAINT'] = '{} {} 1.0 0 0'.format(line[1], line[1])
                
        for EQ_file in glob.glob(s.INPUT_DIR + '/eq*.inp'):
            src = EQ_file
            EQ_file = EQ_file.split('/')
            tgt = self.directory + '/inputfiles/' + EQ_file[-1]
            with open(src) as infile, open(tgt, 'w') as outfile:
                for line in infile:
                    outline = IO.replace(line, self.replacements)
                    outfile.write(outline)
                    
    def write_MD(self):
        for i in range(0, int(self.windows) + 1):
            j = int(self.windows) - i
            lambda1 = self.lambdas[i].replace(".", "")
            lambda2 = self.lambdas[j].replace(".", "")
            
            if self.lambdas[i] == '1.000':
                src = s.INPUT_DIR + '/md_1000_0000.inp'
            else:
                src = s.INPUT_DIR + '/md_XXXX_XXXX.inp'
                self.replacements['FLOAT_LAMBDA1'] = self.lambdas[i]
                self.replacements['FLOAT_LAMBDA2'] = self.lambdas[j]
                
            tgt = self.directory + '/inputfiles/md_{}_{}.inp'.format(lambda1,
                                                                     lambda2)
            self.replacements['FILE'] = 'md_{}_{}'.format(lambda1,
                                                          lambda2)
            
            with open(src) as infile, open(tgt, 'w') as outfile:
                for line in infile:
                    outline = IO.replace(line, self.replacements)
                    outfile.write(outline)
            
            ## Store previous file
            self.replacements['FILE_N'] = 'md_{}_{}'.format(lambda1,
                                                            lambda2)
            
    def write_runfile(self):
        ntasks = getattr(s, self.cluster)['NTASKS']
        src = s.INPUT_DIR + '/run.sh'
        tgt = self.directory + '/inputfiles/run_' + self.cluster + '.sh'
        EQ_files = sorted(glob.glob(self.directory + '/inputfiles/eq*.inp'))
        MD_files = reversed(sorted(glob.glob(self.directory + '/inputfiles/md*.inp')))
        replacements = IO.merge_two_dicts(self.replacements, getattr(s, self.cluster))
        replacements['FEPS'] = ' '.join(self.FEPlist)
            
        with open(src) as infile, open(tgt, 'w') as outfile:
            for line in infile:
                if line.strip() == '#SBATCH -A ACCOUNT':
                    try:
                        replacements['ACCOUNT']
                        
                    except:
                        line = ''
                outline = IO.replace(line, replacements)
                outfile.write(outline)
                
                if line.strip() == '#EQ_FILES':
                    for line in EQ_files:
                        file_base = line.split('/')[-1][:-4]
                        outline = 'time mpirun -np {} $qdyn {}.inp > {}.log\n'.format(ntasks,
                                                                                      file_base,
                                                                                      file_base)
                        outfile.write(outline)
                        
                if line.strip() == '#RUN_FILES':
                    for line in MD_files:
                        file_base = line.split('/')[-1][:-4]
                        outline = 'time mpirun -np {} $qdyn {}.inp > {}.log\n'.format(ntasks,
                                                                                      file_base,
                                                                                      file_base)
                        outfile.write(outline)
                        
    def write_submitfile(self):
        IO.write_submitfile(self.directory, self.replacements)
                        
    def write_singletop_FEPs(self):
        # This piece needs to be cleaner at some point
        RES_list = self.PDB[int(self.PDB2Q[self.mutation[1]])]
        RES_dic = {}
        block = 0
        empty = ['[change_charges]\n', 
                 '[FEP]\n',
                 '[atom_types]\n',
                 '[change_atoms]\n',
                 '[softcore]\n',
                 '[bond_types]\n'                
                ]
        atom_map = {}
        for line in RES_list:
            RES_dic[line[2].strip()] = line
            
        for src in sorted(glob.glob(self.FEPdir + '/FEP*.fep')):
            tgt = self.directory + '/inputfiles/' + src.split('/')[-1]
            with open (src) as infile, open (tgt, 'w') as outfile:
                for line in infile:
                    if line == '[atoms]\n':
                        block = 1
                        
                    if line == '[change_bonds]\n':
                        block = 2
                        
                    if line in empty:
                        block = 0
                    
                    if block == 0:
                        outfile.write(line)
                        
                    if block == 1:
                        if line != '[atoms]\n':
                            line = line.split()
                            if len(line) > 1:
                                atnr = RES_dic[line[3]][1]
                                atom_map[line[0]] = atnr

                                try:
                                    comment = line[4]

                                except:
                                    comment = ' '
                                outline = ('{:7s}{:<7d}{:5s}{:5s}{:5s}\n'.format(line[1],
                                                                              atnr,
                                                                              '!',
                                                                              line[3],
                                                                              comment))
                            
                                outfile.write(outline)
                            else:
                                outfile.write('\n')

                        else:
                            outfile.write(line)
                            
                    if block == 2:
                        if line != '[change_bonds]\n':
                            line = line.split()
                            if len(line) > 1:
                                at_1 = atom_map[line[0]]
                                at_2 = atom_map[line[1]]
                                outline = '{:<7d}{:<7d}{:5s}{:5s}\n'.format(at_1,
                                                                         at_2,
                                                                         line[2],
                                                                         line[3]) 
            
                                outfile.write(outline)
                            else:
                                outfile.write('\n')
                
                        else:
                            outfile.write(line)
                            
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
    
    parser.add_argument('-s', '--sampling',
                        dest = "sampling",
                        required = False,
                        choices = ['linear',
                                   'sigmoidal', 
                                   'exponential', 
                                   'reverse_exponential'],
                        help = "Sampling type to be used, default taken from settings.py.")
    
    parser.add_argument('-w', '--windows',
                        dest = "windows",
                        required = False,
                        help = "Amount of windows to be used per FEP window, default taken from settings.py")    
    
    parser.add_argument('-S', '--system',
                        dest = "system",
                        required = True,
                        choices = ['protein', 'water', 'vacuum'],
                        help = "System type, can be protein, water or vacuum")
    
    parser.add_argument('-T', '--temperature',
                        dest = "temperature",
                        required = False,
                        help = "List of temperatures, default taken from settings.py")
    
    parser.add_argument('-r', '--replicates',
                        dest = "replicates",
                        required = False,
                        help = "List of temperatures, default taken from settings.py")   
    
    parser.add_argument('-C', '--cluster',
                        dest = "cluster",
                        required = True,
                        help = "Cluster to use.")   

    args = parser.parse_args()
    run = Run(mutation = args.mutation,
              cofactor = args.cofactor,
              forcefield = args.forcefield,
              sampling = args.sampling,
              windows = args.windows,
              system = args.system,
              cluster = args.cluster,
              temperature = args.temperature,
              replicates = args.replicates,
              include = ('ATOM','HETATM')
             )
    
    run.checkFEP()
    run.create_environment()
    run.readpdb()
    run.read_input()
    run.merge_prm()
    run.get_lambdas()
    run.write_EQ()
    run.write_MD()
    run.write_submitfile()
    run.write_runfile()
    run.write_singletop_FEPs()