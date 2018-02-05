import argparse
import re
import glob
import os
import shutil
import stat
import shlex
from subprocess import check_output

import functions as f
import settings as s

class Run(object):
    """
    Create FEP files from a common substructure for a given set of
    ligands
    """
    def __init__(self, lig1, lig2, FF, system, cluster, sphereradius, cysbond, *args, **kwargs):
        """
        The init method is a kind of constructor, called when an instance
        of the class is created. The method serves to initialize what you
        want to do with the object.
        """
        self.lig1 = lig1
        self.lig2 = lig2
        self.FF = FF
        self.system = system
        self.rootdir = os.getcwd()
        self.cluster = cluster
        self.sphereradius = sphereradius
        self.cysbond = cysbond
        
        if self.system == 'protein':
            # Get last atom and residue from complexfile!
            with open('protein.pdb') as infile:
                for line in infile:
                    try:
                        resnr = int(line[22:26])
                        atnr = int(line[6:11])
                        
                    except:
                        continue

                self.residueoffset = resnr
                self.atomoffset = atnr
                
        else:
            self.atomoffset = 0
            self.residueoffset = 0
        
    def makedir(self):
        lignames = self.lig1 + '-' +self.lig2
        directory = self.rootdir + '/FEP_' + lignames
        if not os.path.exists(directory):
            os.makedirs(directory)
            
        if not os.path.exists(directory + '/inputfiles'):
            os.makedirs(directory + '/inputfiles')

        return directory
    def replace(self, string, replacements):
        pattern = re.compile(r'\b(' + '|'.join(replacements.keys()) + r')\b')
        replaced_string = pattern.sub(lambda x: replacements[x.group()], string)
        return replaced_string
        
    def read_files(self):
        changes_1 = {}
        changes_2 = {}
        charges = []
        atomtypes = []
        merged_molsize = 0
        
        with open(self.lig1 + '.lib') as infile:
            block = 0
            for line in infile:
                line = line.split()
                if len(line) > 0:
                    if line[0] == '[atoms]':
                        block = 1
                        continue
                    if line[0] == '[bonds]':
                        block = 2
                        
                if block == 1 and len(line) > 0:
                    #construct for FEP file
                    merged_molsize = merged_molsize + 1
                    charges.append([merged_molsize, line[3], '0.000'])
                    atomtypes.append([merged_molsize, line[2], 'DUM'])
                
                if block == 2:
                    break
                    
            molsize_lig1 = len(atomtypes)
            
        with open(self.lig2 + '.lib') as infile:
            block = 0
            for line in infile:
                line = line.split()
                if len(line) > 0:
                    if line[0] == '[atoms]':
                        block = 1
                        continue
                    if line[0] == '[bonds]':
                        block = 2
                        
                if block == 1 and len(line) > 0:
                    #construct for FEP file
                    merged_molsize = merged_molsize + 1
                    charges.append([merged_molsize, '0.000', line[3]])
                    
                    
                    #adjustments to be made for lib and prm files
                    cnt = 0
                    for i in [line[1], line[2]]:
                        cnt = cnt + 1
                        match = re.match(r"([a-z]+)([0-9]+)", i, re.I)
                        if match:
                            items = match.groups()
                            j = str(items[0]) + str(int(items[1]) + int(molsize_lig1))
                            
                        if cnt == 1:
                            changes_1[i] = j
                        if cnt == 2:
                            changes_2[i] = j
                            atomtypes.append([merged_molsize ,'DUM', j])
        
        molsize_lig2 = merged_molsize - molsize_lig1
        return [changes_1, changes_2],                      \
               [charges, atomtypes],                        \
               [molsize_lig1, molsize_lig2]
        

        
    def change_lib(self, replacements, writedir):
        replacements['LIG'] = 'LID'
        pattern = re.compile(r'\b(' + '|'.join(replacements.keys()) + r')\b')
        
        with open(self.lig2 + '.lib', 'r') as infile:
            file_replaced = []
            for line in infile:
                line2 = pattern.sub(lambda x: replacements[x.group()], line)
                file_replaced.append(line2)
                
        with open(writedir + '/' + self.lig2 +'_renumber.lib', 'w') as outfile:
            block = 0
            # temporary fix for charge group switching atoms
            for line in file_replaced:
                outfile.write(line)
                
        shutil.copy(self.lig1 + '.lib', writedir + '/' + self.lig1 + '.lib')
                                
            
    def change_prm(self, replacements, writedir):      
        pattern = re.compile(r'\b(' + '|'.join(replacements.keys()) + r')\b')
        file1 = glob.glob(self.lig1 + '*.prm')[0]
        file2 = glob.glob(self.lig2 + '*.prm')[0]
        prm_file = s.FF_DIR + '/' + self.FF + '.prm'
        prm_merged = {'vdw':[],
                       'bonds':[],
                       'angle':[],
                       'torsion':[],
                       'improper':[]
                       }
        
        for file in [file1, file2]:
            with open(file) as infile:
                block = 0
                for line in infile:
                    if file == file2:
                        line = pattern.sub(lambda x: replacements[x.group()], line)
                    if line == '[atom_types]\n':
                        block = 1
                        continue
                    elif line == '[bonds]\n':
                        block = 2
                        continue
                    elif line == '[angles]\n':
                        block = 3
                        continue
                    elif line == '[torsions]\n':
                        block = 4
                        continue
                    if line == '[impropers]\n':
                        block = 5
                        continue

                    if block == 1:
                        prm_merged['vdw'].append(line)

                    elif block == 2:
                        prm_merged['bonds'].append(line)                

                    elif block == 3:
                        prm_merged['angle'].append(line)                

                    elif block == 4:
                        prm_merged['torsion'].append(line)                

                    elif block == 5:
                        prm_merged['improper'].append(line)

        with open(prm_file) as infile, open(writedir + 
                                            '/' + 
                                            self.FF + 
                                            '_' + 
                                            self.lig1 +  
                                            '_' + 
                                            self.lig2 
                                            +'_merged.prm', 'w') as outfile:
            for line in infile:
                block = 0
                outfile.write(line)
                if len(line) > 1:
                    if line == "! Ligand vdW parameters\n":
                        block = 1
                    if line == "! Ligand bond parameters\n":
                        block = 2     
                    if line == "! Ligand angle parameters\n":
                        block = 3
                    if line == "! Ligand torsion parameters\n":
                        block = 4
                    if line == "! Ligand improper parameters\n":
                        block = 5

                # Create lists (do something with containers?) FIX dic keys
                if block == 1: 
                    for line in prm_merged['vdw']:
                        outfile.write(line)

                if block == 2:
                    for line in prm_merged['bonds']:
                        outfile.write(line)

                if block == 3:
                    for line in prm_merged['angle']:
                        outfile.write(line)

                if block == 4:
                    for line in prm_merged['torsion']:
                        outfile.write(line)

                if block == 5:
                    for line in prm_merged['improper']:
                        outfile.write(line)                  
        
        #AND return the vdW list for the FEP file
        FEP_vdw = []
        for line in prm_merged['vdw']:
            if len(line) > 1 and line[0] != '!':
                line = line.split()
                line2 = "{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}".format(line[0],
                                                                      line[1],
                                                                      line[3],
                                                                      str(0),
                                                                      str(0),
                                                                      line[4],
                                                                      line[5],
                                                                      line[6]
                                                                      )
                FEP_vdw.append(line2)
        return FEP_vdw
        
    def write_FEP_file(self, change_charges, change_vdw, FEP_vdw, writedir):
        with open(writedir + '/FEP1.fep', 'w') as outfile:
            total_atoms = len(change_charges)
            outfile.write('!info: ' + self.lig1 + ' --> ' + self.lig2 + '\n')
            outfile.write('[FEP]\n')
            outfile.write('states 2\n\n')
            
            # defining the atom order taken user given offset into account
            outfile.write('[atoms]\n')
            for i in range(1, total_atoms + 1):
                outfile.write("{:5}{:5}\n".format(str(i),
                                                  str(i + self.atomoffset)))
            outfile.write('\n\n')
            
            # changing charges
            outfile.write('[change_charges]\n')

            for line in change_charges:
                outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                        line[1],
                                                        line[2]))
            outfile.write('\n\n')
            
            #add the Q atomtypes
            outfile.write('[atom_types]\n')
            for line in FEP_vdw:
                outfile.write(line + '\n')
                
            outfile.write('DUM       0.0000    0.0000    0         0         0.0000    0.0000    1.0080')
            outfile.write('\n\n')
            
            # changing atom types
            outfile.write('[change_atoms]\n')
            for line in change_vdw:
                outfile.write('{:<5}{:>10}{:>10}\n'.format(line[0],
                                                        line[1],
                                                        line[2]))
            
            
    def merge_pdbs(self, writedir):
        replacements = {}
        replacements['LIG'] = 'LID'
        pattern = re.compile(r'\b(' + '|'.join(replacements.keys()) + r')\b')
        file_replaced = []
        atnr = self.atomoffset
        with open(self.lig2 + '.pdb') as infile:
            for line in infile:
                if line.split()[0] == 'ATOM':
                    line2 = pattern.sub(lambda x: replacements[x.group()], line)
                    file_replaced.append(line2)
                    
        with open(self.lig1 + '.pdb') as infile,                                            \
             open(writedir + '/' + self.lig1 + '_' + self.lig2 + '.pdb', 'w') as outfile:

            if self.system == 'protein':
                with open('protein.pdb') as protfile:
                    for line in protfile:
                        outfile.write(line)
                        
            for line in infile:
                if line.split()[0] == 'ATOM':
                    resnr = int(line[22:26])
                    # Temporary fix NaN error of overlapping heavy atom in Q, add offset
                    atnr += 1
                    atom1 = line.split()
                    atom1[1] = '{:>5}'.format(int(atom1[1]) + self.atomoffset)
                    atom1[4] = '{:>6}'.format(int(atom1[4]) + self.residueoffset)
                    atom1[5] = float(atom1[5]) + 0.001
                    atom1[6] = float(atom1[6]) + 0.001
                    atom1[7] = float(atom1[7]) + 0.001
                    # FIX string formatting for pdb files
                    line = '{:6}{:>5}  {:<4}{:<3}{:>6}{:>12}{:>8}{:>8}\n'.format(*atom1)
              
                    outfile.write(line)
                    
            self.residueoffset = self.residueoffset + 2
            resnr = '{:4}'.format(self.residueoffset)
            for line in file_replaced:
                atnr = atnr + 1
                atchange = '{:5}'.format(atnr)
                line = line[0:6] + atchange + line[11:22] + resnr + line[26:]
                outfile.write(line)
    
    def write_water_pdb(self, writedir):
        # This function could later include a water removal part
        header = self.sphereradius + '.0 SPHERE\n'
        with open('water.pdb') as infile, open(writedir + '/water.pdb', 'w') as outfile:
            outfile.write(header)
            for line in infile:
                outfile.write(line)
        
        
    def get_lambdas(self, windows, sampling):
        windows = int(windows)
        step = windows/2
        lambdas = []
        lmbda_1 = []
        lmbda_2 = []
        k_dic = {'sigmoidal':-1.1, 
                 'linear':1000,
                 'exponential':-1.1,
                 'reverse_exponential':1.1
                }
        k = k_dic[sampling]

        if sampling == 'sigmoidal': 
            for i in range(0, step + 1):
                lmbda1 = '{:.3f}'.format(0.5 * (f.sigmoid(float(i)/float(step), k) + 1))
                lmbda2 = '{:.3f}'.format(0.5 * (-f.sigmoid(float(i)/float(step), k) + 1))
                lmbda_1.append(lmbda1)
                lmbda_2.append(lmbda2)

            lmbda_2 = lmbda_2[1:]

            for i in reversed(lmbda_2):
                lambdas.append(i)

            for i in lmbda_1:
                lambdas.append(i)

        else:
            for i in range(0, windows + 1):
                lmbda = '{:.3f}'.format(f.sigmoid(float(i)/float(windows), k))
                lambdas.append(lmbda)
        
        lambdas = lambdas[::-1]
        return lambdas            
    
    def write_MD(self, lambdas, writedir, lig_size1, lig_size2):
        lambda_tot = len(lambdas)
        replacements = {}
        file_list1 = []
        file_list2 = []
        file_list3 = []
        lig_total = lig_size1 + lig_size2
        lambda_1 = []
        lambda_2 = []
        block = 0
        index = 0
        
        for line in lambdas:
            if line == '0.500':
                block = 1
                
            if block == 0:
                lambda_1.append(line)
                
            if block == 1:
                lambda_2.append(line)

        lambda_1 = lambda_1[::-1]
        lambda_2 = lambda_2[1:]
        replacements['ATOM_START_LIG1'] =   '{:<6}'.format(self.atomoffset + 1)
        replacements['ATOM_END_LIG1']   =   '{:<7}'.format(self.atomoffset + lig_size1)           
        replacements['ATOM_START_LIG2'] =   '{:<6}'.format(self.atomoffset + lig_size1 + 1)
        replacements['ATOM_END_LIG2']   =   '{:<7}'.format(self.atomoffset + lig_size1 + lig_size2)
        replacements['SPHERE']          =   self.sphereradius
        replacements['ATOM_END']        =   '{:<6}'.format(self.atomoffset + lig_total)        
        
        if self.system == 'water' or self.system == 'vacuum':
            replacements['WATER_RESTRAINT'] = '{:<7}{:<7} 1.0 0 1   '.format(self.atomoffset + 1, 
                                                                          self.atomoffset + lig_size1 + 
                                                                          lig_size2
                                                                         )
            
            
        elif self.system == 'protein':
            replacements['WATER_RESTRAINT'] = ''
        
        for eq_file_in in glob.glob(s.ROOT_DIR + '/INPUTS/eq*.inp'):
            eq_file = eq_file_in.split('/')[-1:][0]
            eq_file_out = writedir + '/' + eq_file

            with open(eq_file_in) as infile, open(eq_file_out, 'w') as outfile:
                for line in infile:
                    line = run.replace(line, replacements)
                    outfile.write(line)

                file_list1.append(eq_file)
                
        file_in = s.INPUT_DIR + '/md_0500_0500.inp'
        file_out = writedir + '/md_0500_0500.inp' 
        with open(file_in) as infile, open(file_out, 'w') as outfile:
            for line in infile:
                line = run.replace(line, replacements)
                outfile.write(line)
        file_list1.append('md_0500_0500.inp')
        
        for lambdas in [lambda_1, lambda_2]:
            index += 1
            filename_N = 'md_0500_0500'
            filenr = -1

            for line in lambdas:
                filenr += 1
                if index == 1:
                    lambda1 = lambda_1[filenr]
                    lambda2 = lambda_2[filenr]
                    
                elif index == 2:
                    lambda1 = lambda_2[filenr]
                    lambda2 = lambda_1[filenr]
                    
                filename = 'md_' + lambda1.replace('.', '') + '_' + lambda2.replace('.', '')
                replacements['FLOAT_LAMBDA1']   =   lambda1 
                replacements['FLOAT_LAMBDA2']   =   lambda2
                replacements['FILE']          =   filename
                replacements['FILE_N'] = filename_N

                # Consider putting this in a function seeing as it is called multiple times
                pattern = re.compile(r'\b(' + '|'.join(replacements.keys()) + r')\b')
                file_in = s.INPUT_DIR + '/md_XXXX_XXXX.inp'
                file_out = writedir + '/' + filename + '.inp'

                with open(file_in) as infile, open(file_out, 'w') as outfile:
                    for line in infile:
                        line = pattern.sub(lambda x: replacements[x.group()], line)
                        outfile.write(line)

                filename_N = filename
                
                if index == 1:
                    file_list2.append(filename + '.inp')
                    
                elif index == 2:
                    file_list3.append(filename + '.inp')
        return [file_list1, file_list2, file_list3]
    
    def write_submitfile(self, writedir):
        replacements = {}
        replacements['TEMP_VAR']    = s.TEMPERATURE
        replacements['RUN_VAR']     = s.REPLICATES
        replacements['RUNFILE']     = 'run' + self.cluster + '.sh'
        submit_in = s.ROOT_DIR + '/INPUTS/FEP_submit.sh'
        submit_out = writedir + ('/FEP_submit.sh')
        with open(submit_in) as infile, open (submit_out, 'w') as outfile:
            for line in infile:
                line = run.replace(line, replacements)
                outfile.write(line)
        
        try:
            st = os.stat(submit_out)
            os.chmod(submit_out, st.st_mode | stat.S_IEXEC)
        
        except:
            print "WARNING: Could not change permission for " + submit_out

        
    def write_runfile(self, writedir, file_list):
        cluster = self.cluster
        run_file_in = s.ROOT_DIR + '/INPUTS/run.sh' 
        run_file_out = writedir + '/run' + cluster + '.sh'
        run1_in = s.ROOT_DIR + '/INPUTS/run_0500-1000.sh'
        run1_out = writedir + '/run_0500-1000.sh'
        run2_in = s.ROOT_DIR + '/INPUTS/run_0500-0000.sh'
        run2_out = writedir + '/run_0500-0000.sh'
        replacements = getattr(s, cluster)
        run_threads = '{}'.format(int(replacements['NTASKS']))
        
        with open(run_file_in) as infile, open(run_file_out, 'w') as outfile:
            block = 0
            for line in infile:
                line = run.replace(line, replacements)
                outfile.write(line)

                if line == '#EQ_FILES\n':
                    for line in file_list[0]:
                        outfile.write('time mpirun -np ' + 
                                      replacements['NTASKS'] + 
                                      ' $qdyn ' +
                                      line + 
                                      ' > ' + 
                                      line[:-4] +
                                      '.log\n')
                
                i = -1
                if line == '#RUN_FILES\n':
                #    for line in file_list[1]:
                #        i += 1
                #        
                #        runline1 = 'time mpirun -np {} $qdyn {} > {}.log &\n'.format(run_threads,
                #                                                                    line,
                #                                                                    line[:-4],)
                #        
                #        runline2 = 'time mpirun -np {} $qdyn {} > {}.log &\n'.format(run_threads,
                #                                                                     file_list[2][i],
                #                                                                     file_list[2][i][:-4])
                #                   
                #        outfile.write(runline1)    
                #        outfile.write(runline2)
                #        outfile.write('wait \n')
                
                    for line in file_list[1]:
                        runline = 'time mpirun -np {} $qdyn {} > {}.log \n'.format(run_threads,
                                                                                     line,
                                                                                     line[:-4],)
                        outfile.write(runline)

                    for line in file_list[2]:
                        runline = 'time mpirun -np {} $qdyn {} > {}.log \n'.format(run_threads,
                                                                                     line,
                                                                                     line[:-4],)
                        outfile.write(runline)
                            
    def write_qprep(self, writedir):
        replacements = {}
        # QUICK fix for protein file, should be COG of lig1/lig2 not lig1 alone
        center = f.COG(self.lig1 + '.pdb')
        center = '{:} {:} {:}'.format(center[0], center[1], center[2])
        qprep_in = s.ROOT_DIR + '/INPUTS/qprep.inp'
        qprep_out = writedir + '/qprep.inp'
        replacements['FF_LIB'] = s.ROOT_DIR + '/FF/' + self.FF + '.lib'
        replacements['LIG1']   = self.lig1 + '.lib'
        replacements['LIG2']   = self.lig2 + '_renumber.lib'
        replacements['LIGPRM'] = self.FF + '_' + self.lig1 + '_' + self.lig2 + '_merged.prm'
        replacements['LIGPDB'] = self.lig1 + '_' + self.lig2 + '.pdb'
        replacements['CENTER'] = center
        replacements['SPHERE'] = self.sphereradius
        if self.system =='vacuum':
            replacements['solvate'] = '!solvate'
        if self.system == 'water':
            replacements['SOLVENT']  = '1 HOH'
        if self.system == 'protein':
            replacements['SOLVENT']  = '4 water.pdb'            
        
        with open(qprep_in) as infile, open(qprep_out, 'w') as outfile:
            for line in infile:
                line = run.replace(line, replacements)
                if line == '!addbond at1 at2 y\n' and self.cysbond != None:
                    cysbond = self.cysbond.split(',')
                    for cys in cysbond:
                        at1 = cys.split(':')[0]
                        at2 = cys.split(':')[1]
                        outfile.write('addbond ' + at1 + ' ' + at2 + '\n' )
                    continue
                outfile.write(line)
        
    def qprep(self, writedir):
        os.chdir(writedir)
        cluster_options = getattr(s, self.cluster)
        qprep = cluster_options['QPREP']
        args = shlex.split(qprep + ' < qprep.inp' )
        out = check_output(args)
        with open('qprep.out', 'w') as outfile:
            for line in out:
                outfile.write(line)
        os.chdir('../../')
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='ligFEP',
        version='1.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Generate FEP files for dual topology ligFEP == ')

    
    parser.add_argument('-l1', '--lig_1',
                        dest = "lig1",
                        required = True,
                        help = "name of ligand 1")
    
    parser.add_argument('-l2', '--lig_2',
                        dest = "lig2",
                        required = True,
                        help = "name of ligand 2")
    
    parser.add_argument('-FF', '--forcefield',
                        dest = "FF",
                        required = True,
                        choices = ['OPLS2005', 'OPLS2015'],
                        help = "Forcefield to be used")
    
    parser.add_argument('-s', '--system',
                        dest = "system",
                        required = True,
                        choices = ['water', 'protein', 'vacuum'],
                        help = "what type of system we are setting up")
    
    parser.add_argument('-c', '--cluster',
                        dest = "cluster",
                        required = True,
                        help = "cluster you want to submit to, cluster specific parameters added to settings"
                       )
    
    parser.add_argument('-r', '--sphereradius',
                        dest = "sphereradius",
                        required = False,
                        default = '15',
                        help = "size of the simulation sphere"
                       ) 
    
    parser.add_argument('-b', '--cysbond',
                        dest = "cysbond",
                        default = None,
                        help = "Temporary function to add cysbonds at1:at2,at3:at4 etc."
                       )    
    
    args = parser.parse_args()
    run = Run(lig1 = args.lig1,
              lig2 = args.lig2,
              FF= args.FF,
              system = args.system,
              cluster = args.cluster,
              sphereradius = args.sphereradius,
              cysbond = args.cysbond
             )

    writedir = run.makedir()
    inputdir = writedir + '/inputfiles'
    a = run.read_files()
    changes_for_libfiles = a[0][1]
    changes_for_prmfiles = a[0][1]
    change_charges       = a[1][0]
    change_vdw           = a[1][1]
    changes_for_pdbfiles = a[0][0]
    lig_size1, lig_size2 = a[2][0], a[2][1]

    # Write the merged files
    run.change_lib(changes_for_libfiles, inputdir)
    FEP_vdw = run.change_prm(changes_for_prmfiles, inputdir)
    run.write_FEP_file(change_charges, change_vdw, FEP_vdw, inputdir)
    run.merge_pdbs(inputdir)
    run.write_water_pdb(inputdir)
    lambdas = run.get_lambdas(s.WINDOWS, s.SAMPLING)
    file_list = run.write_MD(lambdas, inputdir, lig_size1, lig_size2)
    run.write_submitfile(writedir)
    run.write_runfile(inputdir, file_list)
    run.write_qprep(inputdir)
    run.qprep(inputdir)
    
