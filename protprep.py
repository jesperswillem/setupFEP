import argparse
from subprocess import check_output
import os
import collections

import functions as f
import settings as s
import IO

class Run(object):
    """
    Create FEP files from a common substructure for a given set of
    ligands
    """
    def __init__(self, prot, sphereradius, spherecenter, include, water, *args, **kwargs):
        """
        The init method is a kind of constructor, called when an instance
        of the class is created. The method serves to initialize what you
        want to do with the object.
        """
        self.prot = prot
        self.radius = float(sphereradius)
        self.center = spherecenter
        self.include = include
        self.water = water
        self.PDB = {}
        # Log file list to keep track of stuff to write out
        self.log = {'CENTER':None,
                    'DECHARGE':[],
                    'CTERM':[],
                    'CYX':[],
                    'QRESN':{},
                    'QRES_LIST':[]
                   }
        
    def get_center_coordinates(self):
        center = self.center
        center = center.split(':')
        with open(self.prot) as infile:
            for line in infile:
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
                    if center[0] == 'RESN':
                        if line[6] == int(center[1]) \
                        and line[2].strip() == 'CA':
                            self.center = [float(line[8]),
                                           float(line[9]), 
                                           float(line[10])
                                          ]

                    elif center[0] == 'ATN':
                        if line[1] == int(center[1]):
                            self.center = [float(line[8]),
                                           float(line[9]),
                                           float(line[10])
                                          ]
                    
                    elif len(center) == 3:
                        self.center = [float(center[0]),
                                       float(center[1]),
                                       float(center[2])
                                      ]
                        
                    else:
                        print 'Could not get center'
                        
        self.log['CENTER'] = '{} {} {}'.format(*self.center)
    
    def readpdb(self):
        i = 0
        with open(self.prot) as infile:
            header = IO.pdb_parse_in(infile.readline())
            RES_ref = header[6] - 1
            
        with open(self.prot) as infile:
            for line in infile:
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
                    RES = line[6]
                    self.PDB[line[1]] = line
                    if RES != RES_ref:
                        RES_ref = RES
                        i += 1
                        self.log['QRESN'][line[6]] = i
                        if line[4].strip() != 'HOH':
                            self.log['QRES_LIST'].append('{:<10d}{:<10d}{:<10}'.format(i, 
                                                                                  int([line[6]][0]),
                                                                                  line[4]
                                                                                 ))
    
    def decharge(self):
        charged_res = {'GLU':['GLH', 'CD'], 
                       'ASP':['ASH', 'CG'], 
                       'ARG':['ARN', 'CZ'], 
                       'LYS':['LYN', 'NZ'], 
                       'HIP':['HID', 'CG']
                      }
        
        coord1 = self.center
        decharge = []
        # Distance for decharging residues in boundary
        # TO DO, remove charges within salt bridge pair in boundary(?)
        rest_bound = float(self.radius) - 3.0
        for key in self.PDB:
            at = self.PDB[key]
            if at[4] in charged_res:
                if at[2].strip() == charged_res[at[4]][1]:
                    coord2 = [float(at[8]), 
                              float(at[9]), 
                              float(at[10])
                             ]
                    if f.euclidian_overlap(coord1, coord2, rest_bound) == False:
                        decharge.append(at[6])
                        self.log['DECHARGE'].append('{} {}'.format(at[6], 
                                                                   at[4]
                                                                  ))
                        
        # Check if the decharged residue is part of a salt bridge and
        # neutralize this residue as well
        for key in self.PDB:
            at = self.PDB[key]
            if at[6] in decharge and at[2].strip() == charged_res[at[4]][1]:
                coord1 = [float(at[8]), 
                          float(at[9]), 
                          float(at[10])
                         ]
                for key in self.PDB:
                    at_2 = self.PDB[key]
                    if at_2[4] in charged_res:
                        if at_2[2].strip() == charged_res[at_2[4]][1]:
                            coord2 = [float(at_2[8]), 
                                      float(at_2[9]), 
                                      float(at_2[10])
                                     ]
                            if at != at_2 and at_2[6] not in decharge:
                                if f.euclidian_overlap(coord1, coord2, 4.0) == True:
                                    decharge.append(at_2[6])
                                    self.log['DECHARGE'].append('{} {}'.format(at_2[6], 
                                                                               at_2[4]
                                                                              ))

        
        for key in self.PDB:
            at = self.PDB[key]
            if self.PDB[key][6] in decharge:
                self.PDB[key][4] = charged_res[at[4]][0]
                
            else:
                continue
            
    def set_OXT(self):
        CTERM = []
        for key in self.PDB:
            at = self.PDB[key]
            if at[2].strip() == 'OXT':
                CTERM.append(at[6])
                self.log['CTERM'].append('{} {}'.format(at[6], at[4]))
                
        for key in self.PDB:
            at = self.PDB[key]
            if at[6] in CTERM:
                self.PDB[key][4] = 'C' + at[4]
                
    def get_CYX(self):
        cys = {}
        cyx = []
        cys_bond = 2.2
        
        # Reduce coordinate matrix
        for key in self.PDB:
            at = self.PDB[key]
            if at[4] == 'CYS' and at[2].strip() == 'SG':
                cys[at[6]] = [at[8], at[9], at[10]]
        
        # Find S-S bonds
        for SG_1 in cys:
            for SG_2 in cys:
                if SG_1 != SG_2:
                    if SG_1 not in cyx or SG_2 not in cyx:
                        if f.euclidian_overlap(cys[SG_1], cys[SG_2], cys_bond) == True:
                            cyx.append(int(SG_1))
                            SG_1_Q = self.log['QRESN'][SG_1]
                            SG_2_Q = self.log['QRESN'][SG_2]
                            self.log['CYX'].append('{}:SG {}:SG'.format(SG_1_Q, SG_2_Q))

        for key in self.PDB:
            at = self.PDB[key]
            if at[6] in cyx:
                self.PDB[key][4] = 'CYX'
                
    def write_tmpPDB(self):
        with open(self.prot[:-4] + '_tmp.pdb', 'w') as outfile:
            for key in self.PDB:
                outline = IO.pdb_parse_out(self.PDB[key]) + '\n'
                outfile.write(outline)
        
    def write_qprep(self):
        replacements = {'FF_LIB'    :   s.FF_DIR + '/OPLS2015.lib',
                        'FF_PRM'    :   s.FF_DIR + '/OPLS2015.prm',
                        'PROTPDB'   :   self.prot[:-4] + '_tmp.pdb',
                        'CENTER'    :   self.log['CENTER'],
                        'SPHERE'    :   '{:.1f}'.format(self.radius),
                        'SOLVENT'   :   '1 HOH'
                       }
        
        with open (s.INPUT_DIR + '/qprep_protprep.inp') as infile, \
            open ('qprep.inp', 'w') as outfile:
            for line in infile:
                line = IO.replace(line, replacements)
                outfile.write(line)
                if line[0:8] == '!addbond':
                    for line in self.log['CYX']:
                        outline = 'addbond {} y\n'.format(line)
                        outfile.write(outline)
            
        
    def run_qprep(self):
        qprep = s.Q_DIR['LOCAL'] + 'qprep'
        options = ' < qprep.inp > qprep.out'
        # Somehow Q is very annoying with this < > input style so had to implement
        # another function that just calls os.system instead of using the preferred
        # subprocess module....
        IO.run_command(qprep, options, string = True)
        
    ##### THIS PART COMES AFTER THE TEMP .pdb IS WRITTEN ######
    def write_pdb_out(self):
        waters ={'HOH': ['O', 'H1', 'H2'],
                 'SOL': ['OW1', 'HW1', 'HW2'] 
                }
        waters_tokeep = []
        
        if self.water != True:
            return None
        
        with open('top_p.pdb') as infile:
            for line in infile:
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
                    if line[4].strip() in waters and \
                       line[2].strip() == waters[line[4].strip()][0]:
                        coord1 = self.center
                        coord2 = [float(line[8]), 
                                  float(line[9]), 
                                  float(line[10])
                                 ]
                        if f.euclidian_overlap(coord1, coord2, self.radius) == True:
                            waters_tokeep.append(line[6])
                            
        with open('top_p.pdb') as infile, \
             open('water.pdb', 'w') as watout, \
             open('protein.pdb', 'w') as protout:
                    
            watout.write('{:<7.1f}SPHERE\n'.format(self.radius))
            for line in infile:
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
                    if line[6] in waters_tokeep:
                        outline = IO.pdb_parse_out(line) + '\n'
                        watout.write(outline)
                        
                    if line[4] not in waters:
                        outline = IO.pdb_parse_out(line) + '\n'
                        protout.write(outline)                        
                
    def write_log(self):
        with open('protPREP.log', 'w') as outfile:
            outfile.write('The following residues have been decharged:\n')
            for line in self.log['DECHARGE']:
                outfile.write(line + '\n')
                
            outfile.write('\n')
            outfile.write('The following S-S bonds have been found:\n')
            for line in self.log['CYX']:
                outfile.write(line + '\n')
                
            outfile.write('\n')
            outfile.write('The following is a mapping of residue numbers in Q and the input .pdb:\n')
            outfile.write('{:10}{:10}{:10}\n'.format('Q_RESN', 'PDB_IN', 'RESNAME'))
            for line in self.log['QRES_LIST']:
                outfile.write(line + '\n')
                    
    def cleanup(self):
        os.remove(self.prot[:-4] + '_tmp.pdb')
        #os.remove('qprep.inp')
        #os.remove('qprep.out')
        #os.remove('top_p.pdb')
        #os.remove('complexnotexcluded.pdb')
        #os.remove('tmp.top')
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='protPREP',
        version='1.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Generate FEP files for dual topology ligFEP == ')

    
    parser.add_argument('-p', '--prot',
                        dest = "prot",
                        required = True,
                        help = "protein pdb file")
    
    parser.add_argument('-r', '--sphereradius',
                        dest = "sphereradius",
                        required = True,
                        help = "radius of the sphere")    
    
    parser.add_argument('-w', '--nowater',
                        dest = "water",
                        required = False,
                        default = True,
                        action = 'store_false',
                        help = "set no water if crystal waters are NOT to be retained")
    
    parser.add_argument('-c', '--spherecenter',
                        dest = "spherecenter",
                        required = True,
                        help = "center of the sphere, can be residue number (RESN:$)," \
                        "atomnumber (ATN:$) or explicit coordinates (X:Y:Z)")
    
    args = parser.parse_args()
    run = Run(prot = args.prot,
              sphereradius = args.sphereradius,
              spherecenter = args.spherecenter,
              water = args.water,
              include = ('ATOM','HETATM')
             )
    
    run.readpdb()                       # 01
    run.get_center_coordinates()        # 02
    run.decharge()                      # 03
    run.set_OXT()                       # 04
    run.get_CYX()                       # 05
    run.write_tmpPDB()                  # 06
    run.write_qprep()                   # 07
    run.run_qprep()                     # 08
    run.write_pdb_out()                 # 09
    run.write_log()                     # 10
    run.cleanup()                       # 12