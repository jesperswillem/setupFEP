import argparse

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
                    'CYX':[]
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
        with open(self.prot) as infile:
            for line in infile:        
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
                    self.PDB[line[1]] = line
    
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
        # TO DO, pairwise alignment within a salt bridge(?)
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
                
        for SG_1 in cys:
            for SG_2 in cys:
                if SG_1 != SG_2:
                    if SG_1 not in cyx or SG_2 not in cyx:
                        if f.euclidian_overlap(cys[SG_1], cys[SG_2], cys_bond) == True:
                            cyx.append(int(SG_1))
                            self.log['CYX'].append('{}:SG {}:SG'.format(SG_1, SG_2))

        for key in self.PDB:
            at = self.PDB[key]
            if at[6] in cyx:
                self.PDB[key][6] = 'CYX'
        
    def write_qprep(self):
        with open (s.INPUT_DIR + '/qprep_protprep.inp') as infile:
            for line in infile:
                print line
            
        
    def run_qprep(self):
        return None
        
    ##### THIS PART COMES AFTER THE TEMP .pdb IS WRITTEN ######
    def get_water(self):
        waters ={'HOH': ['O1', 'H1', 'H2'],
                 'SOL': ['OW1', 'HW1', 'HW2'] 
                }
        if self.water != True:
            return None
        
        with open(self.prot) as infile:
            for line in infile:
                continue
                
    def print_log(self):
        print self.log
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='ligFEP',
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
    
    run.readpdb()                       # 1
    run.get_center_coordinates()        # 2
    run.decharge()                      # 3
    run.set_OXT()                       # 4
    run.get_water()                     # 5
    run.get_CYX()                       # 6
    run.write_qprep()                   # 7
    run.print_log()                     # 8

    