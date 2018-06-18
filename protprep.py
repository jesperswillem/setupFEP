import argparse

import functions as f
import settings as s
import IO


class Run(object):
    """
    Create FEP files from a common substructure for a given set of
    ligands
    """
    def __init__(self, prot, sphereradius, spherecenter, include, *args, **kwargs):
        """
        The init method is a kind of constructor, called when an instance
        of the class is created. The method serves to initialize what you
        want to do with the object.
        """
        self.prot = prot
        self.radius = float(sphereradius)
        self.center = spherecenter
        self.include = include
        
        
    def get_center_coordinates(self):
        center = self.center
        center = center.split(':')
        with open(self.prot) as infile:
            for line in infile:
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
                    if center[0] == 'RESN':
                        if line[6] == int(center[1]) \
                        and line[2].strip() == 'CB':
                            self.center = [float(line[8]),
                                           float(line[9]), float(line[10])
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
                        
    def decharge(self):
        charged_res = {'GLU':['GLH', 'CD'], 
                       'ASP':['ASH', 'CG'], 
                       'ARG':['ARN', 'CZ'], 
                       'LYS':['LYN', 'NZ'], 
                       'HIP':['HID', 'CG']
                      }
        
        coord1 = self.center
        decharge = []
        
        with open(self.prot) as infile:
            for line in infile:
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
                    if line[4] in charged_res:
                        if line[2].strip() == charged_res[line[4]][1]:
                            coord2 = [float(line[8]), 
                                      float(line[9]), 
                                      float(line[10])
                                     ]
                            if f.eucledian_overlap(coord1, coord2, self.radius) == False:
                                decharge.append(line[6])

        with open(self.prot) as infile, \
        open('protein_decharge.pdb', 'w') as outfile:
            for line in infile:
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
                    if line[6] in decharge:
                        line[4] = charged_res[line[4]][0]
                    line = IO.pdb_parse_out(line) + '\n'
                    outfile.write(line)
                
        
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
    
    parser.add_argument('-c', '--spherecenter',
                        dest = "spherecenter",
                        required = True,
                        help = "center of the sphere can be residue number (RESN:$), atomnumber (ATN:$) or explicit coordinates (X:Y:Z)")    
    
    
    args = parser.parse_args()
    run = Run(prot = args.prot,
              sphereradius = args.sphereradius,
              spherecenter = args.spherecenter,
              include = ('ATOM','HETATM')
             )
    
    run.get_center_coordinates()
    run.decharge()