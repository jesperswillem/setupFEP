import argparse
import shutil
import os
import sys
import string
import math
import time
import glob
import settings as s


class Run(object):
    """
    Create FEP files from a common substructure for a given set of
    ligands
    """
    def __init__(self, lig, *args, **kwargs):
        """
        The init method is a kind of constructor, called when an instance
        of the class is created. The method serves to initialize what you
        want to do with the object.
        """
        self.lig = lig
    
    def make_dir(self):
        """
        Create the main directories for the preparation of the ligands
        """
        dir_make = [str("/" + self.lig),
                    str("/" + self.lig + "/ligand"), 
                    str("/" + self.lig + "/Q_topology"),
                   ]
        for line in dir_make:
            try: 
                os.mkdir(s.ROOT_DIR + line)
            except OSError:
                if not os.path.isdir(s.ROOT_DIR + line):
                    raise

    def get_environment(self):
        ligdir = s.ROOT_DIR + "/" + self.lig + "/ligand"
        q_dir = s.ROOT_DIR + "/" + self.lig + "/Q_topology"
        ff_dir = s.ROOT_DIR + "/FF"

        dir_list = [s.ROOT_DIR, ligdir, q_dir, ff_dir]
        return dir_list    

    def generate_dummy(substructure, superstructure, dir_list, include='ATOM,HETATM'):
        change_atom = []
        substructure_atoms = []
        substructure_dic = {}
        superstructure_atoms = []
        superstructure_dic = {}
        change_atom_dic = {}
        quasi_sub = []
        atcnt = 0
        cnt = 0
        bond_atom = []
        include = tuple(include.split(','))

        #add more later
        heavy_atoms = ['C', 'N', 'O']

        with open(substructure, 'r') as pdb:
            for line in pdb:
                if line.startswith(include):
                    substructure_dic[line[11:16]] = line[0:11]
                    substructure_atoms.append(line)


        with open(superstructure, 'r') as pdb:
            for line in pdb:
                if line.startswith(include):
                    superstructure_dic[line[11:16]] = line[0:11]
                    superstructure_atoms.append(line)

        #First find the atoms that will change to non dummy
        for line in substructure_atoms:
            if line[11:16] not in superstructure_dic:
                change_atom_dic[line[11:16]] = [float(line[30:38]),   # x_coord
                                                float(line[38:46]),   # y_coord
                                                float(line[46:54])]    # z_coord

            else:
                atcnt = atcnt + 1
                quasi_sub.append("HETATM%s%s LIG%s%s%s%s \n"
                                % (str(atcnt).rjust(5, ' '), 
                                   line[11:16].rjust(4, ' '),
                                   str(500).rjust(6, ' '), 
                                   line[30:38].rjust(12, ' '),
                                   line[38:46].rjust( 8, ' '),
                                   line[46:54].rjust( 8, ' ')))

        #Now prepare the quasi substructure atom with dummy atoms etc
        # Dummy atoms
        for line in superstructure_atoms:
            if line[11:16] not in substructure_dic:
                if line[11:16] not in change_atom_dic:
                    atcnt = atcnt + 1
                    quasi_sub.append("HETATM%s  DUM LIG%s%s%s%s \n"
                                    % (str(atcnt).rjust(5, ' '), 
                                       str(500).rjust(6, ' '), 
                                       line[30:38].rjust(12, ' '),
                                       line[38:46].rjust( 8, ' '),
                                       line[46:54].rjust( 8, ' '))) 

            # Find atoms that need to be reconnected
            # This makes the assumption that changing atoms is always heavy atom
            # to H -> potentially dangerous..
            # also smarten up on neighbouring atoms etc....
                elif ((float(line[30:38])-line2[0])**2 + \
                    (float(line[38:46])-line2[1])**2 + \
                    (float(line[46:54])-line2[2])**2)  \
                    < 1.7 **2: 
                    # Change the bonded atoms
                    if line[13] in heavy_atoms:
                        atcnt = atcnt+1
                        quasi_sub.append("HETATM%s%s LIG%s%s%s%s \n"
                                        %   (str(atcnt).rjust(5, ' '), 
                                            line2[3].rjust(4, ' '),
                                            str(500).rjust(6, ' '), 
                                            str(line2[0]).rjust(12, ' '),
                                            str(line2[1]).rjust( 8, ' '),
                                            str(line2[2]).rjust( 8, ' ')))

        with open("test.pdb", 'w') as out_file:
            for line in quasi_sub:
                print line
                out_file.write(line)





    def get_OPLS(self, dir_list):
        """
        Create parameters for the provided ligand pdb's
        """
        os.chdir(dir_list[1])
        shutil.copy2(dir_list[0] + "/" + self.lig + ".pdb", ".")
        os.rename(self.lig + ".pdb", "lig.pdb")
        
        print "========== Preparing ligand =========================="
        print "1) Get OPLS parameters"
        os.system(s.SCHROD_DIR + "/utilities/structconvert " \
                  "-no_renum "                               \
                  "-no_fixelem "                             \
                  "-no_reorder "                             \
                  "-ipdb lig.pdb "                           \
                  "-omae lig.mae "
                 )                                           
                
        # Create the com file
        # Or just use a standard .com file in the FF file
        with open("lig.com", 'w') as outfile:
            outfile.write("lig.mae\n"                                   \
                          "lig-out.mae\n"                               \
                          " MMOD       0      1      0      0     "     \
                          "0.0000     0.0000     0.0000     0.0000 \n"  \
                          " FFLD      14      1      0      0     "     \
                          "1.0000     0.0000     0.0000     0.0000 \n"  \
                          " BDCO       0      0      0      0    "      \
                          "41.5692 99999.0000     0.0000     0.0000 \n" \
                          " READ       0      0      0      0     "     \
                          "0.0000     0.0000     0.0000     0.0000 \n"  \
                          " ELST       1      0      0      0     "     \
                          "0.0000     0.0000     0.0000     0.0000 \n"  \
                          " DEBG      44      0      0      0     "     \
                          "0.0000     0.0000     0.0000     0.0000 \n"
                         )
            
        # Run Bmin minimalisation
        print "Running minimisation"
        os.system(s.SCHROD_DIR + "/bmin -WAIT lig")
        time.sleep(5)
        
    def convert_OPLS(self, dir_list):
        os.chdir(dir_list[1])
        print "Converting OPLS parameters to Q parameters"
        ## This is Marton's code, but adapted for this script
        # init lists
        vdw=[]
        bonds=[]
        angles=[]
        torsions=[]
        torsionsv4=[[]]
        improps=[]
        nonbond=[]
        charges=[]
        elems=dict()
        mass={"H":1.01,"C":12.01,"N":14.01,"O":16.00,"F":19.00,
              "P":30.97,"S":32.07,"Cl":35.45,"Br":79.90,"I":126.90}
        
        def pfloat(s):
            num=s
            if num[0]==".":
                num="0"+num
            if num.count("E"):
                a,b=num.split("E")
                return float(a)*(10.**int(b))
            elif num.count("e"):
                a,b=num.split("e")
                return float(a)*(10.**int(b))
            else:
                return float(num)
        
        block = 0 
        
        with open(dir_list[1] + "/lig-out.mmo", 'r') as infile:
            for line in infile:
                sl=line.split()

                if block==1 and len(sl)>3:
                    if string.digits.count(sl[0][0]):
                        bonds.append([int(sl[0]),
                                      int(sl[1]),
                                      pfloat(sl[2]),
                                      pfloat(sl[3])])
                if block==2 and len(sl)>5:
                    if string.digits.count(sl[0][0]):
                        angles.append([int(sl[0]),
                                       int(sl[1]),
                                       int(sl[2]),
                                       pfloat(sl[3]),
                                       pfloat(sl[5])])
                    if line.count("V4/2") and len(sl)>5:
                        torsionsv4[-1].append(pfloat(sl[3]))
                        
                if block==3 and len(sl)>6:
                        if string.digits.count(sl[0][0]):
                                if len(torsionsv4[-1])<5:
                                        torsionsv4[-1]=[int(sl[0]),
                                                        int(sl[1]),
                                                        int(sl[2]),
                                                        int(sl[3])]
                                else:
                                        torsionsv4.append([int(sl[0]),
                                                           int(sl[1]),
                                                           int(sl[2]),
                                                           int(sl[3])])
                if block==6 and len(sl)>7:
                    if string.digits.count(sl[2][0]):
                        elem=sl[0][0]
                        if sl[0][1]=="l" or sl[0][1]=="r":
                            elem+=sl[0][1]
                        num=int(sl[2].translate(None, ")"))
                        elems[num]=elem
                        charges.append([num,
                                        elem,
                                        pfloat(sl[3]),
                                        pfloat(sl[5]),
                                        pfloat(sl[6]),
                                        pfloat(sl[7])
                                       ])
                if line.count("BOND LENGTHS AND STRETCH ENERGIES"):
                    block=1
                if line.count("ANGLES, BEND AND STRETCH BEND ENERGIES"):
                    block=2
                if line.count("DIHEDRAL ANGLES AND TORSIONAL ENERGIES"):
                    block=3
                if line.count("Improper torsion interaction"):
                    block=4
                if line.count("NONBONDED DISTANCES AND ENERGIES"):
                    block=5
                if line.count("Atomic Charges, Coordinates and Connectivity"):
                    block=6
                if line.count("TIME"):
                    block=0

        if len(torsionsv4[-1])==4:
            del torsionsv4[-1]
        
        # reading log: torsions, impropers, nonbonded
        with open("lig.log", 'r') as infile:
            block = 0
            params=[]
            for line in infile:
                if line.count("quality"):
                    if block==3:
                        v4=0.
                        if len(torsionsv4)>0:
                            for i in torsionsv4:
                                if (at1+at2+at3+at4 == 
                                     str(i[0])+
                                     str(i[1])+
                                     str(i[2])+
                                     str(i[3])):
                                        v4=i[4]

                        torsions.append([int(at1),
                                         int(at2),
                                         int(at3),
                                         int(at4),
                                         pfloat(params[1]),
                                         pfloat(params[2]),
                                         pfloat(params[3]),
                                         v4
                                        ])
                    if block==4:
                        improps.append([int(at1),
                                        int(at2),
                                        int(at3),
                                        int(at4),
                                        pfloat(params[1])
                                       ])
                    if block==5:
                        nonbond.append([int(at1),
                                        int(at2),
                                        pfloat(params[0]),
                                        pfloat(params[1])
                                       ])

                    block = 0
                    params = []

                if block==3:
                    if line.count("atoms"):
                        head,at1,at2,at3,at4=line.split()
                    elif line.count("params"):
                        pass
                    else:
                        params.append(line.strip())
                if block==4:
                    if line.count("atoms"):
                        head,at1,at2,at3,at4=line.split()
                    elif line.count("params"):
                        pass
                    else:
                        params.append(line.strip())
                if block==5:
                    sl=line.split()
                    params.append(sl[1])
                    params.append(sl[2])
                if block==6:
                    sl=line.split()
                    if len(sl)>4:
                        charges[int(sl[0])-1][2]=pfloat(sl[1])
                if line.count("--- torsion"):
                    block=3
                if line.count("--- improper torsion"):
                    block=4
                if line.count("bopls_mmshare_getvdw"):
                    block=5
                    sl=line.split()
                    at1=sl[2]
                    at2=sl[3]
                if line.count("Partial"):
                    block=6


    # guessing vdw parameters from opls force field
        opls = []
        with open(dir_list[3]+ "/f14_oplsaa.vdw",'r') as infile:
            atomtypes=[[] for i in range(len(charges))]
            for line in infile:
                if not line[0]=="#":
                    sl=line.split()
                    elem=sl[0].translate(None, string.digits)
                    opls.append([elem,
                                 pfloat(sl[1]),
                                 pfloat(sl[2]),
                                 sl[0]
                                ])
        for i in nonbond:
            mind=10.**10
            if not elems[i[0]]==elems[i[1]]:
                for j in opls:
                    if elems[i[0]]==j[0]:
                        for k in opls:
                            if elems[i[1]]==k[0]:
                                sigma=math.sqrt(j[1]*k[1])
                                epsilon=math.sqrt(j[2]*k[2])
                                if abs(i[2]-4*epsilon*(sigma**12))*abs(i[3]-4*epsilon*(sigma**6))<mind:
                                    mind=abs(i[2]-4*epsilon*(sigma**12))*abs(i[3]-4*epsilon*(sigma**6))
                                    jmin=j
                                    kmin=k
                atomtypes[i[0]-1].append(jmin)
                atomtypes[i[1]-1].append(kmin)
        for i in atomtypes:
            if len(set([x[3] for x in i]))==1:
                vdw.append(i[0])
            else:
                print "atomtypes not univocal, exiting"
                sys.exit()
                
        with open("lig.lib",'w') as outfile:
            # Fixed this because line splitting was acting up, maybe redo to Marton's code..
            a = sum([entry[2] for entry in charges])
            b = format(a, '.4f')
            outfile.write("*----------------------------------------------------------------------\n"  \
                          "{LIG}                 ! atoms %s charge %s \n\n" % (len(charges), b)        \
                         )
            
            outfile.write("[info]\nSYBYLtype RESIDUE\n\n[atoms]\n")


            for i in charges:
                outfile.write("%i \t %s%i \t T%s%i \t %2.4f\n" % (i[0],     \
                                                                  i[1],     \
                                                                  i[0],     \
                                                                  i[1],     \
                                                                  i[0],     \
                                                                  i[2])
                             ) 

            outfile.write("\n[bonds]\n")

            for i in bonds:
                outfile.write("%s%i %s%i\n" % (elems[i[0]],        \
                                               i[0],elems[i[1]],   \
                                               i[1])
                             )               

            outfile.write("\n[impropers]\n")

            for i in improps:
                outfile.write("%s%i\t%s%i\t%s%i\t%s%i\n" % (elems[i[0]], \
                                                            i[0],        \
                                                            elems[i[2]], \
                                                            i[2],        \
                                                            elems[i[1]], \
                                                            i[1],        \
                                                            elems[i[3]], \
                                                            i[3])
                              )

            outfile.write("\n[charge_groups]\n")

            for i in charges:
                if i[1]!='H':
                    outfile.write("%s%i" % (i[1],i[0]))
                    for j in bonds:
                        if j[0]==i[0] and elems[j[1]]=='H':
                            outfile.write(" H%i" % j[1])
                        if j[1]==i[0] and elems[j[0]]=='H':
                            outfile.write(" H%i" % j[0])
                    outfile.write("\n")
            outfile.write("*----------------------------------------------------------------------\n")

        with open("lig.pdb",'w') as outfile:
            for i in charges:
                if len(i[1]+str(i[0]))<4:
                    outfile.write("ATOM  %5i  %-4sLIG     1    "            \
                                  "%8.3f%8.3f%8.3f\n" % (i[0],              \
                                                         i[1]+str(i[0]),    \
                                                         i[3],              \
                                                         i[4],              \
                                                         i[5])
                                 )
                else:
                    outfile.write("ATOM  %5i %-4s LIG     1    ",           \
                                  "%8.3f%8.3f%8.3f\n" % (i[0],              \
                                                          i[1]+str(i[0]),   \
                                                          i[3],             \
                                                          i[4],             \
                                                          i[5])
                                 )

        with open("lig_w.fep",'w') as outfile:
            outfile.write("[FEP]\nstates 1\n[atoms]\n")
            for i in charges:
                outfile.write("%i %i\n" % (i[0],i[0]))
        with open("Qoplsaa.prm",'w') as outfile:
            # Write vdw parameters
            with open(dir_list[3] + "/Qoplsaa_vdw.prm",'r') as infile:
                for line in infile:
                    outfile.write(line)
            outfile.write("! Ligand vdW parameters\n")
            for i,v in enumerate(vdw,start=1):
                Avdw1=math.sqrt(4.*v[2]*(v[1]**12))
                Bvdw1=math.sqrt(4.*v[2]*(v[1]**6))
                Avdw2=Avdw1/math.sqrt(2.)
                Bvdw2=Bvdw1/math.sqrt(2.)
                outfile.write("T%s%i\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n" %   (v[0], 
                                                                                       i,           
                                                                                       Avdw1,
                                                                                       Avdw1,
                                                                                       Bvdw1,
                                                                                       Avdw2,
                                                                                       Bvdw2,
                                                                                       mass[v[0]]))
            # Write bonds    
            with open(dir_list[3] + "/Qoplsaa_bond.prm",'r') as infile:
                for line in infile:
                    outfile.write(line)

            outfile.write("! Ligand bond parameters\n")
            for i in bonds:
                outfile.write("T%s%i\tT%s%i\t%5.1f\t%5.3f\n" % (elems[i[0]],
                                                                i[0],
                                                                elems[i[1]],
                                                                i[1],
                                                                i[2],
                                                                i[3])
                             )

            # Write angles
            with open(dir_list[3] + "/Qoplsaa_angle.prm",'r') as infile:
                for line in infile:
                    outfile.write(line)

            outfile.write("! Ligand angle parameters\n")
            for i in angles:
                outfile.write("T%s%i\tT%s%i\tT%s%i\t%5.1f\t%5.3f\n" % (elems[i[0]], \
                                                                       i[0],        \
                                                                       elems[i[1]], \
                                                                       i[1],        \
                                                                       elems[i[2]], \
                                                                       i[2],        \
                                                                       i[3],        \
                                                                       i[4])
                             )
            # Write torsions
            with open(dir_list[3] + "/Qoplsaa_torsion.prm",'r') as infile:
                for line in infile:
                    outfile.write(line)

            outfile.write("! Ligand torsion parameters\n")
            for i in torsions:
                if abs(i[4])>0:
                    if abs(i[5])>0 or abs(i[6])>0 or abs(i[7])>0:
                        outfile.write("T%s%i\tT%s%i\tT%s%i\tT%s%i\t%5.3f\t-1.000\t0.000\t1.000\n" %                                                                                 (elems[i[0]],      \
                                                                    i[0],              \
                                                                    elems[i[1]],       \
                                                                    i[1],              \
                                                                    elems[i[2]],       \
                                                                    i[2],              \
                                                                    elems[i[3]],       \
                                                                    i[3],              \
                                                                    i[4])
                               )
                    else:
                        outfile.write("T%s%i\tT%s%i\tT%s%i\tT%s%i\t%5.3f",          \
                                      "\t 1.000\t0.000\t1.000\n" % (elems[i[0]],    \
                                                                   i[0],            \
                                                                   elems[i[1]],     \
                                                                   i[1],            \
                                                                   elems[i[2]],     \
                                                                   i[2],            \
                                                                   elems[i[3]],     \
                                                                   i[3],            \
                                                                   i[4])
                               )

                if abs(i[5])>0:
                            if abs(i[6])>0 or abs(i[7])>0:
                                    outfile.write("T%s%i\tT%s%i\tT%s%i\tT%s%i\t%5.3f\t-2.000\t180.000\t1.000\n"                                                                                % (elems[i[0]],  \
                                                                                  i[0],         \
                                                                                  elems[i[1]],  \
                                                                                  i[1],         \
                                                                                  elems[i[2]],  \
                                                                                  i[2],         \
                                                                                  elems[i[3]],  \
                                                                                  i[3],         \
                                                                                  i[5])
                                                 )
                            else:
                                outfile.write("T%s%i\tT%s%i\tT%s%i\tT%s%i\t%5.3f"               \
                                              "\t 2.000\t180.000\t1.000\n" %  (elems[i[0]],     \
                                                                            i[0],               \
                                                                            elems[i[1]],        \
                                                                            i[1],               \
                                                                            elems[i[2]],        \
                                                                            i[2],               \
                                                                            elems[i[3]],        \
                                                                            i[3],               \
                                                                            i[5])
                                             )
                if abs(i[6])>0:
                    if abs(i[7])>0:
                        outfile.write("T%s%i\tT%s%i\tT%s%i\tT%s%i\t%5.3f",          \
                                      "\t-3.000\t0.000\t1.000\n" % (elems[i[0]],    \
                                                                    i[0],           \
                                                                    elems[i[1]],    \
                                                                    i[1],           \
                                                                    elems[i[2]],    \
                                                                    i[2],           \
                                                                    elems[i[3]],    \
                                                                    i[3],           \
                                                                    i[6])
                                     )
                    else:
                        outfile.write("T%s%i\tT%s%i\tT%s%i\tT%s%i\t%5.3f"           \
                                      "\t 3.000\t0.000\t1.000\n" % (elems[i[0]],    \
                                                                    i[0],           \
                                                                    elems[i[1]],    \
                                                                    i[1],           \
                                                                    elems[i[2]],    \
                                                                    i[2],           \
                                                                    elems[i[3]],    \
                                                                    i[3],           \
                                                                    i[6])
                               )
                if abs(i[7])>0:
                    f.write("T%s%i\tT%s%i\tT%s%i\tT%s%i\t%5.3f" ,           \
                            "\t4.000\t180.000\t1.000\n" % (elems[i[0]],     \
                                                           i[0],            \
                                                           elems[i[1]],     \
                                                           i[1],            \
                                                           elems[i[2]],     \
                                                           i[2],            \
                                                           elems[i[3]],     \
                                                           i[3],            \
                                                           i[7])
                           )
                if (abs(i[4])<0.0001 and 
                    abs(i[5])<0.0001 and 
                    abs(i[6])<0.0001 and 
                    abs(i[7])<0.0001):
                    outfile.write("T%s%i\tT%s%i\tT%s%i\tT%s%i\t%5.3f\t 1.000\t0.000\t1.000\n" % (elems[i[0]],   
                                                                                                i[0],
                                                                                                elems[i[1]],
                                                                                                i[1],
                                                                                                elems[i[2]],
                                                                                                i[2],
                                                                                                elems[i[3]],
                                                                                                i[3],
                                                                                                i[4])
                           )

            # Write impropers
            with open(dir_list[3] + "/Qoplsaa_improper.prm",'r') as infile:
                for line in infile:
                    outfile.write(line)    
            outfile.write("! Ligand improper parameters\n")
            for i in improps:
                outfile.write("T%s%i\tT%s%i\tT%s%i\tT%s%i\t%5.3f\t180.000\n" % (elems[i[0]],
                                                                                i[0],
                                                                                elems[i[2]],
                                                                                i[2],
                                                                                elems[i[1]],
                                                                                i[1],
                                                                                elems[i[3]],
                                                                                i[3],
                                                                                i[4])
                             )


                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='lig_prep',
        version='1.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Generate ligand FEP files for Q. == ')

    
    parser.add_argument('-l', '--ligand',
                        dest = "lig",
                        required = True,
                        help = "name of the docked ligand")
    
    args = parser.parse_args()
    run = Run(lig = args.lig)
    run.make_dir()
    environment = run.get_environment()
    run.get_OPLS(environment)
    run.convert_OPLS(environment)
