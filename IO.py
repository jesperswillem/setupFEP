import re
import shlex
from subprocess import check_output
import os

def pdb_parse_in(line, include=('ATOM','HETATM')):
    """
    Takes a pdb file line and parses it into a list, according to Atomic Coordinate Entry Format 
    v3.3
    """
    at_entry = []

    if line.startswith(include):
        at_entry.append(line[0:6])              #  0 ATOM/HETATM
        at_entry.append(int(line[6:11]))        #  1 ATOM serial number
        at_entry.append(line[12:16])            #  2 ATOM name
        at_entry.append(line[16:17])            #  3 Alternate location indicator
        at_entry.append(line[17:20])            #  4 Residue name
        at_entry.append(line[21:22])            #  5 Chain identifier
        at_entry.append(int(line[22:26]))       #  6 Residue sequence number
        at_entry.append(line[26:27])            #  7 Code for insertion of residue
        at_entry.append(float(line[30:38]))     #  8 Orthogonal coordinates for X
        at_entry.append(float(line[38:46]))     #  9 Orthogonal coordinates for Y
        at_entry.append(float(line[46:54]))     # 10 Orthogonal coordinates for Z
        # These entries can be empty
        try:
            at_entry.append(float(line[54:60])) # 11 Occupancy
            
        except:
            at_entry.append(0.0)                # 11 Empty Occupancy
            
        try:
            at_entry.append(float(line[60:66])) # 12 Temperature factor
            
        except:
            at_entry.append(0.0)                # 12 Empty Temperature factor
            
        try:
            at_entry.append(line[76:78])        # 13 Element symbol
            
        except:
            at_entry.append('  ')               # 13 Empty Element symbol
            
        try:
            at_entry.append(line[78:80])        # 14 Charge on atom
            
        except:
            at_entry.append('  ')               # 14 Empty charge
        
    else:
        at_entry = line
    
    return at_entry
    
def pdb_parse_out(line):
    """
    Takes a list and parses it into a pdb writeable line
    """
    line = '{:6s}{:5d} {:^4s}{:1s}{:4s}{:1s}{:4d}{:1s}   '\
           '{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(*line)
    return line

def replace(string, replacements):
    pattern = re.compile(r'\b(' + '|'.join(replacements.keys()) + r')\b')
    replaced_string = pattern.sub(lambda x: replacements[x.group()], string)
    return replaced_string

def run_command(executable, options, string = False):
    """
    Takes three variables, the executable location and its options as strings and a tag if the
    options need to be split or not (e.g. Q runs with one string), and runs the program.
    Returns the output of that program as an unformatted string.
    """
    if string == False:
        args = shlex.split(executable + options)
        out = check_output(args)

    else:
        os.system(executable + options)
        out = None

    return out

def AA(AA):
    """
    Handy dictionary to convert 3 letter AA code to one and vice versa
    """
    threeAA = {'CYS': 'C', 'CYX': 'C', 'ASH': 'D', 'ASP': 'D', 'SER': 'S', 
               'GLN': 'Q', 'LYN': 'K', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 
               'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HID': 'H', 
               'HIP': 'H', 'HIE': 'H', 'HIS': 'H', 'LEU': 'L', 'ARN': 'R', 
               'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLH': 'E', 
               'GLU': 'E', 'TYR': 'Y', 'MET': 'M',}
    
    fourAA = { 'CCYS': 'C', 'CASP': 'D', 'CASH': 'H', 'CSER': 'S', 
               'CGLN': 'Q', 'CLYN': 'K', 'CLYS': 'K', 'CILE': 'I', 
               'CPRO': 'P', 'CTHR': 'T', 'CPHE': 'F', 'CASN': 'N', 
               'CGLY': 'G', 'CHIE': 'H', 'CHID': 'H', 'CHIP': 'H', 
               'CLEU': 'L', 'CARG': 'R', 'CARN': 'R', 'CTRP': 'W', 
               'CALA': 'A', 'CVAL': 'V', 'CGLU': 'E', 'CGLH': 'E',
               'CTYR': 'Y', 'CMET': 'M'}
    
    oneAA = {v: k for k, v in threeAA.iteritems()}
    
    if len(AA) == 4:
        AA = fourAA[AA]
    
    if len(AA) == 3:
        AA = threeAA[AA]
        
    if len(AA) == 1:
        AA = oneAA[AA]
    return(AA)

    