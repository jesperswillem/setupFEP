#import pymoldyn
import os

#SCRIPT_DIR = os.path.dirname(os.path.realpath(pymoldyn.__file__))

# This is the dir where pymoldyn git repo has been deployed,
# or to be more specific, this file settings.py is located
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

#This dir now has to be the absolute path to the source of templates
# You can refer it relatively like this:
REPO_DIR = os.path.join(ROOT_DIR, "templates")
# Or absolutely like this:
# REPO_DIR = "/path/to/your/templates"
# But this WILL FAIL: * REPO_DIR = "templates" *

# The directories to the input FF and run related input files are given here
FF_DIR = os.path.join(ROOT_DIR, "FF")
INPUT_DIR = os.path.join(ROOT_DIR, "INPUTS")

SCHROD_DIR = '/home/apps/schrodinger2017/'

# Choose a path to the gromacs binaries.
#GROMACS_PATH = "/opt/applications/gromacs/4.0.5/gnu/ib/bin/"
#GROMACS_PATH = "/opt/applications/gromacs/4.0.5/gnu/gige/bin/"
#GROMACS_PATH = "/opt/cesga/gromacs-4.0.7/bin/"
#GROMACS_PATH = "/opt/gromacs405/bin/"                               #cuelebre.inv.usc.es
GROMACS_PATH = "/home/apps/gromacs-4.6.7/bin/"                      #csb.bmc.uu.se
#GROMACS_PATH = "/software/apps/gromacs/4.6.3/g472/bin/"             #Triolith
#GROMACS_PATH = "/sw/bin/"                                           #Standalone in Mac Fink
#GROMACS_PATH = "/Users/esguerra/software/gromacs-4.6.5/bin/"        #Standalone in Mac
#GROMACS_PATH = "/c3se/apps/Glenn/gromacs/4.6.3-p20130821-gcc48/bin" #Glenn at Chalmers
#GROMACS_PATH = "/c3se/apps/Glenn/gromacs/5.0.4-gcc48-cuda/bin/"     #Glenn GPU at Chalmers
#GROMACS_PATH = "/sw/apps/gromacs/4.6.3/tintin/bin"                  #Tintin
#GROMACS_PATH = "/lap/gromacs/4.6.5/bin"                             #Abisko

# FEP related inputs
REPLICATES='10'
WINDOWS='20'
TEMPERATURE='298'

# Lambda sampling options
#SAMPLING='linear'
SAMPLING='sigmoidal'
#SAMPLING='exponential'
#SAMPLING='reverse_exponential'

# CLUSTER INPUTS. To add your own cluster, use the same input as below
CSB = {'NODES'      : '1',
       'NTASKS'     : '16',
       'TIME'       : '0-03:00:00',  # d-hh:mm:ss
       'MODULES'    : 'module load openmpi-x86_64\n', # Add a \n for every added module
       'QDYN'       : 'qdyn=/home/apps/q-5.06/qdynp',
       'QPREP'      : '/home/apps/q-5.06/qprep'
      }

HEBBE = {'NODES'      : '1',
         'NTASKS'     : '20',
         'TIME'       : '0-02:00:00',  # d-hh:mm:ss
         'MODULES'    : 'module load GCC/5.4.0-2.26\nmodule load OpenMPI/1.10.3\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=/c3se/users/jwillem/Hebbe/software/qsource/bin/qdyn5p', #fix qdyn= !!!!!
         'QPREP'      : '/home/apps/q-5.06/qprep', # NOTE: change to where you are setting up, not where you are running!
         'ACCOUNT'    : 'SNIC2017-1-549'
        }

KEBNE = {'NODES'      : '1',
         'NTASKS'     : '28',
         'TIME'       : '0-04:00:00',  # d-hh:mm:ss
         'MODULES'    : 'module load gompi/2017b\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=/home/w/wije/pfs/software/Q5/bin/qdyn5p', #fix qdyn= !!!!!
         'QPREP'      : '/home/apps/q-5.06/qprep', # NOTE: change to where you are setting up, not where you are running!
         'ACCOUNT'    : 'SNIC2017-1-549'
        }
