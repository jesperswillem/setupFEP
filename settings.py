import os
# Root directory of the setup FEP modules
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

# The directories to the input FF and run related input files are given here
FF_DIR = os.path.join(ROOT_DIR, "FF")
INPUT_DIR = os.path.join(ROOT_DIR, "INPUTS")
# Dicionary of locations of Q executables
Q_DIR = {'CSB':'/home/apps/q-5.06/',
         'LOCAL':'/Users/willemjespers/Software/q_510/bin/'
        }
BIN = os.path.join(ROOT_DIR, "bin")
#SCHROD_DIR = '/opt/schrodinger/suites2017-3/'
SCHROD_DIR = '/home/apps/schrodinger2017/'


# CLUSTER INPUTS. To add your own cluster, use the same input as below
CSB = {'NODES'        : '1',
       'NTASKS'       : '16',
       'TIME'         : '0-06:00:00',  # d-hh:mm:ss
       'MODULES'      : 'module load openmpi-x86_64\n module load gcc/6.2.0', 
       'QDYN'         : 'qdyn=/home/apps/q-5.06/qdynp',
       'QPREP'        : '/home/apps/q-5.06/qprep',
       'QFEP'         : '/home/jespers/software/q_510/bin/qfep'
      }

HEBBE = {'NODES'      : '1',
         'NTASKS'     : '20',
         'TIME'       : '0-02:00:00',  # d-hh:mm:ss
         'MODULES'    : 'module load GCC/5.4.0-2.26\nmodule load OpenMPI/1.10.3\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=/c3se/users/jwillem/Hebbe/software/qsource/bin/qdyn5p', #fix qdyn= !!!!!
         'QPREP'      : '/home/jespers/software/q_510/bin/qprep', # NOTE: change to where you are setting up, not where you are running!
         'ACCOUNT'    : 'SNIC2017-1-549'
        }

KEBNE = {'NODES'      : '1',
         'NTASKS'     : '28',
         'TIME'       : '0-04:00:00',  # d-hh:mm:ss
         'MODULES'    : 'module load gompi/2017b\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=/home/w/wije/pfs/software/Q5/bin/qdyn5p', #fix qdyn= !!!!!
         'QPREP'      : '/home/apps/q-5.06/qprep', # NOTE: change to where you are setting up, not where you are running!
         'QFEP'         : '/home/w/wije/pfs/software/q_BAR/bin/qfep5',
         'ACCOUNT'    : 'SNIC2017-12-11'
        }

STALLO = {'NODES'      : '1',
         'NTASKS'     : '20',
         'TIME'       : '0-04:00:00',  # d-hh:mm:ss
         'MODULES'    : 'module load impi/2018.1.163-iccifort-2018.1.163-GCC-6.4.0-2.28\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=/home/jespersw/software/Q6/bin/qdynp', #fix qdyn= !!!!!
         'QPREP'      : '/home/apps/q-5.06/qprep', # NOTE: change to where you are setting up, not where you are running!
         'ACCOUNT'    : 'nn2948k'
        }
