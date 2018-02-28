import glob
import numpy as np

def get_results():
    results = {}
    residues = []
        
    for filename in glob.glob('*/*/FEP_*/FEP*/*/*'):
    #for filename in glob.glob('2.CH3/1.protein/FEP_17-73/FEP1/298/*'):
        line = filename.split('/')
        run = line[0]
        FEP = line[1]
        name = line[2] + '_' + line[1].split('.')[1]
        replicate = line[-1]
        ERROR = False

        with open(filename + '/qfep.out') as infile:
            block = 0
            for line in infile:
                line = line.split()
                if len(line) > 3:
                    if line[0] == 'ERROR:':
                        ERROR = True
                    
                    if line[3] == 'Free':
                        block = 1

                    if line[3] == 'Termodynamic':
                        block = 2

                    if line[3] == 'Overlap':
                        block = 3

                    if line[3] == 'BAR':
                        block = 4

                    if line[3] == 'Reaction':
                        block = 0

                if len(line) > 1:
                    if block == 1:
                        if line[0] == '1.000000':
                            dGr = line[4]

                        elif line[0] == '0.000000':
                            dGf = line[2]
                            
                            if line[5] == '-Infinity':
                                dG = np.nan
                                
                            else:
                                dG = float(line[5])

                    if block == 2 and line[0] == '0.000000':
                        dG_TI = line[2]

                    if block == 3 and line[0] == '0.000000':
                        if line[2] == '-Infinity':
                            dG_overlap = np.nan
                            
                        else:
                            dG_overlap = float(line[2])

                    if block == 4 and line[0] == '0.000000':
                        if line[2] == '-Infinity':
                            dG_BAR = np.nan
                        else:
                            dG_BAR = float(line[2])
        
        if ERROR != True:
            data = [name, replicate, dG, dG_overlap, dG_BAR]
            
        else:
            data = [name, replicate, np.nan, np.nan, np.nan]

        if name in results:
            if len(results[name]) < 10:    # temp fix
                results[name].append(data)
        else:
            results[name] = [data]

    return results

def calc_sum_error(data):
    data = np.array(data)
    mean = np.nanmean(data)
    sem = np.nanstd(np.array(data), ddof =1)
    
    return mean, sem

def calc_ddG(raw_data):
    for name in raw_data:
        scoring = [[], [], []]
        for line in raw_data[name]:
            scoring[0].append(line[2])
            scoring[1].append(line[3])
            scoring[2].append(line[4])
            
        for data in scoring:
            dG = calc_sum_error(data)
            
            print name, dG

raw_data = get_results()
calc_ddG(raw_data)
