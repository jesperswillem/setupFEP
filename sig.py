import settings as s

def sigmoid(t, k):
    return ( k * t / (k -t + 1.0))

def get_lambda(windows, sampling):
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
            lmbda1 = 0.5 * (sigmoid(float(i)/float(step), k) + 1)  
            lmbda2 = 0.5 * (-sigmoid(float(i)/float(step), k) + 1)
            lmbda_1.append(lmbda1)
            lmbda_2.append(lmbda2)

        lmbda_2 = lmbda_2[1:]

        for i in reversed(lmbda_2):
            lambdas.append(i)

        for i in lmbda_1:
            lambdas.append(i)

    else:
        for i in range(0, windows + 1):
            lmbda = sigmoid(float(i)/float(windows), k)
            lambdas.append(lmbda)

    return lambdas