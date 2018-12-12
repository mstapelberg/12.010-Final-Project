import numpy as np

def ni58_total(E):
    x = np.loadtxt("NI58SIGMAT.txt", delimiter = ',')
    #x1 = x[:,0]
    indx = np.searchsorted(x[:,0], E, side = 'left')
    return x[indx, 1]



print(ni58_total(1000))
#f = np.loadtxt("NI58SIGMAT.txt", delimiter = ',')
#print(f.shape)
#print(f)
