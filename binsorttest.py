import numpy as np
from math import sqrt


vacradius = 1.1
vesradius = 1.2
bins = 10
increment = (vesradius - vacradius)/bins
absorb_radius = np.ndarray(shape = (bins, 2), dtype = float)
nalpha_radius = np.ndarray(shape = (bins, 2), dtype = float)
for i in range(bins):
    absorb_radius[i,0] = vacradius + (i+1)*increment
for i in range(bins):
    nalpha_radius[i,0] = vacradius + (i+1)*increment
#kprint(absorb_radius)

def bin_sort(interactions, xneutron, yneutron, zneutron):
    """This function takes the 'killed' neutron and prepares it for sorting."""
#    global bins = 10
#    global absorb_radius = np.array(bins -1, 2)
#    global nalpha_radius = np.array(bins -1, 2)

    radius = sqrt(xneutron**2 + yneutron**2 + zneutron**2)
    print(radius)
    if interactions == 0:
        #THIS PART BELOW IS THROWING ERROR, fixed making radius within bounds
        #of the problem.
        indx = np.searchsorted(absorb_radius[:,0], radius, 'left')
        absorb_radius[indx, 1] = absorb_radius[indx, 1] + 1
        return absorb_radius

    else:
        indx = np.searchsorted(nalpha_radius[:,0], radius, 'left')
        nalpha_radius[indx,1] = nalpha_radius[indx, 1] + 1
        return nalpha_radius

y = bin_sort(0, .6 , .65, .7)
print(absorb_radius)




