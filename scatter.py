import numpy as np
import random

seed = 1
def scatter(energy, A):
    alpha = ((A-1)/(A+1))**2
    newenergy = energy - energy*(1-alpha)*random.random(seed)
    newtheta = calc_theta()
    newphi = calc_phi()
    newr = -np.log(random.random(seed))/ni58_total(energy)
    x = newr * sin(newphi)*cos(newtheta)
    y = newr * sin(newphi)*sin(newtheta)
    z = newr * cos(newtheta)
    bound_checker(x, y, z, newr)
    if bound_checker == 0
        #this means out of bounds, send energy value, and return 0
        bin_sort(2, x, y, z, newenergy)
    #HELP HERE
    newd = collision_distancei(newphi, newtheta, newx, newz)

    return newtheta, newphi, newr, newd, newenergy
