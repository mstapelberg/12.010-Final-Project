import numpy as np
import random

seed = 1
def scatter(energy, A):
    alpha = ((A-1)/(A+1))**2
    newenergy = energy - energy*(1-alpha)*random.random(seed)
    newtheta = calc_theta()
    newphi = calc_phi()
    newr = -np.log(random.random(seed))/ni58_total(energy)
    #HELP HERE
    newd = collision_distancei(newphi, newtheta, newx, newz, radius)

    return newtheta, newphi, newr, newd, newenergy
