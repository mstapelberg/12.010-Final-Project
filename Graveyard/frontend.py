import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import MCdriver
import math
from itertools import product, combinations

maxneutrons = 10
maxinteractions = 10
atomic_number = 58
a = 1.1
#values = mcdriver.MC(maxneutrons, maxinteractions, atomic_number, a)
#print(values)

def mean_scatter_length():
#this can be found by doing the following:
    #

#evaluate the average distance fom the center of the sphere to absorption
#check the index at which the neutron is absorbed, then use that index
#to find the position of the

#evaluate the average position of helium production within the vacuum vessel
#do this by using the neutron_nalpha array indices,
#find the indices that are not NAN,
#plug in the number of that index to find the position of the neutron at that time.
