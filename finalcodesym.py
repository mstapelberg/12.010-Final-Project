import numpy as np
import math


#This is a Monte Carlo Simulator for Fusion Reactor
#
#
#
#The Dimensions of the reactor are the following:
#       inner radius of 1.1m
#       length of the cylinder is 20m
#       there are ten cylinders surrounding the vacuum cylinder, each with thickness of .01m

########GEOMETRY CONFIGURATIONS##############
# (y- y_0)^2 + (z - z_0)^2 = R^2
# (y + dv - y_0)^2 + (z + dw - z_0)^2 = R^2

# ybar = y - y_0
# zbar = z - z_0

#(ybar + dv)^2 + (zbar + dw)^2 = R^2

#(v^2 + w^2) d^2 + 2(ybar v + zbar w)d + (ybar^2 + zbar^2 - R^2) = 0

#a = v^2 + w^2, k = ybar*v + zbar w, and c = ybar^2 + zbar^2 -R^2
#ad^2 + 2kd^2 + c = 0 is the distance to the solution

#d = (-k +/- sqrt(k^2 -a*c))/a

#Set user input here
#All inputs are in meters, except for seed and bins which
#are dimensionless
bins = 10
vacradius = 1.1
vesradius = 1.2
length = 20
seed = 1
absorb_radius = np.ndarray(bins, 2, dtype = float)
nalpha_radius = np.ndarray(bins, 2, dtype = float)


def geometry(bins, vacradius, vesradius, length):
    #bins is the number of shells for tallies
    #vacradius is the radius of the vacuum
    #vesradius is the radius of the vacuum vessel
    #length is the length of the cylinder long the y-axis

    increment = (vesradius - vacradius)/2

    #Defines the left and right bounds of the cylinder
    global ybound = length/2
    global ybound = -length/2

    #this defines the boundaries for the vacuum and the shells
    global xbound = np.array(bins+1)
    global zbound = np.array(bins+1)

    #now we specifically define the conditions of the vacuum
    xbound(0) = vacradius
    zbound(0) = vacradius

    #Now it is time to generate the bounds for the shells that surround the vacuum
    #We will pass the vacradius and the x or z value for the bounds test to see where the neutron
    #is.
    for i in range(1:bins):
        xbound(i) = lambda x, vacradius: x >= vacradius + (i-1)*increment \
                and x < vacradius + i*increment
        zbound(i) = lambda z, vacradius: z >= vacradius + (i-1)*increment \
                and z < vacradius + i*increment

#This fills in the values of the nalpha and absorb radius bins
    for i in range(bins):
        absorb_radius[i,0] = vacradius + (i+1)*increment
    for i in range(bins):
        nalpha_radius[i,0] = vacradius + (i+1)*increment

    return xbound, zbound, absorb_radius, nalpha_radius



def bin_sort(interactions, xneutron, yneutron, zneutron):
        """This function takes the 'killed' neutron and prepares it for sorting."""
        #    global bins = 10
        #    global absorb_radius = np.array(bins -1, 2)
        #    global nalpha_radius = np.array(bins -1, 2)

            radius = sqrt(xneutron**2 + yneutron**2 + zneutron**2)
            print(radius)
            if interactions == 0:
                    #THIS PART BELOW IS THROWING ERROR, fixed making radius within bounds
                    #of the problem
                indx = np.searchsorted(absorb_radius[:,0], radius, 'left')
                absorb_radius[indx, 1] = absorb_radius[indx, 1] + 1
                return absorb_radius
            else:
                indx = np.searchsorted(nalpha_radius[:,0], radius, 'left')
                nalpha_radius[indx,1] = nalpha_radius[indx, 1] + 1
                return nalpha_radius

def  calc_theta(seed):
    """Calculates a random theta value for interactions"""
    theta = 2*np.pi*random.random(seed)
    return theta
def calc_phi(seed):
    """calculates a random phi value for interactions"""
    phi = acos(2*random.random(seed)-1)
    return phi
def collision_distance(phi, theta, xneut, zneut, radius):
    """calculates the distance to the next collision, this equation
    is from OPENMC's documentation to calculate the collision distance
    with an infinite sphere parallel to the y axis.

    xbar is difference between the x value at the boundary of the cyliner
    and the x value of the neutron, the same logic is applied to z bar.

    phi and theta values are taken from calc_phi and calc_theta"""
    xbar = lambda x: x - xneut
    zbar = lambda z: z - zneut
    a = phi**2 + theta**2
    k = xbar*phi + zbar*theta
    c = xbar**2 + zbar**2 - radius**2
    d1 = (-k + sqrt(k**2 - a*c))/a
    d2 = (-k - sqrt(k**2 - a*c))/a
    return d1, d2

def simulator(nparticles, ninteractions, vacradius, vesradius):
    """Simulator that cranks out the Monte Carlo Code in Python"""
    for i in range(nparticles):
        #neutron = neutron_func(i)
        energy = 14E6
        phi = calc_phi()
        theta = calc_theta()
        d = collision_distance(phi, theta, xneut, zneut)
        j = 0
        while (j <= ninteractions and neutron_alive = 1)
            interaction = random.random()
            if interaction <= sigma_ngamma(energy)/sigma_t(energy):
                #here we should check which bin the neutron is in
                #and then add to that bin counter
                bin_sorter(0, )
                break
            elif interaction <= sigma_nalpha(energy)/sigma_t(energy):
                #here we should check which bin the neutron is in
                #and then add to that bin counter
                break
            elif interaction <= sigma_elas(energy)/sigma_t(energy):
                #Here we call the scatter function to calculate
                #the new scatter angle and the distance to the
                #next collision
                j++
                continue


#The following functions are used to find cross sections, via binary searches to minimize
#the overhead in these searchs O(ln(n)) compared to O(n) with a linear search
def ni58_total(E):
    """Cross sections for all interactions that can happen in Nickel 58"""
    x = np.loadtxt("NI58SIGMAT.txt", delimiter = ',')
    indx = np.searchsorted(x[:,0], E, side = 'left')
    return x(indx, 1)

def ni58_ngamma(E):
    """Cross sections for the (n, gamma) reactions that can occur in the vacuum vessel"""
    x = np.loadtxt("NI58SIGMAGAMMA.txt", delimiter = ',')
    indx = np.searchsorted(x[:,0], E, side ='left')
    return x(indx, 1)


def ni58_nalpha(E):
    """Cross sections for the (n, alpha) reactions that can occur in the vacuum vessel"""
    x = np.loadtxt("NI58SIGMANALPHA.txt", delimiter = ',')
    indx = np.searchsorted(x[:,0], E, side = 'left')
    return x(indx, 1)

def ni58_elas(E):
    """Cross sections for the ellastic scattering reactions that can occur in the vacuum vessel"""
    x = np.loadtxt("NI58SIGMAS.txt", delimiter = ',')
    indx = np.searchsorted(x[:,0], E, side = 'left')
    return x(indx, 1)
################################END OF CROSS SECTION FUNCTIONS####################################


