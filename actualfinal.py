import numpy as np
import os
import math
import random
import cProfile, pstats, io


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

#d = (-k +/- math.sqrt(k^2 -a*c))/a

#Set user input here
#Put this stuff in a giant function Or do a class if you dare
#All inputs are in meters, except for s'eed and bins which
#are dimensionless
#bins = 10
#vacradius = 1.1
#vesradius = 1.2
#length = 20
#energy = 14E6
#random.seed(10)
ni58total  = np.loadtxt("NI58SIGMAT.txt", delimiter = ',')
ni58ngamma = np.loadtxt("NI58SIGMAGAMMA.txt", delimiter = ',')
ni58nalpha = np.loadtxt("NI58SIGMANALPHA.txt", delimiter = ',')
ni58elas = np.loadtxt("NI58SIGMAS.txt", delimiter = ',')

def profile(fnc):
    """A decorator that uses cProfile to profile a function"""

    def inner(*args, **kwargs):

        pr = cProfile.Profile()
        pr.enable()
        retval = fnc(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream = s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        return retval
    return inner


def geometry(bins, vacradius, vesradius, length, energy):
    #bins is the number of shells for tallies
    #vacradius is the radius of the vacuum
    #vesradius is the radius of the vacuum vessel
    #length is the length of the cylinder long the y-axis


    global absorb_radius
    absorb_radius = np.zeros((bins, 2), dtype = float)
    global nalpha_radius
    nalpha_radius = np.zeros((bins, 2), dtype = float)
    global energy_out
    energy_out = np.zeros((bins, 2), dtype = float)
    increment = (vesradius - vacradius)/bins

    #Defines the left and right bounds of the cylinder
    global ybound
    ybound = length/2
    global xbound
    xbound = vesradius
    global zbound
    zbound = vesradius
    global rbound
    rbound = vesradius
    #this defines the boundaries for the vacuum and the shells
#    global xbound = np.array(bins+1)
#    global zbound = np.array(bins+1)

    #now we specifically define the conditions of the vacuum
#    xbound(0) = vacradius
#    zbound(0) = vacradius

    #Now it is time to generate the bounds for the shells that surround the vacuum
    #We will pass the vacradius and the x or z value for the bounds test to see where the neutron
    #is.
#    for i in range(1:bins):
#        xbound(i) = lambda x, vacradius: x >= vacradius + (i-1)*increment \
#                and x < vacradius + i*increment
#        zbound(i) = lambda z, vacradius: z >= vacradius + (i-1)*increment \
#                and z < vacradius + i*increment

#This fills in the values of the nalpha and absorb radius bins
    for i in range(bins):
        absorb_radius[i,0] = vacradius + (i+1)*increment
    for i in range(bins):
        nalpha_radius[i,0] = vacradius + (i+1)*increment
    for i in range(bins):
        energy_out[i,0] = energy/bins*i

    #xbound and zbound used to be returned here, is this okay to take out??
    return absorb_radius, nalpha_radius, energy_out, energy


def bound_checker(xneut, yneut, zneut, rneut):
    """FUnction that determines whether or not the
    neutron is within the bounds of the geometry"""
    #all takes in a list of booleans, put in list
    if all([abs(yneut) < ybound, abs(xneut) < xbound, abs(zneut) < zbound, abs(rneut) < rbound]):
        return 1
    else:
        return 0

def find_nearest(array, value):
    indx = np.searchsorted(array, value, side = 'left')
    if indx > 0  and (indx == len(array) or math.fabs(value - array[indx-1]) < math.fabs(value - array[indx])):
        return indx-1 #array[indx-1]
    else:
        return indx #array[indx]

def bin_sort(interactions, xneutron, yneutron, zneutron, energy = None):
        """This function takes the 'killed' neutron and prepares it for sorting."""
        #    global bins = 10
        #    global absorb_radius = np.array(bins -1, 2)
        #    global nalpha_radius = np.array(bins -1, 2)
        if interactions == 2 and energy is not None:
            #we do not calculate the radius here because we only need energy
            #indx = np.searchsorted(energy_out[:,0], energy, 'left')
            indx = find_nearest(energy_out[:,0], energy)
            energy_out[indx,1] += 1
            return energy_out
        else:
            radius = math.sqrt(xneutron**2 + yneutron**2 + zneutron**2)
            #print(radius)
            if interactions == 0:
                    #THIS PART BELOW IS THROWING ERROR, fixed making radius within bounds
                    #of the problem
                #indx = np.searchsorted(absorb_radius[:,0], radius, 'right')
                indx = find_nearest(absorb_radius[:,0], radius)
                #print("\n"*4, ">>>>"*9)
                #print(absorb_radius.shape)
                #print("\n"*4)
                absorb_radius[indx, 1] += 1
                return absorb_radius
            elif interactions == 1:
                #indx = np.searchsorted(nalpha_radius[:,0], radius, 'right')
                indx = find_nearest(nalpha_radius[:,0], radius)
                nalpha_radius[indx,1] += 1
                return nalpha_radius

def calc_theta():
    """Calculates a random theta value for interactions"""
    theta = 2*np.pi*random.random()
    return theta
def calc_phi():
    """calculates a random phi value for interactions"""
    phi = math.cos(2*random.random()-1)
    return phi
def surface_distance(phi, theta, xneut, zneut):
    """calculates the distance to the next collision, this equation
    is from OPENMC's documentation to calculate the collision distance
    with an infinite sphere parallel to the y axis.

    xbar is difference between the x value at the boundary of the cyliner
    and the x value of the neutron, the same logic is applied to z bar.

    phi and theta values are taken from calc_phi and calc_theta"""
    #for now because I am lazy
    vacradius = 1.1
    a = phi**2 + theta**2
    k = xneut*phi + zneut*theta
    c = xneut**2 + zneut**2 - vacradius**2
    if k**2 -a*c < 0:
        return 1E20
    d1 = (-k + math.sqrt(k**2 - a*c))/a
    d2 = (-k - math.sqrt(k**2 - a*c))/a
    dist_3d1 = d1/math.sin(phi)
    #dist_3d2 = d2/math.sin(phi)
    return dist_3d1#, dist_3d2

def scatter(energy, A):
    alpha = ((A-1)/(A+1))**2
    newenergy = energy - energy*(1-alpha)*random.random()
    newtheta = calc_theta()
    newphi = calc_phi()
    newr = -np.log(random.random())/ni58_total(energy)
    newx = newr * math.sin(newphi)*math.cos(newtheta)
    newy = newr * math.sin(newphi)*math.sin(newtheta)
    newz = newr * math.cos(newphi)
    bound_checker(newx, newy, newz, newr)
    if bound_checker == 0:
    #this means out of bounds, send energy value, and return 0
        bin_sort(2, newx, newy, newz, newenergy)
        return None #test if it received a value of none
        # to decide whether or not hte neutron stays
        #HELP HERE
    else:
        newd = surface_distance(newphi, newtheta, newx, newz)
        return newtheta, newphi, newr, newd, newenergy
@profile
def simulator(nparticles, ninteractions, bins, vacradius, vesradius,length, energy_i):
    """Simulator that cranks out the Monte Carlo Code in Python"""
    geometry(bins, vacradius, vesradius, length, energy_i)
    for i in range(nparticles):
        energy = energy_i
        xneut = 0
        yneut = 0
        zneut = 0
        j = 0
        neutron_alive = 1
        while (j <= ninteractions and neutron_alive == 1):
            theta = calc_theta()
            phi = calc_phi()
            d = surface_distance(phi, theta, xneut, zneut)
            #d2 = surface_distance(phi, theta, xneut, zneut)[1]
            r = -np.log(random.random())/ni58_total(energy)
            if r > d and r-d != 0:
                #This if statement checks if the vacuum surface is closer
                #or if the distance to the next collision is
                #The following statements move the particle to the surface
                #of the vacuum

                xneut = d * math.sin(phi)*math.cos(theta)
                yneut = d * math.sin(phi)*math.sin(theta)
                zneut = d * math.cos(phi)
                #CHECK HERE
                continue

            elif r < d:
                xneut = r * math.sin(phi)*math.cos(theta)
                yneut = r * math.sin(phi)*math.sin(theta)
                zneut = r * math.cos(phi)
                interaction = random.random()
                if interaction <= ni58_ngamma(energy)/ni58_total(energy):
                    #here we should check which bin the neutron is in
                    #and then add to that bin counter
                    bin_sort(0, xneut, yneut, zneut )
                    neutron_alive = 0
                elif interaction <= ni58_nalpha(energy)/ni58_total(energy):
                    #here we should check which bin the neutron is in
                    #and then add to that bin counter
                    bin_sort(1, xneut, yneut, zneut)
                    neutron_alive = 0
                elif interaction <= ni58_elas(energy)/ni58_total(energy):
                    #Here we call the scatter function to calculate
                    #scatter(energy, 58) #future implementations will have more elements
                    scatter_output = scatter(energy,58)
                    theta = scatter_output[0]
                    phi = scatter_output[1]
                    r = scatter_output[2]
                    d = scatter_output[3]
                    energy = scatter_output[4]
                    #for now we have nickel 58
                    #the new scatter angle and the distance to the
                    #thenext collision

                    j+= 1
                    continue
    return nalpha_radius, absorb_radius, energy_out

#def plotter(nalpha_radius, absorb_radius, energy_out):
    #import matplotlib here
    """Plot the three arrays, with the fist column as the x values,
    and the second column as the tallies per bin. These are histograms"""
#    return 0

#The following functions are used to find cross sections, via binary searches to minimize
#the overhead in these searchs O(ln(n)) compared to O(n) with a linear search
def ni58_total(E):
    """Cross sections for all interactions that can happen in Nickel 58"""
    indx = np.searchsorted(ni58total[:,0], E, side = 'left')
    return ni58total[indx, 1]

def ni58_ngamma(E):
    """Cross sections for the (n, gamma) reactions that can occur in the vacuum vessel"""
    indx = np.searchsorted(ni58ngamma[:,0], E, side ='left')
    return ni58ngamma[indx, 1]


def ni58_nalpha(E):
    """Cross sections for the (n, alpha) reactions that can occur in the vacuum vessel"""
    indx = np.searchsorted(ni58nalpha[:,0], E, side = 'left')
    return ni58nalpha[indx, 1]

def ni58_elas(E):
    """Cross sections for the ellastic scattering reactions that can occur in the vacuum vessel"""
    indx = np.searchsorted(ni58elas[:,0], E, side = 'left')
    return ni58elas[indx, 1]
################################END OF CROSS SECTION FUNCTIONS####################################


