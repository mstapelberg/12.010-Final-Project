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



def geometry(bins, vacradius, vesradius, length):
    #bins is the number of shells for tallies
    #vacradius is the radius of the vacuum
    #vesradius is the radius of the vacuum vessel
    #length is the length of the cylinder long the y-axis

    #Defines the left and right bounds of the cylinder
    yboundu = length/2
    yboundl = -length/2

    #this defines the boundaries for the vacuum and the shells
    xboundu = np.array(bins+1)
    xboundl = np.array(bins+1)
    zboundu = np.array(bins+1)
    zboundl = np.array(bins+1)

    #now we specifically define the conditions of the vacuum
    xboundu(0) = vacradius
    xboundl(0) = -vacradius
    zboundu(0) = vacradius
    zboundl(0) = -vacradius

    #Now it is time to generate the bounds for the shells that surround the vacuum
    for i in range(1:bins):
        xboundu(i) = vesradius + increment
        xboundl(i) = vesradius - increment



