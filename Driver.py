from actualfinal import simulator, scatter, surface_distance, bin_sort, calc_theta, calc_phi, bound_checker, geometry, ni58_total, ni58_ngamma, ni58_nalpha, ni58_elas, profile
import cProfile, pstats, io
import matplotlib.pyplot as plt


#USER INPUTS
#The inputs are, nparticles, ninteractions, bins, vacradius, vesradius, length, energy_i
nparticles = 1000
ninteractions = 250
bins = 10
vacradius = 1.1
vesradius = 1.2
length = 20
energy_i = 14E6

#END OF USER INPUT
#Runs the simulation code
output = simulator(nparticles, ninteractions, bins, vacradius, vesradius, length, energy_i)

#These prints are for tables or numerical verification.
print("This is the nalpha_radius for neutrons")
nalpha = output[0]
print(nalpha)

print("\n \n This is the absorb_radius for neutrons")
absorb = output[1]
print(absorb)

print("\n \n This is the energy_out for neutrons")
energyo = output[2]
print(energyo)

#Plot for Nalpha
x1 = nalpha[:, 0]
x2 = absorb[:, 0]
x3 = energyo[:,0]

y1 = nalpha[:, 1]
y2 = absorb[:, 1]
y3 = energyo[:,1]

#fig, axs = plt.subplots(1,2, sharey = True, tight_layout = True)
#plt.subplot(2,1,1)
#plt.plot(x1, y1, 'o-')
#plt.title("(n,alpha) Production")
#plt.ylabel("Helium Production (# atoms)")
#plt.xlabel("Radial Distance from Center of Vacuum Vessel")

#plt.subplot(2,1,2)
#plt.plot(x2, y2, '.-')
#plt.title("Neutron Absorption")
#plt.ylabel("Number of Neutron Absorptions")
#plt.xlabel("Radial Distance from Center of Vacuum Vessel")

plt.subplot()
plt.xscale('log')
plt.plot(x3, y3, '--')
plt.title("Exit Energy Distribution of Neutons")
plt.ylabel("Number of Neutrons")
plt.xlabel("Energy (eV)")

#ax.plot(x,y)
#plt.yscale('log')
#ax.set (xlabel='Radial distance from the center of the Cylinder (m)',
#        ylabel = '# of He atoms produced',
#        title ='(n, alpha) production in Vacuum Vessel')
#plt.hist(x,y, bins = 10)
plt.show()

