from finalfunctions import simulator, scatter, surface_distance, bin_sort, calc_theta, calc_phi, bound_checker, geometry, ni58_total, ni58_ngamma, ni58_nalpha, ni58_elas, profile
import cProfile, pstats, io

#The inputs are, nparticles, ninteractions, bins, vacradius, vesradius, length, energy_i
print(simulator(2, 250, 10, 1.1, 1.2,20, 14E6)[2])

