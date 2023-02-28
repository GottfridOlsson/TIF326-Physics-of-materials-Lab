##====================================================##
##     Project: TIF326 LAB HEAT TREATMENT OF Al-7075
##        File: formatting-raw-CSV.py
##      Author: GOTTFRID OLSSON 
##     Created: 2023-02-21
##     Updated: 2023-02-22
##       About: Calculating Delta sigma for Al7075
##====================================================##
import numpy as np

# Equations to estimate increase in stress (sigma)
def Delta_sigma_bend_dislocation_Orowan(M, G, b, L, r):
    return M*G*b/(L-2*r)

def Delta_sigma_shear_dislocation(M, r, gamma, f, G, b):
    return (2*M/r)*np.sqrt(gamma**3*f/(G*b))

def Delta_sigma_Canvas(alpha, M, G, b, N, d):
    return alpha*M*G*b*np.sqrt(N*d)


# Material properties Al 7075
a = 4.04*1e-10      #m, FCC lattice parameter
G = 26.9*1e9        #Pa, shear modulus
M = 2.5             #-, Taylor number (given in lab PM)
f = 0.9             #-, volume fraction Al (estimted from alloy composition)
b = a/2             #m, Burger vector for FCC ({111} planes)
L = 12345           #m, distance between percipitated particles # TODO
gamma = 20*1e-3     #J m^-2, particle surface energy
alpha = 0.15        #-, given by Canvas
d = 8.7*1e-9        # m, diameter of percipitated particles
r = d/2         
N = 12345           # m^-3, number density of percipitated particles # TODO 

# Calculation and write to file
Delta_sigma_Orowan = Delta_sigma_bend_dislocation_Orowan(M, G, b, L, r)*1e-6 #MPa
Delta_sigma_shear = Delta_sigma_shear_dislocation(M, r, gamma, f, G, b)*1e-6 #Mpa
Delta_sigma_new = Delta_sigma_Canvas(alpha, M, G, b, N, d)*1e-6              #MPa

with open('Al7075 calculation/TIF326_Lab_increased_sigma_Al7075.txt', 'w') as file:
    file.write(f"Delta_sigma_bend_dislocation_Orowan: {Delta_sigma_Orowan:.2e} MPa\n")
    file.write(f"Delta_sigma_shear_dislocation:       {Delta_sigma_shear:.2e} MPa")
    file.write(f"Delta_sigma_Canvas:                  {Delta_sigma_new:.2e} MPa")