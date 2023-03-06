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
a = 4.04e-10        #m, FCC lattice parameter
G = 26.9e9          #Pa, shear modulus of Al
M = 2.5             #-, Taylor number (given in lab PM)
x = 0.9             #-, volume fraction Al (estimted from alloy composition)
b = a/2             #m, Burger vector for FCC ({111} planes)
d = 8.7e-9          #m, diameter of percipitated particles (from SAXS)
gamma = 20e-3       #J m^-2, particle surface energy
alpha = 0.15        #-, given by Canvas
r = d/2             #m, radius of percipitated particles

f = 1-x                         #-, volume fraction of percipitated particles
L = d * (np.pi/(6*f))**(1/3)    #m, distance between percipitated particles
N = L**-3                       #m^-3, number density of percipitated particles

# Calculation and write to file
Delta_sigma_Orowan = Delta_sigma_bend_dislocation_Orowan(M, G, b, L, r)*1e-6 #MPa
Delta_sigma_shear = Delta_sigma_shear_dislocation(M, r, gamma, f, G, b)*1e-6 #Mpa
Delta_sigma_new = Delta_sigma_Canvas(alpha, M, G, b, N, d)*1e-6              #MPa

with open('Al7075 calculation/TIF326_Lab_increased_sigma_Al7075.txt', 'w') as file:
    file.write(f"Delta_sigma_bend_dislocation_Orowan: {Delta_sigma_Orowan:.2e} MPa\n")
    file.write(f"Delta_sigma_shear_dislocation:       {Delta_sigma_shear:.2e} MPa\n")
    file.write(f"Delta_sigma_Canvas:                  {Delta_sigma_new:.2e} MPa")