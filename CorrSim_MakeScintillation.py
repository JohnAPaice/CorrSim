## A selection of functions for use with CorrSim.

import numpy as np
from sys import exit

# calculate_scintillation_noise(empirical_value, telescope_diameter, exposure_time, target_altitude, telescope_altitude, turbulence_height):
#     zenith_distance = ( (90-target_altitude) * 2*np.pi) / 360    #Convert to Zenith Distance in Radians

##    ===============================================================
##    ===============         Input Parameters        ===============

print("")
print("-----------------")
print("Welcome to CorrSim_MakeScintillation.py!")
print("Calculates a parameter for atmospheric scintillation.")
print("(Uses Eqn. 7 of Osborn+2015)")
print("")

C_Y = input("Define Empirical Value (C_Y > 0; default = 1.5): ")
if C_Y == "":
    C_Y = '1.5'

D = input("Define Telescope Diameter in metres (D > 0; default = 3.6): ")
if D == "":
    D = '3.6'

t = input("Define Exposure Time in seconds (t > 0; default = 0.09): ")
if t == "":
    t = '0.09'

gamma = input("Define Zenith Distance in degrees (90 > gamma > 0; default = 20): ")
if gamma == "":
    gamma = '20'

h_obs = input("Define Telescope Altitude in metres (h_obs > 0; default = 2500): ")
if h_obs == "":
    h_obs = '2500'

H = input("Define Turbulence Height in metres (H > 0; default = 8000): ")
if H == "":
    H = '8000'


C_Y          = float(C_Y)
D            = float(D)
t            = float(t)
gamma        = float(gamma)
h_obs        = float(h_obs)
H            = float(H)

try:
    scin_noise = 10 * 10**(-6) * C_Y**2 * D**(-4/3) * t**(-1) * \
    np.cos(np.deg2rad(gamma))**(-3) * np.exp((-2*h_obs)/H)
except Exception as e:
    print("")
    print("!WARNING! Calculation Failed:")
    print(e)
    print("Please check your parameters.")
    print("No parameter should be 0 or a negative number.")
    exit()


print("")
print("==========")
print("Scintillation Noise Parameter:")

print(scin_noise)    ## Osborn+2015, Eqn. 7

if scin_noise <= 0:
    print("WARNING: Zero or negative value for scin_noise.")
    print("Did you get all the parameters correct?")
