##############################################################################################
## PROJECT 2: SIMULATION OF TURBIE ##########################################################
##############################################################################################

# This script constructs a function that simulates and builds Turbie - a a simple, 
# two-degree-of-freedom (2DOF) system based of the DTU 10 MW Reference Wind Turbine

"Check for NBs before handing in"


#Variables that needs to come from the function
V = 20 # Wind speed: V(m/s)
TI_cat = 0.1 #TI category = 0.1, 0.05 or 0.15

#Import relevating packages
import numpy as np
import pandas as pd
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

"NB: LOAD THE DATA IS AS DATA INSTEAD OF TEXT"

# We import the Turbie parameters used in the function: 
with open('turbie_parameters.txt', 'rt') as t_parameters:
    content = t_parameters.readlines()
    m_blade = float(content[1].split('#')[0].strip())
    m_nacelle = float(content[2].split('#')[0].strip())
    m_tower = float(content[3].split('#')[0].strip())
    m_hub = float(content[4].split('#')[0].strip())
    c1 = float(content[5].split('#')[0].strip())
    c2 = float(content[6].split('#')[0].strip())
    k1 = float(content[7].split('#')[0].strip())
    k2 = float(content[8].split('#')[0].strip())
    D_rotor = float(content[13].split('#')[0].strip())
    rho = float(content[14].split('#')[0].strip())
    t_parameters.closed
True

# We define m1(mass of the three blades) and m2(combined effects of the nacelle, hub and tower)
m1 = 3*m_blade
m2 = m_nacelle + m_hub + m_tower

# We define the Mass Matrix (M), the Damping Matrix (C), and the Stiffness Matrix (K) for the Turbie system:
M = np.array([[m1, 0],[0, m2]])
C = np.array([[c1, -c1],[-c1, c1 + c2]])
K = np.array([[k1, -k1],[-k1, k1 + k2]])

# We import the Turbie parameters used in the function: 
ct_df = pd.read_csv('CT.txt', sep='\t')

#We look-up Ct from the user input wind speed
Ct = ct_df[ct_df['# V(m/s)'] == V]['CT'].values[0]

#We estimate the the Swept Area A
radius = D_rotor/2
A = radius**2*math.pi

# We define number of freedom
N = 2

#We look up the wind profile based on TI and wind, and load it into the wind data frame
wind_df = pd.read_csv(f'wind_files/wind_TI_{TI_cat}/wind_{V}_ms_TI_{TI_cat}.txt', sep='\t')
print(wind_df)

#We define the system martrix (vector A)
ZeroN = np.zeros((N, N))
I = np.eye(N)
Minv = np.linalg.inv(M)

A = np.array([[ZeroN, I], [-Minv@K, -Minv@C]])

print(A)




# We define the aerodynamic forcing on the blades: 
f_aero = 1/2*rho*Ct*A
print(f_aero)

# We define the forcing vector F(t)
Ft = np.array([[f_aero],[0]])

