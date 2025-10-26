##############################################################################################
## PROJECT 2: SIMULATION OF TURBIE ##########################################################
##############################################################################################

# This script constructs a function that simulates and builds Turbie - a a simple, 
# two-degree-of-freedom (2DOF) system based of the DTU 10 MW Reference Wind Turbine

"Check for NBs before handing in"


#Variables that needs to come from the function - What we should be looping over
V = 20 # Wind speed: V(m/s)
TI_cat = 0.1 #TI category = 0.1, 0.05 or 0.15
t = 0.010

#Import relevating packages
import numpy as np
import pandas as pd
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

"NB: LOAD THE DATA IS AS DATA INSTEAD OF TEXT"

"I need to loop over the different files using f-string to load and save the results"
"""for ws in [3, 4, 5]:
     file_name = f"results_for_ws_{ws}.txt"
     print(file_name) """

"""We need the forloop to run the simulation looping over different wind speeds"""

"""I need to use the cypy function to solve this problem"""

"""We have to simulate for each of the wind farms / data set"""

### IMPORTING DATASETS

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

# We import the Turbie parameters used in the function: 
ct_df = pd.read_csv('CT.txt', sep='\t')

#We look up the wind profile based on TI and wind, and load it into the wind data frame
wind_df = pd.read_csv(f'wind_files/wind_TI_{TI_cat}/wind_{V}_ms_TI_{TI_cat}.txt', sep='\t')

## DEFINE AND BUILD THE M, C, AND K MATRICES

# We define m1(mass of the three blades) and m2(combined effects of the nacelle, hub and tower)
m1 = 3*m_blade
m2 = m_nacelle + m_hub + m_tower

# We define the Mass Matrix (M), the Damping Matrix (C), and the Stiffness Matrix (K) for the Turbie system:
M = np.array([[m1, 0],[0, m2]])
C = np.array([[c1, -c1],[-c1, c1 + c2]])
K = np.array([[k1, -k1],[-k1, k1 + k2]])

### BUILD A FUNCTION TO DETERMINE Y'(T)

## Define aerodynamic parameter not defined in parameters dataset
radius = D_rotor/2
SW_A = radius**2*math.pi #Swept Area
C_t = ct_df[ct_df['# V(m/s)'] == V]['CT'].values[0] #We look-up Ct from the user input wind speed
u_t = wind_df[wind_df['Time(s)'] == 0.01]['V(m/s)'].values[0]

## Set up function to determine y'(t)
def y_2dof(t, y, M, C, K):
    x1, x2, x1_dot, x2_dot = y
    
    # Aerodynamic force on mass 1
    v_rel = u_t - x1_dot
    f_aero = 1/2 * rho * C_t * SW_A * v_rel * abs(v_rel)
    # No external force on mass 2
    f_2 = 0
        
    # Current positions and velocities
    x = np.array([x1, x2])
    x_dot = np.array([x1_dot, x2_dot])
    F_t = np.array([f_aero, f_2])
    
    # Accelerations: M*ddot{x} + C*dot{x} + K*x = F  =>  ddot{x} = M^-1*(F - C*dot{x} - K*x)
    x_ddot = np.linalg.inv(M) @ (F_t - C @ x_dot - K @ x)

    ## We define parameters to define y'(t)
    N = 2 # Number of freedom

    ZeroN = np.zeros((N, N))
    I = np.eye(N)
    Minv = np.linalg.inv(M)

    #We define the system martrix (vector A) and input matrix (vector B)
    A = np.block([[ZeroN, I], [-Minv@K, -Minv@C]])
    B = np.block([0],-Minv@F_t)

    y_prime = A@y + B
    
    return [x1_dot, x2_dot, x_ddot[0], x_ddot[1]]

# We set the timespan
t_span = (0, 100) #659.990
t_eval = np.linspace(*t_span, 100) #66000

# We set initial conditions for y'(t)
y0 = [0,0,0,0]

res = solve_ivp(y_2dof, t_span, y0, t_eval=t_eval, args=(M,C,K))

#Plot the results: 
import matplotlib.pyplot as plt
plt.plot(res.t, res.y[0], label='x1 (mass 1)')
plt.plot(res.t, res.y[1], label='x2 (mass 2)')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.title('2-DOF Mass-Spring-Damper with Aerodynamic Force')
plt.legend()
plt.grid(True)
plt.show()




