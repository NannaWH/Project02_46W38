##############################################################################################
## PROJECT 2: SIMULATION OF TURBIE ##########################################################
##############################################################################################

# This script constructs a function that simulates and builds Turbie - a a simple, 
# two-degree-of-freedom (2DOF) system based of the DTU 10 MW Reference Wind Turbine

#Import relevant packages
import numpy as np
import pandas as pd
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

### IMPORTING DATASETS
## We create a function to import the Turbie parameters: 
def load_turbie_parameters(filepath:str):
    """Load turbie parameters""" 
    with open(filepath, 'rt') as t_parameters:
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
    return m_blade, m_nacelle, m_tower, m_hub, c1, c2, k1, k2, D_rotor, rho

## We create a function to import CT parameters: 
def load_CT_data(filepath:str): 
    """Load thrust coefficient data."""
    ct_df = pd.read_csv(filepath, sep='\t')
    return ct_df

## We create a function to import Wind Profile Data for a given mean wind and turbulence: 
def load_wind_data(V: float, TI_cat: float):   
    """Load wind data based on assumed Wind Speed (V) and Turbulence (TI)""" 
    wind_df = pd.read_csv(f'inputs/wind_files/wind_TI_{TI_cat}/wind_{V}_ms_TI_{TI_cat}.txt', sep='\t')
    
    #We convert the wind data to an array
    wind_array = wind_df[['Time(s)', 'V(m/s)']].to_numpy()
    
    return wind_array

### BUILDING MATRICES
##  We create a function that builds M, C, and K matrices
def build_matrices(m_blade, m_nacelle, m_tower, m_hub, c1, c2, k1, k2):
    """Build turbie system matrices""" 
    
    m1 = 3*m_blade # mass of the three blades
    m2 = m_nacelle + m_hub + m_tower #combined mass of nacelle, hub and tower

    # We define the Mass Matrix (M), the Damping Matrix (C), and the Stiffness Matrix (K) for the Turbie system:
    M = np.array([[m1, 0],[0, m2]])
    C = np.array([[c1, -c1],[-c1, c1 + c2]])
    K = np.array([[k1, -k1],[-k1, k1 + k2]])
    
    return M, C, K

### BUILDING A FUNCTION TO DETERMINE Y'(T)
def calc_yprime(t, y, M, C, K, D_rotor, rho, ct_df, wind_array):
    x1, x2, x1_dot, x2_dot = y

    # We Build interpolation function for wind speed over time
    wind_time = wind_array[:, 0]
    wind_speed = wind_array[:, 1]
    u_t_func = interp1d(wind_time, wind_speed, kind='linear', fill_value='extrapolate')

    # We define u_t as the instantaneous wind speed at any time t
    u_t = float(u_t_func(t))

    ## We define aerodynamic parameters not defined in parameters dataset
    radius = D_rotor/2
    SW_A = radius**2*math.pi #Swept Area
    U = np.mean(wind_array[:, 1]) # We determine the mean wind speed

    # We build a linear interpolation function for C_t
    Ct_linear = interp1d(ct_df['# V(m/s)'], ct_df['CT'], kind='linear', fill_value='extrapolate') 
    C_t = Ct_linear(U) #We estimate C_t based on our linear interpolation function
    
    # We calculate aerodynamic force on mass 1
    v_rel = u_t - x1_dot
    f_aero = 1/2 * rho * C_t * SW_A * v_rel * abs(v_rel)
    
    # We define no external force on mass 2
    f_2 = 0
        
    # We define current positions and velocities
    x = np.array([x1, x2])
    x_dot = np.array([x1_dot, x2_dot])
    F_t = np.array([f_aero, f_2])
    
    # We calculate Accelerations
    Minv = np.linalg.inv(M) #The inverse of the mass matrix
    x_ddot = Minv @ (F_t - C @ x_dot - K @ x)

    #We define the system martrix (vector A) and input matrix (vector B)
    N = 2 # Number of freedom
    ZeroN = np.zeros((N, N))
    I = np.eye(N)
    A = np.block([[ZeroN, I], [-Minv@K, -Minv@C]])
    B = np.concatenate((np.zeros((2,)), Minv @ F_t))

    # We define y_prime (y')
    y_prime = A@y + B
    
    return y_prime




