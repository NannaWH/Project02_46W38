##############################################################################################
## PROJECT 2: SIMULATION OF TURBIE ##########################################################
##############################################################################################

# This scripts calls the function in the turbie module to determine how the blade and 
# tower deflections vary with the wind speed and turbulence intensity 

"NB: LOAD THE DATA IS AS DATA INSTEAD OF TEXT"

"I need to loop over the different files using f-string to load and save the results"
"""for ws in [3, 4, 5]:
     file_name = f"results_for_ws_{ws}.txt"
     print(file_name) """

"""We need the forloop to run the simulation looping over different wind speeds"""
"""We have to simulate for each of the wind farms / data set"""

"Check for NBs before handing in"

#We import relevant packages
import itertools
import numpy as np
import pandas as pd
from pathlib import Path
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import itertools

#Variables that needs to come from the function statement - What we should be looping over
V = 20 # Wind speed: V(m/s)
TI_cat = 0.05 #TI category = 0.1, 0.05 or 0.15

#Folder for output
output_folder = Path("output")
output_folder.mkdir(parents=True, exist_ok=True)

#Defining the different wind speeds (V) and Turbilence categories (TI)
Vs = list(range(4, 26))
TI_cats = [0.05, 0.10, 0.15]

for V, TI_cat in itertools.product(Vs, TI_cats):

    #We import functions from "turbine_mod.py"
    from turbie_mod import (load_turbie_parameters, 
                            load_CT_data, 
                            load_wind_data, 
                            build_matrices, 
                            calc_yprime)

    # We load data into the model
    m_blade, m_nacelle, m_tower, m_hub, c1, c2, k1, k2, D_rotor, rho = load_turbie_parameters('inputs/turbie_parameters.txt')
    ct_df = load_CT_data('inputs/CT.txt')
    wind_array = load_wind_data(V=V, TI_cat=TI_cat)

    # We call the function to build the M, C, K matrices
    M, C, K = build_matrices(m_blade, m_nacelle, m_tower, m_hub, c1, c2, k1, k2)

    # We set the timespan
    t_span = (0, 659.990) 
    t_eval = np.linspace(*t_span, 66000)

    # We set initial conditions for y'(t)
    y0 = [0,0,0,0]

    res = solve_ivp(calc_yprime, t_span, y0, t_eval=t_eval, args=(M, C, K, D_rotor, rho, ct_df, wind_array))

    # Convert results to DataFrame
    results_df = pd.DataFrame(res.y.T, columns=[f"y{i}" for i in range(res.y.shape[0])])
    results_df.insert(0, "time", res.t)

    # Save results
    file_name = output_folder / f"results_for_ws_V={V}_TI={TI_cat}.txt"
    results_df.to_csv(file_name, index=False, sep="\t")

    print(f"Saved results for V={V}, TI={TI_cat} to {file_name}")


#Plotting the Blade and Tower Displacements: 
plt.plot(res.t, res.y[0], label='x1 (Blade Displacement)')
plt.plot(res.t, res.y[1], label='x2 (Tower Displacement)')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.title('2-DOF Mass-Spring-Damper with Aerodynamic Force')
plt.legend()
plt.grid(True)
plt.show()

#Plotting the Relative Blade and Tower Displacements: 
plt.plot(res.t, res.y[0]-res.y[1], label='x1-x2 (Relative Displacement)')
plt.plot(res.t, res.y[1], label='x2 (Tower Displacement)')
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.title('2-DOF Mass-Spring-Damper with Aerodynamic Force')
plt.legend()
plt.grid(True)
plt.show()

#Plot the C_t curve
plt.plot(ct_df['# V(m/s)'], ct_df['CT'])
plt.xlabel('Wind Speed (m/s)')
plt.ylabel('CT')
plt.show()

#Calculate the means and standard deviations
