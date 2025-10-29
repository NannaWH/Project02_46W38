##############################################################################
## PROJECT 2: SIMULATION OF TURBIE ## ##############################################################################

# This scripts calls the function in the turbie module to determine how the blade and tower deflections vary with the wind speed and turbulence intensity 

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

#We import functions from "turbine_mod.py"
from turbie_mod import (load_turbie_parameters, 
                        load_CT_data, 
                        load_wind_data, 
                        build_matrices, 
                        calc_yprime)


# We load the general data into the model
if __name__ == "__main__": 
    m_blade, m_nacelle, m_tower, m_hub, c1, c2, k1, k2, D_rotor, rho = load_turbie_parameters('inputs/turbie_parameters.txt')
    ct_df = load_CT_data('inputs/CT.txt')

    # We call the function to build the M, C, K matrices
    M, C, K = build_matrices(m_blade, m_nacelle, m_tower, m_hub, c1, c2, k1, k2)

#We define a folder for output from the Turbie model
output_folder = Path("output")
output_folder.mkdir(parents=True, exist_ok=True)

#We define a dataframe for the mean and standard deviations data
mean_std_data = []

## RUNNING THE TURBIE MODEL FOR ALL DATASETS

#We define the different wind speeds (V) and Turbulence categories (TI)
Vs = list(range(4, 26))
TI_cats = [0.05, 0.10, 0.15]

# We run the turbie model for each of of wind speed and turbulence datasets
if __name__ == "__main__": 
    for V, TI_cat in itertools.product(Vs, TI_cats):

        #We load the wind and turbulence specific datasets into the model
        wind_array = load_wind_data(V=V, TI_cat=TI_cat)

        # We set the timespan
        t_span = (0.000, 659.990) 
        t_eval = np.linspace(*t_span, 66000)

        # We set initial conditions for y'(t)
        y0 = [0,0,0,0]

        #We run the solve_ivp function to solve for displacements
        res = solve_ivp(calc_yprime, t_span, y0, t_eval=t_eval, args=(M, C, K, D_rotor, rho, ct_df, wind_array))

        # Convert results to DataFrame
        results_df = pd.DataFrame(res.y.T, columns=[f"y{i}" for i in range(res.y.shape[0])])
        results_df.insert(0, "Time(s)", res.t)
        # We remove the results of the first minute of the data to remove the effect of the initial value which is unknown
        results_df = results_df[6000:]

        # Save results
        file_name = output_folder / f"results_for_ws_V={V}_TI={TI_cat}.txt"
        results_df.to_csv(file_name, index=False, sep="\t")
        print(f"Saved results for V={V}, TI={TI_cat} to {file_name}")

        #Calculate mean and standard deviations of displacement
        # and save it to the "mean_st_data" dataframe
        mean_std_data.append({
        "Turbulence" : TI_cat,
        "Wind_speed" : V,
        "mean_blades" : np.mean(results_df["y0"]),
        "std_blades" : np.std(results_df["y0"]),
        "mean_tower" : np.mean(results_df["y1"]),
        "std_tower" : np.std(results_df["y1"])
        })

#Saving the mean and std dataframe to a file in the output folder
summary_df = pd.DataFrame(mean_std_data)
file_name_m_std = output_folder / f"mean_std_data.txt"
summary_df.to_csv(file_name_m_std, index=False, sep="\t")

### PLOTTING THE RESULTS
## Plotting the Relative Blade and Tower Displacements for a given V and TI: 
V= 10
TI_cat = 0.10

#Weload the data
displacement_data = pd.read_csv(f"output/results_for_ws_V={V}_TI={TI_cat}.txt", sep='\t')

wind_data = pd.read_csv(f'inputs/wind_files/wind_TI_{TI_cat}/wind_{V}_ms_TI_{TI_cat}.txt', sep='\t')

#We merge the wind speeds to the data
merged_data = pd.merge(displacement_data, wind_data, on = 'Time(s)', how ='left')

fig, ax1 = plt.subplots()

# Plot the first two displacement series on ax1
line1, = ax1.plot(merged_data["Time(s)"], merged_data["y0"] - merged_data["y1"],
                   'r-', label='x1 - x2 (Relative Displacement)')
line2, = ax1.plot(merged_data["Time(s)"], merged_data["y1"],
                   'g-', label='x2 (Tower Displacement)')

# Plot wind speed on the secondary y-axis
ax2 = ax1.twinx()
line3, = ax2.plot(merged_data["Time(s)"], merged_data["V(m/s)"],
                  'b--', label='V(m/s)')

# Axis labels and title
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Displacement (m)')
ax2.set_ylabel('V(m/s)')
plt.title(f"2-DOF Mass-Spring-Damper - Relative Displacements, V={V} TI={TI_cat}")

# Combine legends from both axes
lines = [line1, line2, line3]
labels = [line.get_label() for line in lines]
ax1.legend(lines, labels, loc='upper right')

ax1.grid(True)
plt.tight_layout()

# We save the plot
plt.savefig("output/Relative_Displacements.png")


## Plot the means and standard deviations of displacements
# We load the data
summary_df = pd.read_csv("output/mean_std_data.txt", sep="\t")

# We prepare the plot to be a 2x2
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)

# We define turbulence categories and colors
turb_values = summary_df["Turbulence"].unique()
colors = plt.cm.tab10.colors  # up to 10 distinct colors

# We plot each graph in the 2x2 graph frame
for t, c in zip(turb_values, colors):
    subset = summary_df[summary_df["Turbulence"] == t]
    axes[0, 0].plot(subset["Wind_speed"], subset["mean_blades"], marker="o", label=f"Turb={t}", color=c)
    axes[0, 1].plot(subset["Wind_speed"], subset["mean_tower"], marker="o", label=f"Turb={t}", color=c)
    axes[1, 0].plot(subset["Wind_speed"], subset["std_blades"], marker="o", label=f"Turb={t}", color=c)
    axes[1, 1].plot(subset["Wind_speed"], subset["std_tower"], marker="o", label=f"Turb={t}", color=c)

# We add titles, labels, and legends
axes[0, 0].set_title("Blade Displacement (mean) by Mean Wind Speed")
axes[0, 0].set_xlabel("Mean Wind Speed")
axes[0, 0].set_ylabel("Mean Displacement (Blades)")

axes[0, 1].set_title("Tower Displacement (mean) by Mean Wind Speed")
axes[0, 1].set_xlabel("Mean Wind Speed")
axes[0, 1].set_ylabel("Mean Displacement (Tower)")

axes[1, 0].set_title("Blade Displacement (std) by Wind Speed (SD)")
axes[1, 0].set_xlabel("Mean Wind Speed")
axes[1, 0].set_ylabel("Standard Deviation of Displacement (Blades)")

axes[1, 1].set_title("Tower Displacement (std) by Mean Wind Speed")
axes[1, 1].set_xlabel("Mean Wind Speed")
axes[1, 1].set_ylabel("Standard Deviation of Displacement (Tower)")

# We apply legend and grid to all subplots
for ax in axes.flat:
    ax.legend()
    ax.grid(True)

# We save the figure
plt.tight_layout()
plt.savefig("output/Mean_std.png", dpi=300)





