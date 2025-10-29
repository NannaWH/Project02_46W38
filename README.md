# Project02_46W38
This repository includes the scripts used to solve project 2 in the course Scientific Programming in Wind Energy at DTU, Fall 2025.

The scripts simulates Turbie, a simple two-degree-of-freedom (2DOF) system based of the DTU 10 MW Reference WindTurbine. She is equivalent to the forced 2DOF mass-spring-damper system shownbelow, which is defined by a mass matrix, stiffness matrix, and damping matrix.

The repository consist of two python files, an input folder, and an output folder:
    iputs: this is the input folder that contains the data inputs used in the simulation of turbie.
    output: this is the output that contain the output of each turbie simulation as well as graphs used to interpret the results. 
    turbie_mod.py: defines the functions that are used in main to simulate turbie.
    main: runs the turbie simulation for each of the wind and turbulence cases, and plots graphs with some of the main results.

The simulations of turbie shows that the relative displacement of blades and towers fluctuates almost equally with the wind speeds which makes sense as the fluctuations in the aerodynamic forcing on the blades is determined by fluctuations in wind only, everything else in the aerodynamic forcing equation is constant for a given simulation.

From the model we find that the relationship between mean displacement and mean wind speed is almost unaffected by turbulence while the displacement's standard deviations are propotional to the turbulence, i.e., higher turbulence results in a higher standard deviation of displacement. This makes sense as the turbulence is defined as short-term fluctuations around the mean wind speed and the size of the displacement of the blades is driven by the wind speeds - thereby will the tower displacement also be affected by the wind speeds only through the displacement of the blades. The standard deviations of displacement for both the blades and tower increases as more turbulence results in higher fluctuations in the displacement around its mean value. 
We also find that both the means and standard deviations of the blades and tower displacements a non-linear functions with a maximum displacement and standard deviation around a wind speed of 11 m/s. This is because the torque on the blades affecting the aerodynamic forcing on the blades decreases dramatically for wind speeds above 11 m/s which is probably due to controlling systems in the wind turbine pitching the blades for wind speeds that are higher than optimal for the turbine.

One things to note is that we have simulated turbie based on the total mass off the tower/hub/nacelle and the blades, respectively. Thereby, ee assume that it is all of the wind blade mass and tower/hub/nacelle mass that will be displaced by wind speeds which would not be true in the real world. It will only be the top part of the tower incl. hub and nacelle and the flexible part of the blades that will get into movement and affect the displacement.

