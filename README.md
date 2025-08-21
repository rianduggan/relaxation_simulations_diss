# relaxation_simulations_diss
This repository contains scripts accompanying my MSc Dissertation. For information on how to carry out the experiments, please read this README file.

## Section 2 - Lotka-Volterra and Hare-Lynx
The following two scripts deal with the simulations in Section 2. One can run them to reproduce the results shown, or adjust parameters to experiment with different initial conditions.

### LV_Sim.py
Run the script and adjust the initial conditions to explore different predator-prey orbits in the Lotka-Volterra system.
A time series is also plotted to see the interaction of populations over time.

### HareLynx_Sim.py
Run the script to plot hare and lynx populations over time, and adjust the year filter to explore different predator–prey cycles.

## Section 3 - Aziz Model
The following two scripts are related to the figures created in the Numerical Simulation (Section 3.10). 

### Aziz_Relax_Sim.py
This script is initially set to simulation relaxation oscillation and an accompanying time series. Following the figure captions in Section 3.10, one can reproduce all of the figures, or feel free to experiment with other conditions. 

### Aziz_Canard_Sim.py
This code is exactly the same as Aziz_Relax_Sim.py, except it has an zoomed inset to show the emergence of canards at the fold point. One can adjust parameters to see how sensitive these solutions really are.

## Section 4 - Wang Model
The following two scripts are related to the figures created in the Numerical Simulation (Section 4.9). 

### Wang_Relax_Sim.py
This script is initially set to simulation relaxation oscillation and an accompanying time series. Following the figure captions in Section 4.9, one can reproduce all of the figures, or feel free to experiment with other conditions. 
We also note here that this code is used for figures in Section 4.11 relating to fear. 
All varibles are defined at the beginning of the code.

### Aziz_Canard_Sim.py
This code is exactly the same as Wang_Relax_Sim.py, except it has an zoomed inset to show the emergence of canards at the fold point. One can adjust parameters to see how sensitive these solutions really are.
