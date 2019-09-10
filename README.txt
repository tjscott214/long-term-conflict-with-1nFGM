### Supplemental code for Scott and Queller (2019) Ecology and Evolution ###

Python files in this folder simulate a one dimensional Fisher’s geometric model with conflict and simulates abiotic change with fixations from conflict simulations.

FILES

FGMconflict.py - simulates conflict and standard geometric model
FGMabiotic.py - takes mut.csv output from FGMconflict.py and simulates adaptation to abiotic change of equal magnitude to change due to conflict
run_simulations.sh - runs example simulations and plots output
plot_results.R - plots the output of example simulations

RUNNING THE CODE

Below is example command line code for running both programs:
python3 FGMconflict.py -x 1 -z 500 -c 0.2
python3 FGMabiotic.py -x 1 -z 500 -c 0.2 -mut mut.csv
This code runs conflict, standard, and abiotic change simulations with 500 iterations, average mutations sizes of 0.1 (default), and conflict equal to 0.2

To run these example simulations and plot the results:
sh run_simulations.sh

For parameter options, use the -h flag
python3 FGMconflict.py -h
python3 FGMabiotic.py -h 

CONTACT

For questions or suggestions, contact Trey J Scott at tjscott@wustl.edu


