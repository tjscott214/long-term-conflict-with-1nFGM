#!/bin/sh


python3 FGMconflict.py -x 1 -z 500 -c 0.2
python3 FGMabiotic.py -x 1 -z 500 -c 0.2 -mut mut.csv
Rscript plot_results.R FGM_data.csv abiotic_data.csv 0.1 0.2

exit 0