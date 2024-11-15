########################################################################################
# hpc.sh - Runs the Polytunnel-PV model in an HPC array job.                           #
#                                                                                      #
# Author(s): Ben Winchester                                                            #
# Copyright: Ben Winchester, 2024                                                      #
# Date created: 15/11/2024                                                             #
# License: Open source                                                                 #
# Most recent update: 15/11/2024                                                       #
#                                                                                      #
# For more information, please email:                                                  #
#     benedict.winchester@gmail.com                                                    #
########################################################################################

# Depending on the environmental variable, run the appropriate HPC job.
module load anaconda3/personal
eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate py310

# Change to the submission directory
cd $PBS_O_WORKDIR

# Determine the scenario to run
SCENARIO="$(head -n $PBS_ARRAY_INDEX scenarios.txt | tail -n 1)"

python -m src.polytunnelpv --scenario $SCENARIO -st 0 -i 8760 --operating-mode hourly_mpp \
    --timestamps-file hourly_december_data_martyn.csv
