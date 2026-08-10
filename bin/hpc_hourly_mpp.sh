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
source activate py312

# Change to the submission directory
cd $PBS_O_WORKDIR

echo "PBS ARRAY INDEX $PBS_ARRAY_INDEX"
START_TIME=$(($PBS_ARRAY_INDEX * 24 + 7224))
echo "Launching with start time $START_TIME"

python -u -m src.polytunnelpv -wds hadlow1625 -o hourly_mpp \
    -s circular_polysolar_kent_bypassed_2 -st $START_TIME -i 24 --on-hpc

