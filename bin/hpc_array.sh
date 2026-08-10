########################################################################################
# hpc.sh - Script for executing Polytunnel-PV as an array job on the HPC.              #
#                                                                                      #
# Author: Ben Winchester                                                               #
# Copyright: Ben Winchester, 2024                                                      #
# Date created: 15/11/2024                                                             #
#                                                                                      #
# For more information, please email:                                                  #
#   benedict.winchester@gmail.com                                                      #
########################################################################################
#PBS -J 1-2527
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=8Gb
#PBS -N ppv

echo -e "HPC array script executed"

# Load the anaconda environment
eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate py312

cd $PBS_O_WORKDIR

echo -e "Running Polytunnel-PV launch script."
if python -m ./bin/hpc.sh ; then
    echo -e "Polytunnel-PV model successfully run."
else
    echo -e "FAILED. See logs for details."
fi
