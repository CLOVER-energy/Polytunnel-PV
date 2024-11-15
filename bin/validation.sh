########################################################################################
# validation.sh - Script for executing Polytunnel-PV validation a job on the HPC.      #
#                                                                                      #
# Author: Ben Winchester                                                               #
# Copyright: Ben Winchester, 2024                                                      #
# Date created: 15/11/2024                                                             #
#                                                                                      #
# For more information, please email:                                                  #
#   benedict.winchester@gmail.com                                                      #
########################################################################################
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=8Gb
#PBS -N polytunnel-pv-validation

python -m src.polytunnelpv --scenario circular_polysolar_kent_bypassed_2 \
    --operating-mode validation --timestamps-file hourly_december_data_martyn.csv
