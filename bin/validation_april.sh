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
#PBS -N ppv-validation

eval "$(~/anaconda3/bin/conda shell.bash hook)"
source activate py312
conda activate py312

# Change to the submission directory
cd $PBS_O_WORKDIR

python -m src.polytunnelpv -s circular_polysolar_kent_bypassed_2 \
    -o validation -t hourly_april_map_data_martyn_rehourd.csv \
    -st {START_TIME} -i 24 -wds hadlow1625
