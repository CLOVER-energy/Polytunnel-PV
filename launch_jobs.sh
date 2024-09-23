#!/bin/bash
#PBS -J 1-10
#PBS -lwalltime=8:00:00
#PBS -lselect=1:ncpus=64:mem=50000Mb

echo -e "HPC array script executed"

# Load the anaconda environment
module load anaconda3/personal
# Ensure conda is initialized
source /rds/general/user/bs921/home/anaconda3/etc/profile.d/conda.sh

# Activate the environment
conda activate myenv

cd $PBS_O_WORKDIR

# Install requirements
pip install --user -r requirements.txt

echo -e "Running HPC python script."
if python run_simulation.py; then
    echo -e "Polytunnel MPP Simulation successfully run."
else
    echo -e "FAILED. See logs for details."
fi
