#!/bin/bash
#PBS -J 1-10  
#PBS -l walltime=03:00:00  
#PBS -l select=1:ncpus=8:mem=12gb  
#PBS -q v1_general72  

echo -e "HPC array script executed"

# Load the anaconda environment
module load anaconda3/personal
source activate py310

cd $PBS_O_WORKDIR

# Install requirements
pip install --user -r requirements.txt

echo -e "Running HPC python script."
if python run_simulation.py; then
    echo -e "Polytunnel MPP Simulation successfully run."
else
    echo -e "FAILED. See logs for details."
fi
