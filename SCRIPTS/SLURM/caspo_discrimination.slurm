#!/bin/bash

#SBATCH --job-name=caspo_big2

#SBATCH --time=1-00:00:00       # Time limit
#SBATCH --nodes=1               # How many nodes to run on
#SBATCH --ntasks=1              # How many tasks per node
#SBATCH --cpus-per-task=1       # Number of CPUs per task
#SBATCH --mem-per-cpu=100g       # Memory per CPU
#SBATCH --qos=long

cd ..
python readout_discrimination.py
