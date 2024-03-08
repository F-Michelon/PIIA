#!/bin/bash

#SBATCH --job-name=caspo_big

#SBATCH --time=0-24:00:00       # Time limit
#SBATCH --nodes=1               # How many nodes to run on
#SBATCH --ntasks=1              # How many tasks per node
#SBATCH --cpus-per-task=1       # Number of CPUs per task
#SBATCH --mem-per-cpu=100g       # Memory per CPU
#SBATCH --qos=short             # priority/quality of service

caspo --out ../../RESULTATS/out_late learn ../../DONNEES/reduced_transformed_pkn.sif ../../DONNEES/dataset_late_midas.csv 10 --fit 0.001
caspo --out ../../RESULTATS/out_late/viz visualize --pkn ../../DONNEES/reduced_transformed_pkn.sif --networks ../../RESULTATS/out_late/networks.csv --setup ../../DONNEES/setup.json

caspo --out ../../RESULTATS/out_medium learn ../../DONNEES/reduced_transformed_pkn.sif ../../DONNEES/dataset_medium_midas.csv 10 --fit 0.001
caspo --out ../../RESULTATS/out_medium/viz visualize --pkn ../../DONNEES/reduced_transformed_pkn.sif --networks ../../RESULTATS/out_medium/networks.csv --setup ../../DONNEES/setup.json
