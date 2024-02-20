#!/bin/bash

#SBATCH --job-name=myjobPKN2023inputinter        # Name for your job
#SBATCH --comment="Run My Job"  # Comment for your job

#SBATCH --time=0-24:00:00       # Time limit
#SBATCH --nodes=1               # How many nodes to run on
#SBATCH --ntasks=1              # How many tasks per node
#SBATCH --cpus-per-task=2       # Number of CPUs per task
#SBATCH --mem-per-cpu=100g       # Memory per CPU

cd ../../

#Avec inter, input etc
#clingo ./DONNEES/toy_datasets/D_3_classes_traite_4_asp_data.lp SCRIPTS/v10.2_no_parents_k_classes-v2.lp > ./RESULTATS/K_CLASSES_WITHOUT_ZEROS/resultat_D_3_classes_traite_4_v10.2_no_parents_k_classes-v2.txt
#clingo DONNEES/toy_datasets/D_4_classes_traite_2_asp_data.lp SCRIPTS/v10.2_no_parents_k_classes-v2.lp > RESULTATS/K_CLASSES_WITHOUT_ZEROS/resultat_D_4_classes_traite_2_v10.2_no_parents_k_classes-v2.txt
#clingo DONNEES/toy_datasets/C_4_classes_traite_9_asp_data.lp SCRIPTS/v10.2_no_parents_k_classes-v2.lp > RESULTATS/K_CLASSES_WITHOUT_ZEROS/resultat_C_4_classes_traite_9_v10.2_no_parents_k_classes-v2.txt
clingo ./DONNEES/real_datasets/2023_PKN_earlyAndMediumAndLate_traite_2_asp_data.lp SCRIPTS/v10.2_no_parents_k_classes-v2.lp > ./RESULTATS/K_CLASSES_WITHOUT_ZEROS/resultat_2023_PKN_earlyAndMediumAndLate_traite_2_v10.2_no_parents_k_classes-v2.txt

#Programme 2023 Sans inter, input etc
#clingo ./DONNEES/toy_datasets/D_3_classes_traite_4_asp_data2.lp SCRIPTS/k-classes-v2.lp > ./RESULTATS/K_CLASSES_WITHOUT_ZEROS/resultat_D_3_classes_traite_4_k-classes-v2.txt
#clingo DONNEES/toy_datasets/D_4_classes_traite_2_asp_data2.lp SCRIPTS/k-classes-v2.lp > RESULTATS/K_CLASSES_WITHOUT_ZEROS/resultat_D_4_classes_traite_2_k-classes-v2.txt
#clingo DONNEES/toy_datasets/C_4_classes_traite_9_asp_data2.lp SCRIPTS/k-classes-v2.lp > RESULTATS/K_CLASSES_WITHOUT_ZEROS/resultat_C_4_classes_traite_9_k-classes-v2.txt

#Programme 2024 Sans inter, input etc
#clingo ./DONNEES/toy_datasets/D_3_classes_traite_4_asp_data2.lp SCRIPTS/v10.2_no_parents_k_classes.lp > ./RESULTATS/K_CLASSES_WITHOUT_ZEROS/resultat_D_3_classes_traite_4_v10.2_no_parents_k_classes.txt
#clingo DONNEES/toy_datasets/D_4_classes_traite_2_asp_data2.lp SCRIPTS/v10.2_no_parents_k_classes.lp > RESULTATS/K_CLASSES_WITHOUT_ZEROS/resultat_D_4_classes_traite_2_v10.2_no_parents_k_classes.txt
#clingo DONNEES/toy_datasets/C_4_classes_traite_9_asp_data2.lp SCRIPTS/v10.2_no_parents_k_classes.lp > RESULTATS/K_CLASSES_WITHOUT_ZEROS/resultat_C_4_classes_traite_9_v10.2_no_parents_k_classes.txt
