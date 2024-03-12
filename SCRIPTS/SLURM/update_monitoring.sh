#!/bin/bash

#SBATCH --job-name=myjobPKN2023bash   # Name for your job
#SBATCH --comment="Run My Job"  # Comment for your job

#SBATCH --time=1-00:00:00       # Time limit
#SBATCH --qos=long              # Length of task
#SBATCH --nodes=1               # How many nodes to run on
#SBATCH --ntasks=1              # How many tasks per node
#SBATCH --cpus-per-task=2       # Number of CPUs per task
#SBATCH --mem-per-cpu=100g       # Memory per CPU

cd ../..

# Spécifiez le chemin du fichier à surveiller
FILE="./RESULTATS/RESULTATS_NAUTILUS/K_CLASSES_WITHOUT_ZEROS/resultat_2023_PKN_earlyAndMediumAndLate_traite_2_v10.2_no_parents_k_classes-v2test.txt"

# Commande d'exécution du programme
clingo -n 0 ./DONNEES/real_datasets/2023_PKN_earlyAndMediumAndLate_traite_2_asp_data.lp SCRIPTS/v10.2_no_parents_k_classes-v2.lp > "$FILE" &

# Boucle infinie pour surveiller les modifications du fichier
while true; do
    # Vérifier si le fichier a été modifié depuis la dernière sauvegarde
    if [ -f "$FILE" ]; then
        # Calculer le hachage MD5 du fichier
        MD5SUM_OLD=$(md5sum "$FILE")
	sleep 0.001
        MD5SUM_NEW=$(md5sum "$FILE")

        # Vérifier si le fichier a été modifié
        if [ "$MD5SUM_OLD" != "$MD5SUM_NEW" ]; then
            # Créer une copie de sauvegarde avec un timestamp
            TIMESTAMP=$(date +"%Y%m%d_%H%M%S%N")
            cp "$FILE" "${FILE}_backup_$TIMESTAMP"
            echo "Sauvegarde créée : ${FILE}_backup_$TIMESTAMP"
        fi
    fi

    # Attendre un certain temps avant de vérifier à nouveau
    sleep 0.000001
done

