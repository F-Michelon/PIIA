#!/bin/bash

# Spécifiez le chemin du fichier à surveiller
FILE="/home/acazaudumec/PIIA7/PIIA/RESULTATS/K_CLASSES_WITH_ZEROS/TESTBASH.txt"

# Commande d'exécution du programme
clingo /home/acazaudumec/PIIA7/PIIA/SCRIPTS/v10.2_no_parents_k_classes.lp /home/acazaudumec/PIIA7/PIIA/DONNEES/toy_datasets/B_3_classes_asp_data.lp > "$FILE" &

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

