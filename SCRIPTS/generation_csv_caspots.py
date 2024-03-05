"""
Script permettant de traiter le fichier : PIIA/RESULTATS/CASPO/2023_PKN_earlyAndMediumAndLate_traite_2_discrimination.csv
Principe du traitement :
Le fichier comprend 3 colonnes, une pour chaque classe (early, medium et late)
chaque ligne contient un nombre qui correspond à une ligne (= une cellule) dans le fichier : PIIA/DONNEES/real_datasets/raw_mtx.csv
il faut donc d'abord récupérer les valeurs des cellules correspondant à ces nombres
pour ensuite les regrouper par vecteur booléen et enfin leur donner une valeur temporelle dans un fichier midas
"""

import pandas as pd
import json
import numpy as np

PATH = "../DONNEES/real_datasets/2023_PKN_earlyAndMediumAndLate_traite_2_readouts.csv"

#with open('../DONNEES/setup.json') as setup_file:
#    setup = json.load(setup_file)
#with open('../DONNEES/perturbations.json') as perturbation_file:
#    perturbation = json.load(perturbation_file)
    
D = pd.read_csv(PATH)
input_genes = D.keys()[1:11]
readouts_genes = D.keys()[11:-1]

interesting_cells = pd.read_csv("../RESULTATS/CASPO/2023_PKN_earlyAndMediumAndLate_traite_2_discrimination.csv")

D = D.loc[list(np.array(interesting_cells).flatten())]

dico_vect_bool = dict()

for i in D.index :
    
    #On récupère le vecteur booléen de la ligne 
    vect_bool = ""
    for gene in input_genes :
        vect_bool = vect_bool + str(D[gene].loc[i])
    
    #On regarde dans quelle classe se trouve la cellule
    classe = D["classes"].loc[i]
    
    #On teste s'il on l'a déjà en mémoire
    if vect_bool in dico_vect_bool.keys():
        #Si oui on l'ajoute 
        tmp_dict = dico_vect_bool[vect_bool]
        tmp_dict[classe] = i
    else :
        #Sinon on crée dans le dico
        new_dict = dict()
        new_dict[classe] = i
        dico_vect_bool[vect_bool] = new_dict

#Création du dictionnaire qui va servir à la création du fichier csv pour caspots
Dataset4Caspots = dict()

#On commence par les noms des colonnes
Dataset4Caspots["TR:mock:CellLine"] = [1]*len(D)
for key in input_genes :
    Dataset4Caspots["TR:" + key] = []
for key in readouts_genes :
    Dataset4Caspots["DA:" + key] = []
for key in readouts_genes :
    Dataset4Caspots["DV:" + key] = []

#Le dico de vecteur booléen contient tous les vecteurs booléens différents
#pour chacun il y a n cellules (une de chaque classe)
#il faut rajouter dans le fichier d'abord tous les vecteurs d'une même classe
#puis les mêmes vecteurs dans le même order pour la classe suivante
for numero_classe in range(D["classes"].max()+1) :
    for vect_bool in dico_vect_bool :
        
        #D'abord les gènes input (ie vecteurs booléens)
        for key,val_bool in zip(input_genes,vect_bool) :
            Dataset4Caspots["TR:" + key].append(val_bool)
        
        #Ensuite la classe pour les gènes de readouts
        for key in readouts_genes :
            Dataset4Caspots["DA:" + key].append(numero_classe)
        
        #Enfin les valeurs des readouts pour les gènes de readouts
        numero_cell = dico_vect_bool[vect_bool][numero_classe]
        cell_data = D[readouts_genes].loc[numero_cell]
        for key,val_readout in zip(readouts_genes,cell_data):
            Dataset4Caspots["DV:" + key].append(val_readout)

Dataset4Caspots = pd.DataFrame(Dataset4Caspots)
Dataset4Caspots.to_csv("../DONNEES/Dataset4Caspots.csv", index=False)



