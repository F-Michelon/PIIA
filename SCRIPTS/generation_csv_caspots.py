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

with open('../DONNEES/setup.json') as setup_file:
    setup = json.load(setup_file)
with open('../DONNEES/perturbations.json') as perturbation_file:
    perturbation = json.load(perturbation_file)
    
D = pd.read_csv(PATH, usecols=['Name', 'clusterUmap'] + setup["stimuli"] + setup["inhibitors"] + setup['readouts'])[['Name', 'clusterUmap'] + setup["stimuli"] + setup["inhibitors"] + setup['readouts']]

interesting_cells = pd.read_csv("../RESULTATS/CASPO/2023_PKN_earlyAndMediumAndLate_traite_2_discrimination.csv")

D = D.iloc[list(np.array(interesting_cells).flatten())]
D = D.sort_values(by=['clusterUmap'], ascending=False)

D[setup["stimuli"] + setup["inhibitors"]].replace(1,0)
D[setup["stimuli"] + setup["inhibitors"]] = D[setup["stimuli"] + setup["inhibitors"]].mask(D[setup["stimuli"] + setup["inhibitors"]]>1,1)
D[setup["readouts"]] = (2/np.pi) * np.arctan(D[setup['readouts']])



dico_vect_bool = dict()

for i in D.index :
    
    #On récupère le vecteur booléen de la ligne 
    vect_bool = ""
    for gene in setup["stimuli"] + setup["inhibitors"] :
        vect_bool = vect_bool + str(D[gene].loc[i])
    
    #On regarde dans quelle classe se trouve la cellule
    classe = D["clusterUmap"].loc[i]
    
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
        
    

    
        
    

dataset_caspots = dict()



