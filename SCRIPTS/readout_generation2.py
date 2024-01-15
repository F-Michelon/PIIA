'''
Ce script génère des readouts réels (non fictif) à partir du fichier raw_mtx.csv qui est trop volumineux pour être partagé dans le git
A chacun de l'ajouter dans son répertoire DONNEES mais sans le versionné sur git ! Donc fini les "git add --all". 
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder, normalize


PATH = "../DONNEES/real_datasets/raw_mtx.csv"

D = pd.read_csv(PATH, usecols=np.concatenate((np.arange(8),np.random.randint(9,30000, 10)))).drop(columns={'Stage', 'Dataset', 'Pseudotime', 'Embryo', 'EmbryoDay', 'clusterEmbryoCells', 'Name'})

labelEncode = LabelEncoder()
D['classes'] = labelEncode.fit_transform(D['clusterUmap'])
D = D.drop(columns='clusterUmap')

# On défini des vecteurs booléens fictifs
readouts = (2/np.pi) * np.arctan(D[D.keys()[:10]].values)

# On crée des vecteurs booléens de taille 10 (4 types différents présents dans l'ensemble des classes)

vect_bool = np.zeros((len(D),10))
vect_bool_types = np.random.randint(low=0, high=2, size=(4, 10))
for i in D['classes'].unique():
    index_i = D[D['classes'] == i].index
    for k, idx in enumerate(index_i):
        # On fait en sorte que chaque vecteur booléens soit présent dans chaque classe
        if k < (len(index_i) // 4):
            chosen_vect_bool = 0
        elif k < (len(index_i) // 2) and k >= (len(index_i) // 4):
            chosen_vect_bool = 1
        elif k < (len(index_i) // 4) * 3 and k >= (len(index_i) // 2):
            chosen_vect_bool = 2
        else:
            chosen_vect_bool = 3
        vect_bool[idx] = vect_bool_types[chosen_vect_bool]
data = np.concatenate((vect_bool, readouts), axis=1)

D_new = pd.DataFrame(data)
D_new['classes'] = D['classes']
D_new.to_csv("../DONNEES/toy_datasets/readout_fictifs_raw_mtw.csv", index=False)

print(D_new.head(5))
for i in range(10):
    print(np.sum(readouts[i]))
print(D_new['classes'].value_counts())
print(D_new[D_new.keys()[:10]])