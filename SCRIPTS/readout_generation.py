# Ce script génère des readouts fictifs pour trois classes cellulaires

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder

PATH = "../DONNEES/toy_datasets/D_3_classes.csv"

D = pd.read_csv(PATH).drop(columns={'Stage', 'Dataset', 'Pseudotime', 'Embryo', 'EmbryoDay', 'clusterEmbryoCells', 'Name'})

labelEncode = LabelEncoder()
D['classe'] = labelEncode.fit_transform(D['clusterUmap'])
D = D.drop(columns='clusterUmap')

# Création de gènes fictifs
N_readouts = 20 # nombre de gène
N_cells_per_pseudoperturbation = 30
Nb_classes = 3

readouts_list = [str(i) for i in range(N_readouts)]

# cells chosen to represent a behaviour for pseudo perturbation
random_index = np.random.randint(0, 100, 4)

# On crée un nouveau dataframe 
columns_list = D.keys()[:-1]

data = []
for i in range(Nb_classes):
    for j in random_index:
        for k in range(N_cells_per_pseudoperturbation):
            data.append(np.concatenate((D.iloc[j].values[:-1], np.random.random(size=N_readouts), [list(random_index).index(j)])))

D_new = pd.DataFrame(np.array(data))
D_new = D_new.rename(columns={140 : 'vect_bool_alike'})
D_new['classes'] = [0 for i in range(120)] + [1 for i in range(120)] + [2 for i in range(120)]
D_new.to_csv("../DONNEES/toy_datasets/readout_fictifs_D_3", index=False)