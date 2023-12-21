import pandas as pd

def enlever_les_zeros(name_file):
    #On récupère le fichier csv et les données qui le concerne
    dataframe = pd.read_csv(name_file)
    long_avant_modif = len(dataframe.columns)
    premieres_colonnes = dataframe[dataframe.columns[0:8].values]
    #print(premieres_colonnes)
    
    #On enlève les colonnes qui ne sont pas à traiter
    for elt in dataframe.columns[0:8].values:
        dataframe = dataframe.drop(columns=elt)
    
    #On fait la somme cumulative sur l'axe des ordonnées
    cumsum = dataframe.cumsum(axis = 0).iloc[-1]
    liste_colonne_a_enlever = []
    for i in range(len(cumsum)):
        if cumsum[i] == 0:
            liste_colonne_a_enlever.append(dataframe.columns[i])
    for name in liste_colonne_a_enlever:
        dataframe = dataframe.drop(columns = name)
    dataframe = pd.concat([premieres_colonnes,dataframe],axis = 1)
    dataframe.to_csv(name_file[:-4]+f"_traite_{len(liste_colonne_a_enlever)}.csv", index=False)
    print(len(liste_colonne_a_enlever)/long_avant_modif)
    return dataframe

enlever_les_zeros("../DONNEES/toy_datasets/C_3_classes.csv")
enlever_les_zeros("../DONNEES/toy_datasets/C_4_classes.csv")