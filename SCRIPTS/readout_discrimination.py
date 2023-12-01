# Ce script développe des methodes de discrimination de readout pour du multi classe

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

long_vect_bool = 120

PATH = "../DONNEES/toy_datasets/readout_fictifs_D_3"

D = pd.read_csv(PATH)

# method 1
""" Choisir séquentiellement les readouts les plus différents
"""

# method 2
"""
"""

def find_same_vect_bool(df):
    """
    Fonction qui trouve les lignes du
    tableau pandas qui ont le même vecteur booléen
    
    Renvoie
    -------
    List_index_same_vect_bool : List[List[int]]
        Une liste contenant des sous listes.
        Chaque sous-liste correspond à un vecteur booléen différent
        Ces sous-listes contiennent les indexes (pour le dataframe)
        des cellules qui ont le même vecteur booléen
    List_bool_vect : List[List[bool]]
        Une liste contenant des sous listes.
        Chaque sous liste correspond à un vecteur booléen différent
    """
    N = len(df)
    List_bool_vect = []
    List_index_same_vect_bool = []
    for i in range(N):
        vect_bool = [[df.iloc[i][k] for k in range(long_vect_bool)]] #Le vecteur booléen de la ligne i du tableau
        if not(vect_bool in List_bool_vect): #S'il existe pas déja
            List_bool_vect.append(vect_bool) #On l'ajoute à la liste
            List_index_same_vect_bool.append([i]) #Et on ajoute son indice à la liste des numéros de chasue vecteurs boléens
        
        else : #S'il est dans la liste
            for j in range(len(List_bool_vect)): #On cherche quel vecteur c'est
                if List_bool_vect[j] == vect_bool: #Si c'est le vecteur n° j de la liste
                    List_index_same_vect_bool[j].append(i) #On ajoute son indice à la liste
    
    return List_index_same_vect_bool,List_bool_vect

# #Test
# test_find = find_same_vect_bool(D)

def define_3_cel_for_each_vect(df,long_vect_bool,List_index_same_vect_bool,List_bool_vect):
    """
    Fonction qui, pour chaque vecteur booléen différent,
    prend au hasard trois cellules qui partagent ce même vecteur booléen
    mais qui appartiennent à des classes différentes
    
    Renvoie
    -------
    List_3_cels_for_each_vect : List[int]
        La liste des indexes dans le dataframe des cellules choisies,
        il y en a 12 au total, 3 cellules pour chaque vecteurs
        et il y a 4 vecteurs différents
    """
    List_3_cels_for_each_vect = []
    nb_vectors_diff = len(List_bool_vect)
    for v in range(nb_vectors_diff):
        
        List_3_cels = [0, 0, 0] #La liste des 3 cellules qu'on va choisir
        List_class_of_List_3_cel = [] #La liste des classes correspondantes à ces trois cellules
        add = 0
        
        while add != 3:
            #On prend une cellule aléatoirement dans la liste des cellules qui ont le même vecteur booléen 
            r = np.random.randint(0,len(List_index_same_vect_bool[v])) #On la choisi parmis la liste qu'on connait
            indice_cel = List_index_same_vect_bool[v][r] #On prend son indice dans le dataframe
            cel_class = df['classes'][indice_cel] #on prend aussi sa classe
            
            if not(cel_class in List_class_of_List_3_cel): #Si on a pas déja une cellule de cette classe
                add += 1
                List_3_cels[cel_class] = indice_cel #On l'ajoute à la liste
                List_class_of_List_3_cel.append(cel_class) #On ajoute aussi sa classe à la liste
        
        #Une fois qu'on a choisi nos 3 cellules (une pour chaque classe)
        #On peut l'ajouter à la liste globale
        List_3_cels_for_each_vect.append(List_3_cels)
    
    return List_3_cels_for_each_vect
# #Test
# test_init = define_3_cel_for_each_vect(D,long_vect_bool,test_find[0],test_find[1])
    
def choix3cel(index_vect_bool,List_3_cels,List_index_same_vect_bool,df,cel_to_change=None):
    """
    Fonction qui prend en argument 3 cellules ayant le même vecteur booléen
    mais qui sont chacune de classe différentes.
    La fonction renvoie un triplet de cellules dont 2 sont les mêmes que celles
    prises en argument et la troisième fait partie de la même classe que celle
    qui a été supprimé du triplet précédent.
    """
    if cel_to_change is None:
        cel_to_change = np.random.randint(0, 3)
    index_cel_to_change = List_3_cels[cel_to_change]
    class_to_change = df['classes'][index_cel_to_change]
    
    change = False
    while not(change):
        index_new_cel = random.sample(List_index_same_vect_bool[index_vect_bool], 1)[0]
        class_new_cel = df['classes'][index_new_cel]
        change = index_new_cel != index_cel_to_change and class_to_change == class_new_cel
    
    List_3_cels[cel_to_change] = index_new_cel

# #Test
# print(test_init[0])
# choix3cel(0,test_init[0],test_find[0],D)
# print(test_init[0])
    

def calculate_score(gene,vecClass1,vecClass2,df):
    """
    Fonction qui, pour un gène,
    calcule le score de différence entre les vecteur booléen.
    """
    return np.linalg.norm(df.iloc[vecClass1][gene].values - df.iloc[vecClass2][gene].values)
# #Test
# print(test_init[1])
# print(calculate_score('139', test_init[0], test_init[1], D))

def compute_score(gene, test_init, df):
    score = 0
    for i in range(3):
        for j in range(i + 1, 4):
            score += calculate_score(gene, test_init[i], test_init[j], df)
    return score / 6

# Test
# print(compute_score('139', test_init, D))

def compute_score_global(test_init, df):
    score = 0
    for gene in df.keys()[long_vect_bool:-2]:
        score += compute_score(gene, test_init, df)
    return score / len(df.keys()[long_vect_bool:-2])

# Test
# print(compute_score_global(test_init, D))

def maximize_score(df):
    """
    Fonction qui cherche le maximum du score
    """
    vect_find = find_same_vect_bool(df)
    init = define_3_cel_for_each_vect(df,long_vect_bool,vect_find[0],vect_find[1])
    iter = 0
    while iter < 1000:
        iter += 1
        init_copy = init.copy()
        score = compute_score_global(init_copy, df)
        for i in range(4):
            init_copy2 = init_copy.copy()
            choix3cel(i,init_copy2[i],vect_find[0],df)
            new_score = compute_score_global(init_copy2, df)
            if new_score > score:
                init_copy = init_copy2.copy()
                score = new_score
        if new_score > score:
            init = init_copy
        if iter%10 == 0: print(score, iter)
    return init

maximize_score(D)