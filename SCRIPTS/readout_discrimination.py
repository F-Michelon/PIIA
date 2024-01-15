# Ce script développe des methodes de discrimination de readout pour du multi classe

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

# TODO décrire la tête des données formattage
# paramétriser l'ensemble des données numérique
'''
Les données possèdent la forme suivante :
    - des vecteurs booléens de 120 gènes ;
    - puis 20 gènes de readouts ;
    - vect_bool_alike : attribu décrivant l'appartenance aux différents vecteurs booléens (pas forcémment utile pour l'instant) ;
    - classe : attribu détaillant l'appartenance aux différentes classes (entre 0 et nclasse - 1).
'''

class Discrimination:
    def __init__(self, data, optim, long_vect_bool):
        self.data = data
        self.score = 0
        self.optimization = optim
        self.same_vect_bool = None
        self.bool_vect = None
        self.init = None
        self.long_vect_bool = long_vect_bool
        self.nb_classe = len(self.data['classes'].unique())
        
    def plot(self) -> None:
        nb_gene = 5
        #Partie récupération des données
        values_readouts = []
        for list_3_cels in self.init:
            readouts_genes = [[] for i in range(nb_gene)]
            for numero_cel in list_3_cels:
                data_cel = self.data.loc[numero_cel][self.long_vect_bool:self.long_vect_bool+5]
                for numero_gene in range(nb_gene):
                    readouts_genes[numero_gene].append(data_cel[numero_gene])
            values_readouts.append(readouts_genes)
        print(values_readouts)
        
        #Partie dessin avec matplotlib
        classes = [i for i in range(1,self.nb_classe+1)]
        fig, axs = plt.subplots(nb_gene, 1, figsize=(5, 10))
    
        for readouts_same_vect_bool in values_readouts:
            for i,(ax,data_1_gene) in enumerate(zip(axs.flatten(),readouts_same_vect_bool)):
                ax.plot(classes,data_1_gene)
                ax.set_xticks([1,2,3])
                ax.set_yticks([0,0.25,0.5,0.75,1])
                ax.set_title(f"Readouts du gène n°{i}")
            #ax.text(0.5, 0.5, f"Trace des readouts pour tous les vecteurs booléen différents", ha='center', va='bottom', transform=ax.transAxes, fontsize=12, color='red')
        plt.tight_layout()
        plt.show()


    def find_same_vect_bool(self):
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
        N = len(self.data)
        List_bool_vect = []
        List_index_same_vect_bool = []
        for i in range(N):
            vect_bool = [[self.data.iloc[i][k] for k in range(self.long_vect_bool)]] #Le vecteur booléen de la ligne i du tableau
            if not(vect_bool in List_bool_vect): #S'il existe pas déja
                List_bool_vect.append(vect_bool) #On l'ajoute à la liste
                index_bool_classes = [[] for i in range(self.nb_classe)]
                index_bool_classes[self.class_cell(i)].append(i)
                List_index_same_vect_bool.append(index_bool_classes) #Et on ajoute son indice à la liste des numéros de chaque vecteurs boléens
            
            else : #S'il est dans la liste
                for j in range(len(List_bool_vect)): #On cherche quel vecteur c'est
                    if List_bool_vect[j] == vect_bool: #Si c'est le vecteur n° j de la liste
                        List_index_same_vect_bool[j][self.class_cell(i)].append(i) #On ajoute son indice à la liste
        
        self.same_vect_bool = List_index_same_vect_bool
        self.bool_vect = List_bool_vect

    def class_cell(self, index_cell):
        # Renvoie le numéro de la classe de la cellule
        return self.data['classes'][index_cell]

    def define_3_cel_for_each_vect(self):
        """
        Fonction qui, pour chaque vecteur booléen différent,
        prend au hasard trois cellules qui partagent ce même vecteur booléen
        mais qui appartiennent à des classes différentes
        
        Renvoie
        -------
        List_3_cels_for_each_vect : List[List[int]]
            La liste des indexes dans le dataframe des cellules choisies,
            il y en a 12 au total, 3 cellules pour chaque vecteurs
            et il y a 4 vecteurs différents
        """
        List_3_cels_for_each_vect = []
        nb_vectors_diff = len(self.bool_vect)
        for v in range(nb_vectors_diff):
            
            List_3_cels = [0, 0, 0] #La liste des 3 cellules qu'on va choisir
            i = 0
            
            while i < self.nb_classe:
                #On prend une cellule aléatoirement dans la liste des cellules qui ont le même vecteur booléen 
                indice_cel = random.sample(self.same_vect_bool[v][i], 1)[0] #On prend son indice dans le dataframe
                List_3_cels[i] = indice_cel
                i = i + 1
                
            #Une fois qu'on a choisi nos 3 cellules (une pour chaque classe)
            #On peut l'ajouter à la liste globale
            List_3_cels_for_each_vect.append(List_3_cels)
        
        self.init = List_3_cels_for_each_vect

    def choix3cel(self, index_vect_bool, nb_cell_to_change=1):
        ### TODO ajout argument pour changer plus de 1 cellules
        
        """
        Fonction qui prend en argument 3 cellules ayant le même vecteur booléen
        mais qui sont chacune de classe différentes.
        La fonction renvoie un triplet de cellules dont 2 sont les mêmes que celles
        prises en argument et la troisième fait partie de la même classe que celle
        qui a été supprimé du triplet précédent.
        
        nb_cell_to_change = Nombre de cellules à changer, supérieur à 1, FIXME add max value depend nb cell per bool vect
        """
        
        cells_to_change = np.random.randint(0, self.nb_classe, nb_cell_to_change)
        for cell_to_change in cells_to_change:
            index_cel_to_change = self.init[index_vect_bool][cell_to_change]
            class_to_change = self.class_cell(index_cel_to_change)
            
            change = False
            while not(change):
                index_new_cel = random.sample(self.same_vect_bool[index_vect_bool][class_to_change], 1)[0]
                change = index_new_cel != index_cel_to_change
            
            self.init[index_vect_bool][cell_to_change] = index_new_cel

    def calculate_score(self, id_gene, vecClass1, vecClass2):
        """
        Fonction qui, pour un gène,
        calcule le score de différence entre les vecteurs booléens.
        """
        return np.linalg.norm(self.data.iloc[vecClass1][id_gene].values - self.data.iloc[vecClass2][id_gene].values)

    def compute_score(self, gene):
        score = 0
        for i in range(self.nb_classe):
            for j in range(i + 1, len(self.bool_vect)):
                score += self.calculate_score(gene, self.init[i], self.init[j])
        return score / 6

    def compute_score_global(self):
        score = 0
        for gene in self.data.keys()[self.long_vect_bool:-2]:
            score += self.compute_score(gene)
        self.score = score / len(self.data.keys()[self.long_vect_bool:-2])

    def maximize_score(self,max_iter):
        """
        Fonction qui cherche le maximum du score
        """
        self.find_same_vect_bool()
        self.define_3_cel_for_each_vect()
        iter = 0
        no_change_iter = 0
        self.compute_score_global()
        while iter < max_iter and no_change_iter < 100:
            iter += 1
            old_init = self.init.copy()
            old_score = self.score.copy()
            for i in range(len(self.bool_vect)):
                self.choix3cel(i, nb_cell_to_change=int(1 + (no_change_iter > 10) + (no_change_iter > 30)))
                self.compute_score_global()
                if old_score >= self.score:
                    self.init = old_init
                    self.score = old_score
                    no_change = True
                else:
                    no_change_iter = 0
                    no_change = False
            no_change_iter += 1 * (no_change)
            if iter%10 == 0:
                print(self.score, iter, no_change_iter)
                #self.plot()
            
    def recherche_parallele(self,nb_parallele,max_iter_global,max_iter_tmp):
        """
        Fonction qui fait plusieurs recherche en partant de points différents et
        qui regarde qui avance le plus vite pour les faire tous progresser
        """
        Liste_discrimination = [Discrimination(self.data,None,120) for i in range(nb_parallele)]
        
        max_score,id_max_score = 0,0
        for i in range(max_iter_global):
            for j in range(len(Liste_discrimination)):
                print(f"Maximisation du n°{j} en cours")
                D = Liste_discrimination[j]
                D.maximize_score(max_iter_tmp)
                if D.score > max_score:
                    max_score = D.score
                    id_max_score = j
            for D in Liste_discrimination:
                D.init = Liste_discrimination[id_max_score].init
            print(f"Le score maximal a été trouvé par le numéro {id_max_score} et vaut {max_score}")
        
        Liste_discrimination[id_max_score].plot()
                    
            
        

PATH = "../DONNEES/toy_datasets/readout_fictifs_D_3.csv"

D = pd.read_csv(PATH)

# #Test
discrimination = Discrimination(D, None, 120)
discrimination.recherche_parallele(10,50,5)
#discrimination.maximize_score(10)
#discrimination.plot()

