# Ce script développe des methodes de discrimination de readout pour du multi classe

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import copy

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
    ###FIXME###
    ###Ajouter descriptiuon de l'ensemble des paramètres###
    
    def __init__(self, data, optim, long_vect_bool,genes_to_optim = None, children=None):
        """
        Parameters
        ----------
        data : pandas.dataframe
            Les données possèdent la forme suivante :
            - des vecteurs booléens ;
            - puis les gènes de readouts ;
            - classes : attribu détaillant l'appartenance aux différentes classes (entre 0 et nclasse - 1).
        
        optim : str or None
            Si optim vaut 'p' : on applique la recherche en parallèle
            SI optim vaut 's' : on 
        
        genes_to_optim : List[str] or int
            La liste des gènes sur lesquels on calculs les scores
            si un entier n est choisi on prendra les n premiers
        """
        self.data = data
        self.score = 0
        self.optimization = optim
        self.same_vect_bool = None
        self.bool_vect = None
        self.init = None
        self.long_vect_bool = long_vect_bool
        self.nb_classes = len(self.data['classes'].unique())
        if genes_to_optim is None:
            genes_to_optim = [str(i+self.long_vect_bool) for i in range(len(data.keys())-long_vect_bool-1)]
        elif type(genes_to_optim) == int:
            genes_to_optim = [str(i+self.long_vect_bool) for i in range(genes_to_optim)]
        self.genes_to_optim = genes_to_optim
        self.children = children
    
    def __repr__(self) -> str:
        print(f"data = {self.data.head(5)}")
        print(f"score = {self.score}")
        print(f"optimization  = {self.optimization}")
        print(f"same_vect_bool = {self.same_vect_bool}")
        print(f"bool_vect = {self.bool_vect}")
        print(f"init = {self.init}")
        print(f"long_vect_bool = {self.long_vect_bool}")
        print(f"nb_classes = {self.nb_classes}")
        return ""
        
    def plot(self) -> None:
        nb_gene = len(self.genes_to_optim)
        #Partie récupération des données
        values_readouts = []
        for list_cells in self.init:
            readouts_genes = [[] for i in range(nb_gene)]
            for numero_cel in list_cells:
                data_cel = self.data.loc[numero_cel][self.genes_to_optim]
                for numero_gene in range(nb_gene):
                    readouts_genes[numero_gene].append(data_cel[numero_gene])
            values_readouts.append(readouts_genes)
        print(values_readouts)
        
        #Partie dessin avec matplotlib
        classes = [i for i in range(1,self.nb_classes+1)]
        shape = (nb_gene,1)
        for i in range(2,6):
            if nb_gene%i == 0 and nb_gene!=i:
                shape = (nb_gene//i,i)
        
        if shape == (nb_gene,1):
            fig, axs = plt.subplots(shape[0], shape[1], figsize=(shape[1]*6, shape[0]*3))
        else :
            fig, axs = plt.subplots(shape[0], shape[1], figsize=(shape[0]*5, shape[1]*3))
        
        for i,ax in enumerate(axs.flatten()):
            for j,readouts_same_vect_bool in enumerate(values_readouts):
                data_1_gene = readouts_same_vect_bool[i]
                ax.plot(classes,data_1_gene,label = f"vb n°{j}")
            ax.set_xticks([i+1 for i in range(self.nb_classes)])
            ax.set_yticks([0,0.25,0.5,0.75,1])
            ax.set_title(f"Readouts du gène n°{self.genes_to_optim[i]}")
        ax.legend()

        # for readouts_same_vect_bool in values_readouts:
        #     for i,(ax,data_1_gene) in enumerate(zip(axs.flatten(),readouts_same_vect_bool)):
        #         ax.plot(classes,data_1_gene,label=f"v b n°{i}")
        #         ax.set_xticks([i+1 for i in range(self.nb_classes)])
        #         ax.set_yticks([0,0.25,0.5,0.75,1])
        #         ax.set_title(f"Readouts du gène n°{self.genes_to_optim[i]}")
        #         ax.legend()
        
            #ax.text(0.5, 0.5, f"Trace des readouts pour tous les vecteurs booléen différents", ha='center', va='bottom', transform=ax.transAxes, fontsize=12, color='red')
        plt.tight_layout()
        plt.legend()
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
                index_bool_classes = [[] for i in range(self.nb_classes)]
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

    def define_cells_for_each_vect(self):
        """
        Fonction qui, pour chaque vecteur booléen différent,
        prend au hasard trois cellules qui partagent ce même vecteur booléen
        mais qui appartiennent à des classes différentes
        
        Renvoie
        -------
        List_cells_for_each_vect : List[List[int]]
            La liste des indexes dans le dataframe des cellules choisies,
            ###FIXME####
            il y en a 12 au total, 3 cellules pour chaque vecteurs
            et il y a 4 vecteurs différents
        """
        List_cells_for_each_vect = []
        nb_vectors_diff = len(self.bool_vect)
        for v in range(nb_vectors_diff):
            
            List_cells = [0 for i in range(self.nb_classes)] #La liste des nb_classes cellules qu'on va choisir
            i = 0
            
            while i < self.nb_classes:
                #On prend une cellule aléatoirement dans la liste des cellules qui ont le même vecteur booléen 
                indice_cel = random.sample(self.same_vect_bool[v][i], 1)[0] #On prend son indice dans le dataframe
                List_cells[i] = indice_cel
                i = i + 1
                
            #Une fois qu'on a choisi nos 3 cellules (une pour chaque classe)
            #On peut l'ajouter à la liste globale
            List_cells_for_each_vect.append(List_cells)
        
        self.init = List_cells_for_each_vect

    def choixcells(self, index_vect_bool, nb_cell_to_change=1):
        ### TODO ajout argument pour changer plus de 1 cellules
        
        """
        Fonction qui prend en argument 3 cellules ayant le même vecteur booléen
        mais qui sont chacune de classe différentes.
        La fonction renvoie un triplet de cellules dont 2 sont les mêmes que celles
        prises en argument et la troisième fait partie de la même classe que celle
        qui a été supprimé du triplet précédent.
        
        nb_cell_to_change = Nombre de cellules à changer, supérieur à 1, FIXME add max value depend nb cell per bool vect
        """
        
        cells_to_change = np.random.randint(0, self.nb_classes, nb_cell_to_change)
        for cell_to_change in cells_to_change:
            index_cel_to_change = self.init[index_vect_bool][cell_to_change]
            class_to_change = self.class_cell(index_cel_to_change)
            
            change = False
            while not(change):
                index_new_cel = random.sample(self.same_vect_bool[index_vect_bool][class_to_change], 1)[0]
                change = index_new_cel != index_cel_to_change
            
            self.init[index_vect_bool][cell_to_change] = index_new_cel

    def calculate_score(self, gene, vecClass1, vecClass2):
        """
        Fonction qui, pour un gène,
        calcule le score de différence entre les vecteurs booléens.
        """
        return np.linalg.norm(self.data.iloc[vecClass1][gene].values - self.data.iloc[vecClass2][gene].values)

    def compute_score(self, gene):
        score = 0
        cpt = 0
        for i in range(self.nb_classes):
            for j in range(i + 1, len(self.bool_vect)):
                cpt += 1
                score += self.calculate_score(gene, self.init[i], self.init[j])
        return score / cpt

    def compute_score_global(self):
        score = 0
        for gene in self.genes_to_optim:
            score += self.compute_score(gene)
        self.score = score / len(self.genes_to_optim)

    def maximize_score(self,max_iter):
        """
        Fonction qui cherche le maximum du score
        """
        if self.bool_vect is None :
            self.find_same_vect_bool()
        if self.init is None :
            self.define_cells_for_each_vect()
        if self.score is None :
            self.compute_score_global()
        iter = 0
        no_change_iter = 0
        while iter < max_iter and no_change_iter < 100:
            iter += 1
            old_init = self.init.copy()
            old_score = copy.copy(self.score)
            for i in range(len(self.bool_vect)):
                self.choixcells(i, nb_cell_to_change=int(1 + (no_change_iter > 10) + (no_change_iter > 30)))
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
        Liste_discrimination = [Discrimination(self.data,None,self.long_vect_bool,self.genes_to_optim) for i in range(nb_parallele)]
        self.children = Liste_discrimination
        max_score,id_max_score = 0,0
        for i in range(max_iter_global):
            for j in range(len(Liste_discrimination)):
                print(f"Maximisation du n°{j} en cours")
                Liste_discrimination[j].maximize_score(max_iter_tmp)
                if Liste_discrimination[j].score > max_score:
                    max_score = Liste_discrimination[j].score
                    id_max_score = j
            for j in range(len(Liste_discrimination)):
                Liste_discrimination[j].init = copy.deepcopy(Liste_discrimination[id_max_score].init)
                Liste_discrimination[j].score = max_score
            print(f"Le score maximal a été trouvé par le numéro {id_max_score} et vaut {max_score}")
        
        Liste_discrimination[id_max_score].plot()
                    
            
        

PATH = "../DONNEES/toy_datasets/readout_fictifs.csv"

D = pd.read_csv(PATH)

# #Test
discrimination = Discrimination(D, None, 10,['10','11','13','14','17','19'])
discrimination.recherche_parallele(3,2,20)
#discrimination.maximize_score(10)
#discrimination.plot()

