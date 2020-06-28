# Example of biopython use
import numpy as np
from Bio import Phylo
import ML_TFK81 as fels

# reading files from sanson
sanson_tree = Phylo.read('Sanson.tre', 'nexus')
print(sanson_tree)

# trying to encode the sequence data from the tips and estimating the priori
from Bio import SeqIO
seq_list = []
id_list = []
for seq_record in SeqIO.parse("SansonLeaves.nex", "nexus"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    id_list.append(seq_record.id)
    seq_list.append(str(seq_record.seq))

# making lists of arrays encoding the step above
# encoding 
def dict_encoder(seq_list, id_list):   
    array_dict = {}
    np.random.seed(0)
    del_list = []
    for i in range(len(seq_list)):
        array = np.zeros(len(seq_list[0]))
        for j in range(len(seq_list[0])):
            if seq_list[i][j] == 'T':
                array[j] = 0
            elif seq_list[i][j] == 'C':
                array[j] = 1
            elif seq_list[i][j] == 'A':
                array[j] = 2
            elif seq_list[i][j] == 'G':
                array[j] = 3
            else:
                del_list.append(j)
                
        array_dict[id_list[i]] = array
    return array_dict, del_list

array_dict, del_list = dict_encoder(seq_list, id_list)[0], dict_encoder(seq_list, id_list)[1]
  
# deleting the missing values in all arrays
def priori_maker(array_dict):
    count_T = 0
    count_C = 0
    count_A = 0
    count_G = 0
    for key in array_dict:
        array_dict[key] = np.delete(array_dict[key], del_list)
        count_T += (array_dict[key] == 0).sum()
        count_C += (array_dict[key] == 1).sum()
        count_A += (array_dict[key] == 2).sum()
        count_G += (array_dict[key] == 3).sum() 
    priori_T = count_T/(count_T + count_C + count_A + count_G)
    priori_C = count_C/(count_T + count_C + count_A + count_G)
    priori_A = count_A/(count_T + count_C + count_A + count_G)
    priori_G = count_G/(count_T + count_C + count_A + count_G)
    
    priori = np.array([[priori_T, priori_C, priori_A, priori_G]]).transpose()
    return priori

priori = priori_maker(array_dict)
np.random.uniform(0, 1, 3)
# drawing the tree
Phylo.draw(sanson_tree)
n_pos = array_dict[id_list[0]].size

u = 1
# itering and making the felsenstein tree
def iter_tree(tree, specie_bio, parents_fels):
    if specie_bio == None and parents_fels == None:
        S0 = fels.Especie(None, np.full(n_pos, None), minha_priori =  priori, taxa_subst = u)
        specie_bio1 = tree.clade[0]
        specie_bio2 = tree.clade[1]
        iter_tree(tree, specie_bio1, S0)
        iter_tree(tree, specie_bio2, S0)
        return S0
    
    elif specie_bio.name == None:
        unif = np.random.uniform(3,5)
        S = fels.Especie(parents_fels, np.full(n_pos, None), meu_tempo = unif, taxa_subst = u)
        specie_bio1 = specie_bio.clades[0]
        specie_bio2 = specie_bio.clades[1]
        iter_tree(tree, specie_bio1, S)
        iter_tree(tree, specie_bio2, S)
        
    elif specie_bio.name != None:
        unif = np.random.uniform(0,2)
        fels.Especie(parents_fels, array_dict[specie_bio.name], meu_tempo = unif, taxa_subst = u)

# testing now
np.random.seed(0)
S0 = iter_tree(sanson_tree, None, None)

# Now that we have our felsenstein tree estabilished, lets compute the likelihood and 
# the optimum lengths more further
filho = S0.filhos[0].filhos[0].filhos[0].filhos[0].valor
# creating 
S0.cria_L_condicional_vetor()
print(S0.L_arvore)
print(np.prod(S0.L_arvore))
print(np.sum(np.log(S0.L_arvore)))
veros_1 = S0.L_arvore

S0.modifica_L_condicional()
print(np.sum(np.log(S0.L_arvore)))

# making non rooted tree
g = fels.Grafo()
g.transforma_grafo(S0)
g.maxim_L_vetor((10**(-4)))
g.muda_tempo_arv()
S0.modifica_L_condicional()
print(np.prod(S0.L_arvore))
print(np.sum(np.log(S0.L_arvore)))
veros_2 = S0.L_arvore
        
