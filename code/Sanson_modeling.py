# Example of biopython use
import numpy as np
from Bio import Phylo
import ML_TFK81 as fels

n_pos = 2236
# reading files from sanson
sanson_tree = Phylo.read('Sanson.tre', 'nexus')
print(sanson_tree)

# trying to encode the sequence data from the tips and estimating the priori
count_T = 0
count_C = 0
count_A = 0
count_G = 0
count_empty = 0
from Bio import SeqIO
seq_list = []
id_list = []
for seq_record in SeqIO.parse("SansonLeaves.nex", "nexus"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    id_list.append(seq_record.id)
    seq_list.append(str(seq_record.seq))
    count_T  += seq_record.seq.count("T")
    count_C  += seq_record.seq.count("C")
    count_A  += seq_record.seq.count("A")

priori_T = count_T/(count_T + count_C + count_A + count_G)
priori_C = count_C/(count_T + count_C + count_A + count_G)
priori_A = count_A/(count_T + count_C + count_A + count_G)
priori_G = count_G/(count_T + count_C + count_A + count_G)

priori = np.array([[priori_T, priori_C, priori_A, priori_G]]).transpose()
# making lists of arrays encoding the step above
# encoding 
array_dict = {}
np.random.seed(0)
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
            u = np.random.uniform(0,1)
            cumulative = np.cumsum(priori)
            if 0 < u <=  cumulative[0]:
                array[j] = 0
            elif cumulative[0] < u <= cumulative[1]:
                array[j] = 1
            elif cumulative[1] < u <= cumulative[2]:
                array[j] = 2
            elif cumulative[2] < u <= cumulative[3]:
                array[j] = 3
    array_dict[id_list[i]] = array


np.random.uniform(0, 1, 3)
# drawing the tree
Phylo.draw(sanson_tree)

# itering and making the felsenstein tree
def iter_tree(tree, specie_bio, parents_fels):
    if specie_bio == None and parents_fels == None:
        S0 = fels.Especie(None, np.full(n_pos, None), minha_priori =  priori)
        specie_bio1 = tree.clade[0]
        specie_bio2 = tree.clade[1]
        iter_tree(tree, specie_bio1, S0)
        iter_tree(tree, specie_bio2, S0)
        return S0
    
    elif specie_bio.name == None:
        unif = np.random.uniform(0,16)
        S = fels.Especie(parents_fels, np.full(n_pos, None), meu_tempo = unif)
        specie_bio1 = specie_bio.clades[0]
        specie_bio2 = specie_bio.clades[1]
        iter_tree(tree, specie_bio1, S)
        iter_tree(tree, specie_bio1, S)
        
    elif specie_bio.name != None:
        unif = np.random.uniform(0,8)
        fels.Especie(parents_fels, array_dict[specie_bio.name], meu_tempo = unif)

# testing now
S0 = iter_tree(sanson_tree, None, None)



        
