from __future__ import division
import numpy.matlib
import numpy as np

n_base = 2
p_stay = 0.8
p_leave = (1-p_stay)/(n_base-1)

# priori para a raiz
priori = np.matrix(np.full((1, n_base), 1.0/n_base)).transpose()
# matriz de transição da caracteristica
p_matrix = (p_stay-p_leave)*np.matlib.identity(n_base) + p_leave

class Especie:
    def __init__(self, meu_pai, meu_valor):
        self.filhos = []
        self.pai = meu_pai
        self.valor = meu_valor
        self.L_condicional = None
        if(meu_pai):
            meu_pai.filhos.append(self)
    
    def cria_L_condicional(self):
        if(not(self.filhos)): # self eh uma folha
            self.L_condicional = p_matrix[self.valor, ].transpose()
        elif(not(self.pai)): # self eh a raiz
            for filho in self.filhos:
                filho.cria_L_condicional()
            L_filhos = np.multiply(self.filhos[0].L_condicional, self.filhos[1].L_condicional)
            self.L_condicional = np.multiply(priori, L_filhos)
        else: # self eh no interno
            for filho in self.filhos:
                filho.cria_L_condicional()
            L_filhos = np.multiply(self.filhos[0].L_condicional, self.filhos[1].L_condicional)
            self.L_condicional = numpy.dot(p_matrix, L_filhos)

S0 = Especie(None, None)
S6 = Especie(S0, None)
S1 = Especie(S6, 1)
S2 = Especie(S6, 0)
S8 = Especie(S0, None)
S3 = Especie(S8, 0)
S7 = Especie(S8, None)
S4 = Especie(S7, 0)
S5 = Especie(S7, 1)
S0.cria_L_condicional()
np.sum(S0.L_condicional)

# Achar a arvore de maxima verossimilhança e
# ler The base substitution probabilities
# Bonus: adaptar o que fizemos para o caso em que
# valor pode ser "A", "T", "C" ou "G".

# Estudar numpy, classes e o codigo
# Bonus: pulley principle
