from __future__ import division
import numpy.matlib
import numpy as np
import math


n_base = 4
u = 1  #taxa de substituição de base por tempo
# priori para a raiz
priori = np.matrix(np.full((n_base, 1), 1.0/n_base))


class Especie:
    def __init__(self, meu_pai, meu_valor, minha_priori = None, meu_tempo = 0):
        self.filhos = []
        self.pai = meu_pai
        self.valor = meu_valor
        self.L_condicional = None
        if(meu_pai):
          meu_pai.filhos.append(self)
          self.priori = meu_pai.priori
          self.tempo = meu_tempo
          self.trs = self.cria_transicao()
        else:
          self.priori = minha_priori
          
    def cria_transicao(self):
      P = np.identity(n_base)*(math.exp(-u*self.tempo))
      for i in range(P[0,].size):
        P[i,] += ((1 - math.exp(-u*self.tempo))*self.priori[i,0])
      return P.transpose()
      
    def cria_L_condicional(self):
        if(not(self.filhos)): # self eh uma folha
            self.L_condicional = self.trs[self.valor, ].transpose()
        elif(not(self.pai)): # self eh a raiz
            for filho in self.filhos:
                filho.cria_L_condicional()
            L_filhos = np.multiply(self.filhos[0].L_condicional, self.filhos[1].L_condicional)
            self.L_condicional = np.multiply(self.priori, L_filhos)
            self.L_arvore = np.sum(self.L_condicional)
        else: # self eh no interno
            for filho in self.filhos:
                filho.cria_L_condicional()
            L_filhos = np.multiply(self.filhos[0].L_condicional, self.filhos[1].L_condicional)
            self.L_condicional = np.dot(self.trs, L_filhos)

S0 = Especie(None, None, minha_priori = priori)
S6 = Especie(S0, None, meu_tempo = 6)
S1 = Especie(S6, 1, meu_tempo = 1)
S2 = Especie(S6, 0, meu_tempo = 0.7 )
S8 = Especie(S0, None, meu_tempo = 0.4)
S3 = Especie(S8, 0, meu_tempo = 0.6)
S7 = Especie(S8, None, meu_tempo = 0.9)
S4 = Especie(S7, 0, meu_tempo = 0.4)
S5 = Especie(S7, 1, meu_tempo = 0.3)
S0.cria_L_condicional()
S0.L_arvore

# the base substitution probabilities
# prox reuniao: implementar o p_matrix que depende do tempo.
# pulley principle
