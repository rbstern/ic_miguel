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

class No:
    def __init__(self, indice, vizinhos = None, atributos = None):
        self.indice = indice
        
        if vizinhos is None:
            self.dici_vizinhos = {}
        else:
            self.dici_vizinhos = vizinhos
        
        if atributos is None:
            self.dici_atributos = {}
        else:
            self.dici_atributos = atributos
            
    def adiciona_vizinho(self, indice_viz, peso=1):
        self.dici_vizinhos[indice_viz] = peso
    
    def retorna_atributos(self):
        return self.dici_atributos
    
    def indices_vizinhos(self):
        return list(self.dici_vizinhos.keys())
    
    def retorna_peso(self, indice_viz):
        if indice_viz in self.dici_vizinhos:
            return self.dici_vizinhos[indice_viz]
        else:
            print('No {} e {} não são vizinhos'.format(self.indice, indice_viz))


class Grafo:
    def __init__(self, raiz):
        self.dici_nos = {}
        self.transforma_arv(raiz)
    
    def transforma_arv(self,especie):
        if (not(especie.pai)):  #eh raiz
            for filho in especie.filhos:
                no = No(filho, None, None)
                self.dici_nos [filho] = no
            self.dici_nos[especie.filhos[0]].adiciona_vizinho(especie.filhos[1],
                                            especie.filhos[0].tempo + especie.filhos[1].tempo)
            self.dici_nos[especie.filhos[1]].adiciona_vizinho(especie.filhos[0],
                                            especie.filhos[0].tempo + especie.filhos[1].tempo)
            self.transforma_arv(especie.filhos[0])
            self.transforma_arv(especie.filhos[1])
            
        elif(not(especie.filhos)):   #eh folha
            self.dici_nos[especie].dici_atributos["valor"] = especie.valor
            
        elif(especie.filhos is not None):     #eh interno
            for filho in especie.filhos:
                no = No(filho, None, None)
                self.dici_nos[filho] = no
            self.dici_nos[especie].adiciona_vizinho(especie.filhos[0],especie.filhos[0].tempo)
            self.dici_nos[especie].adiciona_vizinho(especie.filhos[1], especie.filhos[1].tempo)
            
            if (especie.pai.pai is not None):
                self.dici_nos[especie.filhos[0]].adiciona_vizinho(especie, especie.filhos[0].tempo)
                self.dici_nos[especie.filhos[1]].adiciona_vizinho(especie, especie.filhos[1].tempo)
                
            self.transforma_arv(especie.filhos[0])
            self.transforma_arv(especie.filhos[1])
            
    def no(self, indice_no):
        return self.dici_nos[indice_no]
    
    def nos_grafo(self):
        return list(self.dici_nos.values())
                

g = Grafo(S0)
print(g.dici_nos[S6].retorna_peso(S8))
print(g.dici_nos[S6].indices_vizinhos())





# pulley principle: provar que vale sempre
# começar leitura de finding maximum likelihood tree...
# ir ate Finding Optimal Segment Lengths
