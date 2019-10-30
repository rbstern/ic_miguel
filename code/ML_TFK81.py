from __future__ import division
import numpy.matlib
import numpy as np
import math


n_base = 4
u = 1  #taxa de substituiÃ§Ã£o de base por tempo
# priori para a raiz
priori = np.matrix(np.full((n_base, 1), 1.0/n_base))



class Especie:
    def __init__(self, meu_pai, meu_valor, minha_priori = None, meu_tempo = 0):
        self.filhos = []
        self.pai = meu_pai
        self.valor = meu_valor
        self.num_codon = self.valor.size
        self.L_condicional = np.full(self.num_codon, None)
        self.P_trs = np.full(self.num_codon, None)
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
    
    def cria_L_condicional_vetor(self):
      self.L_arvore = np.zeros(self.num_codon)
      for i in range(self.num_codon):
        self.cria_L_condicional(i)
    
    def cria_L_condicional(self,i):
      if(not(self.filhos)):# self eh uma folha
          valor = np.zeros(n_base)
          valor[self.valor[i]] =  1
          self.L_condicional[i] = valor
      elif(not(self.pai)): # self eh a raiz
          for filho in self.filhos:
              filho.cria_L_condicional(i)
              filho.P_trs[i] = np.dot(filho.L_condicional[i],filho.trs).transpose()
          self.L_condicional[i] = np.multiply(self.filhos[0].P_trs[i], self.filhos[1].P_trs[i])
          self.P_trs[i] = np.multiply(self.priori, self.L_condicional[i])
          self.L_arvore[i] = np.sum(self.P_trs[i])
      else: # self eh no interno
          for filho in self.filhos:
              filho.cria_L_condicional(i)
              filho.P_trs[i] = np.dot(filho.L_condicional[i],filho.trs).transpose()
          self.L_condicional[i] = np.multiply(self.filhos[0].P_trs[i], self.filhos[1].P_trs[i])

S0 = Especie(None, np.full(2,None), minha_priori = priori)
dici = {S0 : 1}
S6 = Especie(S0, np.full(2,None), meu_tempo = 6)
S1 = Especie(S6, np.array([1,2]), meu_tempo = 1)
S2 = Especie(S6, np.array([0,2]), meu_tempo = 0.7 )
S8 = Especie(S0, np.full(2,None), meu_tempo = 0.4)
S3 = Especie(S8, np.array([0,3]), meu_tempo = 0.6)
S7 = Especie(S8, np.full(2,None), meu_tempo = 0.9)
S4 = Especie(S7, np.array([0,2]), meu_tempo = 0.4)
S5 = Especie(S7, np.array([1,3]), meu_tempo = 0.3)
S0.cria_L_condicional_vetor()
print(S6.L_condicional)
print(S0.L_condicional)
print(S0.L_arvore)



class No:
    def __init__(self, indice, vizinhos = None, atributos = None, priori = None):
        self.indice = indice
        self.valor = atributos
        self.num_codon = self.valor.size
        self.L_condicional = np.full(self.num_codon, None)
        self.tempo = None
        self.P_trs = np.full(self.num_codon, None)
        self.priori = priori
        self.trs = None
        
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
            print('No {} e {} nÃ£o sÃ£o vizinhos'.format(self.indice, indice_viz))
            
    def cria_transicao_no(self):
      P = np.identity(n_base)*(math.exp(-u*self.tempo))
      for i in range(P[0,].size):
        P[i,] += ((1 - math.exp(-u*self.tempo))*self.priori[i,0])
      return P.transpose()
        


class Transforma_arv:
    def __init__(self, raiz):
        self.dici_nos = {}
        
    
    def transforma_grafo(self,raiz):
        if (not(raiz.pai)):  #eh raiz
            for filho in raiz.filhos:
                no = No(filho, None, filho.valor, filho.priori)  #adicionando valores da arvore no nï¿½
                self.dici_nos [filho] = no
            self.dici_nos[raiz.filhos[0]].adiciona_vizinho(raiz.filhos[1],
                                            raiz.filhos[0].tempo + raiz.filhos[1].tempo)
            self.dici_nos[raiz.filhos[1]].adiciona_vizinho(raiz.filhos[0],
                                            raiz.filhos[0].tempo + raiz.filhos[1].tempo)
            self.transforma_grafo(raiz.filhos[0])
            self.transforma_grafo(raiz.filhos[1])
            
        elif(not(raiz.filhos)):   #eh folha
            self.dici_nos[raiz].valor = raiz.valor
            
        elif(raiz.filhos is not None):     #eh interno
            for filho in raiz.filhos:
                no = No(filho, None, filho.valor, filho.priori)   #adicionando valores da arvore no nï¿½
                self.dici_nos[filho] = no
            self.dici_nos[raiz].adiciona_vizinho(raiz.filhos[0],raiz.filhos[0].tempo)
            self.dici_nos[raiz].adiciona_vizinho(raiz.filhos[1], raiz.filhos[1].tempo)
            
            if (raiz.pai is not None):
                self.dici_nos[raiz.filhos[0]].adiciona_vizinho(raiz, raiz.filhos[0].tempo)
                self.dici_nos[raiz.filhos[1]].adiciona_vizinho(raiz, raiz.filhos[1].tempo)
                
            self.transforma_grafo(raiz.filhos[0])
            self.transforma_grafo(raiz.filhos[1])
            
    def no(self, indice_no):
        return self.dici_nos[indice_no]
    
    def nos_grafo(self):
        return list(self.dici_nos.values())
        
    def L_condicional_g_vetor(self,especie_1, especie_aresta, especie_2 = None):
        for i in range(self.no(especie_1).num_codon):
            self.L_condicional_g(especie_1, especie_aresta, especie_2,i)
    
    def L_condicional_g(self,especie_1, especie_aresta,especie_2,i):  
        if especie_2 == None:  #especie_1 ï¿½ a que vai, especie_2 ï¿½ a que vem, especie_aresta 
            for especie in self.no(especie_1).indices_vizinhos():
                self.L_condicional_g(especie, especie_aresta,especie_1,i)
                self.no(especie).P_trs[i] = np.dot(self.no(especie).L_condicional[i],self.no(especie).trs).transpose()
            viz = np.array(self.no(especie_1).indices_vizinhos())
            viz_temp = viz[viz == especie_aresta]     #aresta de interesse isolada
            viz = viz[viz != especie_aresta]
            if not (viz.size):
              self.no(especie_1).L_condicional[i] = especie_1.L_condicional[i]
            else:
              self.no(especie_1).L_condicional[i] = np.multiply(self.no(viz[0]).P_trs[i],self.no(viz[1]).P_trs[i])
        #completar com a aresta que serï¿½ isolada
        else:
            viz = np.array(self.no(especie_1).indices_vizinhos())
            viz_temp = list((viz[viz == especie_2]))
            self.no(especie_1).tempo = self.no(especie_1).retorna_peso(viz_temp[0])
            self.no(especie_1).trs = self.no(especie_1).cria_transicao_no()
            viz = viz[viz != especie_2]
            if not(viz.size):  #eh folha
                self.no(especie_1).L_condicional[i] = especie_1.L_condicional[i]
            else:   #eh interno
                for especie in viz:
                    self.L_condicional_g(especie, especie_aresta, especie_1,i)
                    self.no(especie).P_trs[i] = np.dot(self.no(especie).L_condicional[i],self.no(especie).trs).transpose()
                viz = viz[viz != especie_2]
                self.no(especie_1).L_condicional[i] = np.multiply(self.no(viz[0]).P_trs[i],self.no(viz[1]).P_trs[i])
                           

g = Transforma_arv(S0)
g.transforma_grafo(S0)
g.L_condicional_g_vetor(S8,S7)
print(g.no(S1).L_condicional)
print(g.no(S7).L_condicional)
print(g.no(S8).L_condicional)
