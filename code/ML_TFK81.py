import numpy as np
import math



np.random.seed(0)
n_base = 4
# u = 1  
u = np.random.uniform(0.24e-04, 0.42e-04)
#taxa de substituição de base por tempo
# priori para a raiz
priori = np.full((n_base, 1), 1/n_base)
# priori = np.array([[0.123, 0.210, 0.3, 0.367]]).transpose()


class Especie:
    def __init__(self, meu_pai, meu_valor, minha_priori = None, meu_tempo = 0):
        self.filhos = []
        self.pai = meu_pai
        self.valor = meu_valor
        self.num_codon = self.valor.size
        self.prim_calc = True
        self.recalc = np.full(self.num_codon, False)
        self.change_recalc = False
        self.mudou_tempo = False
        self.L_condicional = None
        self.P_trs = None
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
      self.L_arvore = np.zeros((1, self.num_codon))
      self.cria_L_condicional()
    
    def cria_L_condicional(self):
        if(not(self.filhos)):# self eh uma folha
            if (self.prim_calc == True):
                valor = np.zeros((n_base, self.num_codon))
                for i in range(self.num_codon):
                    valor[int(self.valor[i]), i] =  1
                self.L_condicional = valor
                self.P_trs = np.dot(self.trs, self.L_condicional)
                self.prim_calc = False
                self.L_condicional_copy,self.P_trs_copy = self.L_condicional.copy(),self.P_trs.copy()
            elif np.count_nonzero(self.recalc) > 0 and self.prim_calc == False:
                A = list(np.where(self.recalc)[0])
                self.L_condicional_copy[:, A] = np.zeros((n_base, len(A)))
                self.L_condicional_copy[int(self.valor[A]), A] = 1
                self.P_trs_copy[:, A] = np.dot(self.trs, 
                               self.L_condicional_copy[:, A])
                if self.change_recalc == True:
                    self.recalc = np.full(self.num_codon, False)
                
        elif(not(self.pai)): # self eh a raiz
            if (self.prim_calc == True):
                for filho in self.filhos:
                    filho.cria_L_condicional()
                self.L_condicional = np.multiply(self.filhos[0].P_trs, self.filhos[1].P_trs)
                self.L_condicional_copy = self.L_condicional.copy()
                self.L_arvore = np.dot(self.priori.transpose(), self.L_condicional)
                self.L_arvore_copy = self.L_arvore
                self.prim_calc = False
            elif np.count_nonzero(self.recalc) > 0 and self.prim_calc == False:
                A = list(np.where(self.recalc)[0])
                for filho in self.filhos:
                    filho.cria_L_condicional()
                self.L_condicional_copy[:, A] = np.multiply(self.filhos[0].P_trs_copy[:, A], 
                                  self.filhos[1].P_trs_copy[:, A])
                self.L_avore_copy = self.L_arvore.copy()
                self.L_arvore_copy[0, A] = np.dot(self.priori.transpose(), 
                                  self.L_condicional_copy[:, A])
                self.recalc = np.full(self.num_codon, False)
                
        else: # self eh no interno
            if (self.prim_calc == True):
                for filho in self.filhos:
                    filho.cria_L_condicional()
                self.L_condicional = np.multiply(self.filhos[0].P_trs, self.filhos[1].P_trs)
                self.P_trs = np.dot(self.trs, self.L_condicional)
                self.prim_calc = False   
                self.L_condicional_copy,self.P_trs_copy = self.L_condicional.copy(),self.P_trs.copy()
            elif np.count_nonzero(self.recalc) > 0 and self.prim_calc == False:
                A = list(np.where(self.recalc)[0])
                for filho in self.filhos:
                    if self.change_recalc == True:
                        filho.change_recalc = True
                    filho.cria_L_condicional()
                self.L_condicional_copy[:, A] = np.multiply(self.filhos[0].P_trs_copy[:, A],
                                  self.filhos[1].P_trs_copy[:, A])
                self.P_trs_copy[:, A] = np.dot(self.trs, self.L_condicional_copy[:, A])
                if self.change_recalc == True:
                    self.recalc = np.full(self.num_codon, False)
                    
    def prob_condicional(self, val_codon, pos_codon):
        valor_original = self.valor[pos_codon]
        pai_folha = self.pai
        self.recalc[pos_codon] = True
        self.change_recalc = False
        while pai_folha.pai != None:
            no_recalc = pai_folha
            no_recalc.change_recalc = False
            no_recalc.recalc[pos_codon] = True
            pai_folha = pai_folha.pai
        raiz = pai_folha
        self.valor[pos_codon] =  val_codon
        if raiz.prim_calc == True:
            raiz.cria_L_condicional_vetor()
            P_conj = raiz.L_arvore[0, pos_codon]
            print(P_conj)
        else:
            raiz.recalc[pos_codon] = True
            raiz.cria_L_condicional_vetor()
            P_conj = raiz.L_arvore_copy[0, pos_codon]
            print(P_conj)
        P_soma = 0
    # ja temos todas os nos com exceção do que tem desconhecido calculado
    # e não precisamos recalcular para diferentes valores de S1
    #vetor de combinações
    #combin = np.tile(h,folha.num_codon).reshape(folha.num_codon,n_base), combinatorias extensas
        for i in range(n_base):
            self.valor[pos_codon] = i    #para mais de um codon, complica-se a escolha
            if i == (n_base - 1):
                no_recalc.change_recalc = True  
                no_recalc.cria_L_condicional()
            else:
                no_recalc.cria_L_condicional()
            raiz.L_condicional_copy[:, pos_codon] = np.multiply(
                    raiz.filhos[0].P_trs_copy[:, pos_codon], 
                                         raiz.filhos[1].P_trs_copy[:, pos_codon])
            P_soma += np.dot(raiz.priori.transpose(), raiz.L_condicional_copy[:, pos_codon])
        P_cond = P_conj/P_soma
    # voltando para o valor original do começo, antes de se calcular todas as probabilidades
        self.valor[pos_codon] = valor_original
        return(P_cond)
                
    def modifica_L_condicional(self): # Ja calculou L_condicional antes
        if self.prim_calc == False and self.mudou_tempo == True:
            if(not(self.pai)):
                for filho in self.filhos:
                    filho.modifica_L_condicional()
                self.L_condicional = np.multiply(self.filhos[0].P_trs, self.filhos[1].P_trs)
                self.L_arvore = np.dot(self.priori.transpose(), self.L_condicional)
                print(self.L_arvore)
            else:
                for filho in self.filhos:
                    if (not(filho.filhos)):
                        filho.P_trs = np.dot(filho.trs, filho.L_condicional)
                    else:
                        filho.modifica_L_condicional
                        filho.P_trs = np.dot(filho.trs, filho.L_condicional)
                self.L_condicional = np.multiply(self.filhos[0].P_trs, self.filhos[1].P_trs)
                self.P_trs = np.dot(self.trs, self.L_condicional)
        else:
            if(not(self.pai)):
                self.L_arvore = np.zeros(self.num_codon)
                self.cria_L_condicional_vetor()
            else:
                self.cria_L_condicional
        
                


S0 = Especie(None, np.full(8,None), minha_priori = priori)
S6 = Especie(S0, np.full(8,None), meu_tempo = 10)
S1 = Especie(S6, np.array([1, 2, 0, 3, 1, 0, 2, 0]), meu_tempo = 2)
S2 = Especie(S6, np.array([0, 3, 1, 1, 2, 1, 0, 0]), meu_tempo = 4)
S8 = Especie(S0, np.full(8,None), meu_tempo = 15)
S3 = Especie(S8, np.array([0, 1, 1, 2, 2, 3, 3, 0]), meu_tempo = 8)
S7 = Especie(S8, np.full(8,None), meu_tempo = 14)
S4 = Especie(S7, np.array([1, 2, 3, 0, 1, 1, 0, 2]), meu_tempo = 4)
S5 = Especie(S7, np.array([2, 3, 0, 1, 0, 2, 2, 1]), meu_tempo = 5)
S0.cria_L_condicional_vetor()
print(np.sum(np.log(S0.L_arvore)))
print(S1.L_condicional)
print(S1.P_trs)
print(S2.P_trs)
print(S0.P_trs)
print(S6.L_condicional)
print(S6.P_trs)
print(S0.L_condicional)
print(S0.L_arvore)
print(S8.L_condicional)

# Após modificação
S0.modifica_L_condicional()
print(np.sum(np.log(S0.L_arvore)))

 #testando o caso condicional com S1 = 'None' para apenas um codon
#S0 = Especie(S0, np.full(1,None), minha_priori = priori)
#S6 = Especie(S0, np.full(1,None), meu_tempo = 6)
#S1 = Especie(S6, np.array([None]), meu_tempo = 1)
#S2 = Especie(S6, np.array([0]), meu_tempo = 0.7 )
#S8 = Especie(S0, np.full(1,None), meu_tempo = 0.4)
#S3 = Especie(S8, np.array([0]), meu_tempo = 0.6)
#S7 = Especie(S8, np.full(1,None), meu_tempo = 0.9)
#S4 = Especie(S7, np.array([0]), meu_tempo = 0.4)
#S5 = Especie(S7, np.array([3]), meu_tempo = 0.3)



# =============================================================================
# S0 = Especie(None, np.full(2,None), minha_priori = priori)
# S6 = Especie(S0, np.full(2,None), meu_tempo = 15)
# S1 = Especie(S6, np.array([1,2]), meu_tempo = 12)
# S2 = Especie(S6, np.array([0,2]), meu_tempo = 10)
# S8 = Especie(S0, np.full(2,None), meu_tempo = 10)
# S3 = Especie(S8, np.array([0,3]), meu_tempo = 5)
# S7 = Especie(S8, np.full(2,None), meu_tempo = 20)
# S4 = Especie(S7, np.array([None,2]), meu_tempo = 10)
# S5 = Especie(S7, np.array([1,3]), meu_tempo = 3)
# print(S4.prob_condicional(0,0))
# print(S4.prob_condicional(1,0))
# print(S4.prob_condicional(2,0))
# print(S4.prob_condicional(3,0))
# print(S4.prob_condicional(0,0) + S4.prob_condicional(1,0) +
#       S4.prob_condicional(2,0) + S4.prob_condicional(3,0))
# =============================================================================




class No:
    def __init__(self, indice, vizinhos = None, atributos = None, priori = None):
        self.indice = indice
        self.valor = atributos
        self.num_codon = self.valor.size
        self.L_condicional = None
        self.tempo = None
        self.P_trs = None
        self.priori = priori
        self.trs = None
        
        if vizinhos is None:
            self.dici_vizinhos = {}
        else:
            self.dici_vizinhos = vizinhos
            
    def adiciona_vizinho(self, indice_viz, peso=1):
        self.dici_vizinhos[indice_viz] = peso
    
    def indices_vizinhos(self):
        return list(self.dici_vizinhos.keys())
    
    def retorna_peso(self, indice_viz):
        if indice_viz in self.dici_vizinhos:
            return self.dici_vizinhos[indice_viz]
        else:
            print('No {} e {} não são vizinhos'.format(self.indice, indice_viz))
            
    def cria_transicao_no(self):
      P = np.identity(n_base)*(math.exp(-u*self.tempo))
      for i in range(P[0,].size):
        P[i,] += ((1 - math.exp(-u*self.tempo))*self.priori[i,0])
      return P.transpose()
        


class Grafo:
    def __init__(self):
        self.dici_nos = {}
        self.num_codon = None
    def transforma_grafo(self,raiz):
        if (not(raiz.pai)):  #eh raiz
            for filho in raiz.filhos:
                no = No(filho, None, filho.valor, filho.priori)  #adicionando valores da arvore no n�
                self.dici_nos [filho] = no
            self.dici_nos[raiz.filhos[0]].adiciona_vizinho(raiz.filhos[1],
                                            raiz.filhos[0].tempo + raiz.filhos[1].tempo)
            self.dici_nos[raiz.filhos[1]].adiciona_vizinho(raiz.filhos[0],
                                            raiz.filhos[0].tempo + raiz.filhos[1].tempo)
            self.transforma_grafo(raiz.filhos[0])
            self.transforma_grafo(raiz.filhos[1])
            
        elif(not(raiz.filhos)):   #eh folha
            self.dici_nos[raiz].valor = raiz.valor
            if self.num_codon == None:
                self.num_codon = self.dici_nos[raiz].num_codon
            
        elif(raiz.filhos is not None):     #eh interno
            for filho in raiz.filhos:
                no = No(filho, None, filho.valor, filho.priori)   #adicionando valores da arvore no n�
                self.dici_nos[filho] = no
            self.dici_nos[raiz].adiciona_vizinho(raiz.filhos[0],raiz.filhos[0].tempo)
            self.dici_nos[raiz].adiciona_vizinho(raiz.filhos[1], raiz.filhos[1].tempo)
            self.dici_nos[raiz.filhos[0]].adiciona_vizinho(raiz, raiz.filhos[0].tempo)
            self.dici_nos[raiz.filhos[1]].adiciona_vizinho(raiz, raiz.filhos[1].tempo)
                
            self.transforma_grafo(raiz.filhos[0])
            self.transforma_grafo(raiz.filhos[1])
            
    def no(self, indice_no):
        return self.dici_nos[indice_no]
    
    def nos_grafo(self):
        return list(self.dici_nos.keys())
    
    def get_num_aresta(self):
      num_nos = len(self.nos_grafo())
      num_aresta = num_nos -  1
      return(num_aresta)
    
    def muda_peso(self,indice_1,indice_2,peso_novo):
      self.no(indice_1).dici_vizinhos[indice_2] = peso_novo
      self.no(indice_2).dici_vizinhos[indice_1] = peso_novo
      
    
    def maxim_L_vetor(self,e, it = 0):
      erros = np.zeros(self.get_num_aresta())
      lista_nos = self.nos_grafo()
      lista_repetidos = []
      contador = 0
      # inicializando o laço
      while True:
          lista_repetidos = []
          contador = 0
          for no in lista_nos:
            viz = self.no(no).indices_vizinhos()
            tam = len(viz)
            if lista_repetidos == []:
              for i in range(tam):
                erros[contador] = self.maxim_L(no,viz[i])
                contador += 1
            else:
              for i in range(tam):
                if viz[i] not in lista_repetidos:
                  erros [contador] = self.maxim_L(no,viz[i])
                  contador += 1
            lista_repetidos.append(no)
          it += 1
          if np.max(erros) <= e:
              return([erros, it])
    def maxim_L_vetor2(self, e):
        erros = np.zeros(self.get_num_aresta())
        lista_nos = self.nos_grafo()
        lista_repetidos = []
        contador = 0
        for no in lista_nos:
            viz = self.no(no).indices_vizinhos()
            tam = len(viz)
            if lista_repetidos == []:
                for i in range(tam):
                    erros[contador] = self.maxim_L(no, viz[i])
                    while erros[contador] > e:
                        erros[contador] = self.maxim_L(no, viz[i])
                    contador += 1
            else:
                for i in range(tam):
                    if viz[i] not in lista_repetidos:
                        erros[contador] = self.maxim_L(no, viz[i])
                        while erros[contador] > e:
                            erros[contador] = self.maxim_L(no, viz[i])
                        contador += 1
            lista_repetidos.append(no)
            
                            
      
    def maxim_L(self,especie_1,especie_aresta):
      v = self.no(especie_1).retorna_peso(especie_aresta)
      q = math.exp(-v)
      p = 1 - q
      K = self.num_codon
      self.L_condicional_g_vetor(especie_1,especie_aresta)
      A = np.dot(self.no(especie_1).priori.transpose(), np.multiply(self.no(especie_1).L_condicional,
                             self.no(especie_aresta).L_condicional))
      B = np.multiply(np.dot(self.no(especie_1).priori.transpose(), self.no(especie_1).L_condicional),
      np.dot(self.no(especie_aresta).priori.transpose(), self.no(especie_aresta).L_condicional))
      p_atual = (1/K)*(np.sum((p*B)/((A*q) + (B*p))))
      v_atual = -np.log(1 - p_atual)
      self.muda_peso(especie_1,especie_aresta,v_atual)
      erro = abs(p_atual - p)
      return(erro)
      
    
    def L_condicional_g_vetor(self,especie_1, especie_aresta, especie_2 = None):
        self.L_condicional_g(especie_1, especie_aresta, especie_2)
    
    def L_condicional_g(self,especie_1, especie_aresta, especie_2):  
        if especie_2 == None:  #especie_1 � a que vai, especie_2 � a que vem, especie_aresta 
            for especie in self.no(especie_1).indices_vizinhos():
                self.L_condicional_g(especie, especie_aresta,especie_1)
                if especie != especie_aresta:
                    self.no(especie).P_trs = np.dot(self.no(especie).trs, 
                            self.no(especie).L_condicional)
            viz = np.array(self.no(especie_1).indices_vizinhos())
            viz_temp = viz[viz == especie_aresta]     #aresta de interesse isolada
            viz = viz[viz != especie_aresta]
            if not (viz.size):
              self.no(especie_1).L_condicional = especie_1.L_condicional
            else:
              self.no(especie_1).L_condicional = np.multiply(self.no(viz[0]).P_trs,
                      self.no(viz[1]).P_trs)
        else:
            viz = np.array(self.no(especie_1).indices_vizinhos())
            viz_temp = list((viz[viz == especie_2]))
            self.no(especie_1).tempo = self.no(especie_1).retorna_peso(viz_temp[0])
            self.no(especie_1).trs = self.no(especie_1).cria_transicao_no()
            viz = viz[viz != especie_2]
            if not(viz.size):  #eh folha
                self.no(especie_1).L_condicional = especie_1.L_condicional
            else:   #eh interno
                for especie in viz:
                    self.L_condicional_g(especie, especie_aresta, especie_1)
                    self.no(especie).P_trs = np.dot(self.no(especie).trs, 
                            self.no(especie).L_condicional)
                viz = viz[viz != especie_2]
                self.no(especie_1).L_condicional = np.multiply(self.no(viz[0]).P_trs,
                        self.no(viz[1]).P_trs)
        
    def muda_tempo_arv(self):
        lista_nos = self.nos_grafo()
        i = 0
        while i < len(lista_nos):
            if (lista_nos[i]).pai.pai == None:
                if (lista_nos[i]) == (lista_nos[i]).pai.filhos[0]:
                    lista_nos[i].tempo = (self.no(lista_nos[i]).retorna_peso(
                            lista_nos[i].pai.filhos[1]))*1/2
                    lista_nos[i].pai.mudou_tempo = True
                    lista_nos[i].mudou_tempo = True
                    lista_nos[i].trs = lista_nos[i].cria_transicao()
                    i += 1
                else:
                    lista_nos[i].tempo = (self.no(lista_nos[i]).retorna_peso(
                            lista_nos[i].pai.filhos[0]))*1/2
                    lista_nos[i].mudou_tempo = True
                    lista_nos[i].trs = lista_nos[i].cria_transicao()
                    i += 1
            
            else:
                lista_nos[i].tempo = self.no(lista_nos[i]).retorna_peso(lista_nos[i].pai)
                lista_nos[i].mudou_tempo = True
                lista_nos[i].trs = lista_nos[i].cria_transicao()
                i += 1
            


g = Grafo()
g.transforma_grafo(S0)
g.maxim_L_vetor2((10**(-4)))
g.muda_tempo_arv()
# g.no(S5).retorna_peso(S7)

### posteriormente criar classe do tipo "arvore"
### predição só use valores temporarios (val.temp)
### resolver problema na predição (slice e recalc na raiz) e testar para caso com 2 e 3 codons
### foco em apresentar
# =============================================================================
# a = np.array([1,2,3,4,5])
# z = a.copy()
# z[0] = 12
# A = list(np.where(a < 3)[0])
# valor = np.zeros((4, 5))
# for i in range(5):
#     valor[(0,i)] = 1
# valor[:,[0, 2]] = np.zeros((4, 2))
# =============================================================================


