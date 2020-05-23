from __future__ import division
import numpy.matlib
import numpy as np
import math




n_base = 4
u = 1  #taxa de substituição de base por tempo
# priori para a raiz
#priori = np.matrix(np.full((n_base, 1), 1.0/n_base))
priori = np.array([[0.123, 0.210, 0.3, 0.367]]).transpose()


class Especie:
    def __init__(self, meu_pai, meu_valor, minha_priori = None, meu_tempo = 0):
        self.filhos = []
        self.pai = meu_pai
        self.valor = meu_valor
        self.num_codon = self.valor.size
        self.prim_calc = True
        self.recalc = np.full(self.num_codon, False)
        self.L_condicional = np.full((self.num_codon,n_base), None)
        self.P_trs = np.full((self.num_codon,n_base), None)
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
      self.cria_L_condicional()
    
    def cria_L_condicional(self):
        if(not(self.filhos)):# self eh uma folha
            if np.count_nonzero(self.recalc) > 0 or np.all(self.L_condicional) == None:
                valor = np.zeros((self.num_codon,n_base))
                for i in range(self.num_codon):
                    if list(self.recalc)[i] == True or self.prim_calc == True:
                        valor[(i,self.valor[i])] =  1
                self.L_condicional = valor
                self.prim_calc = np.full(self.num_codon, False)
        elif(not(self.pai)): # self eh a raiz
            for filho in self.filhos:
                filho.cria_L_condicional()
            self.L_condicional = np.multiply(self.filhos[0].P_trs, self.filhos[1].P_trs)
            self.P_trs = np.multiply(self.priori.transpose(), self.L_condicional)
            self.L_arvore = np.dot(self.L_condicional, self.priori)
        else: # self eh no interno
            for filho in self.filhos:
                filho.cria_L_condicional()
                filho.P_trs = np.dot(filho.L_condicional,filho.trs)
            if np.all(self.L_condicional) == None and self.prim_calc == True:
                self.L_condicional = np.multiply(self.filhos[0].P_trs, self.filhos[1].P_trs)
                self.P_trs = np.dot(self.L_condicional, self.trs)
                self.prim_calc = False    
            elif np.count_nonzero(self.recalc) > 0 and self.prim_calc == False:
                A = np.where(self.recalc)[0]
                self.L_condicional[A] = np.multiply(self.filhos[0].P_trs[A], 
                                  self.filhos[1].P_trs[A])
                self.P_trs[A] = np.dot(self.L_condicional[A], self.trs)



S0 = Especie(None, np.full(2,None), minha_priori = priori)
S6 = Especie(S0, np.full(2,None), meu_tempo = 20)
S1 = Especie(S6, np.array([2,2]), meu_tempo = 0.12)
S2 = Especie(S6, np.array([1,2]), meu_tempo = 0.6)
S8 = Especie(S0, np.full(2,None), meu_tempo = 30)
S3 = Especie(S8, np.array([0,3]), meu_tempo = 0.4)
S7 = Especie(S8, np.full(2,None), meu_tempo = 15)
S4 = Especie(S7, np.array([3,2]), meu_tempo = 0.4)
S5 = Especie(S7, np.array([2,1]), meu_tempo = 0.6)
S0.cria_L_condicional_vetor()
print(S1.trs)
print(S1.P_trs)
print(S0.P_trs)
print(S6.L_condicional)
print(S6.P_trs)
print(S0.L_condicional)
print(S0.L_arvore)


# testando o caso condicional com S1 = 'None' para apenas um codon
#S0 = Especie(S0, np.full(1,None), minha_priori = priori)
#S6 = Especie(S0, np.full(1,None), meu_tempo = 6)
#S1 = Especie(S6, np.array([None]), meu_tempo = 1)
#S2 = Especie(S6, np.array([0]), meu_tempo = 0.7 )
#S8 = Especie(S0, np.full(1,None), meu_tempo = 0.4)
#S3 = Especie(S8, np.array([0]), meu_tempo = 0.6)
#S7 = Especie(S8, np.full(1,None), meu_tempo = 0.9)
#S4 = Especie(S7, np.array([0]), meu_tempo = 0.4)
#S5 = Especie(S7, np.array([3]), meu_tempo = 0.3)


def prob_condicional(raiz, folha, val_codon, pos_codon): #um codon de cada vez
    # definimos fx1,x2,x3,x4 como S0.cria_L_condicional_vetor()
    # criando fx1,x2,x3,x4 para o valor fixado em val_codon
    pai_folha = folha.pai
    folha.recalc = np.full(folha.num_codon, True)
    while pai_folha is not raiz:
        no_recalc = pai_folha
        no_recalc.recalc[pos_codon] = True
        pai_folha = pai_folha.pai
    folha.valor[pos_codon] =  val_codon
    raiz.cria_L_condicional_vetor()
    P_conj = raiz.L_arvore
    P_conj = P_conj[0, pos_codon]
    P_soma = 0
    # ja temos todas os nos com exceção do que tem desconhecido calculado
    # e não precisamos recalcular para diferentes valores de S1
    #vetor de combinações
    #combin = np.tile(h,folha.num_codon).reshape(folha.num_codon,n_base), combinatorias extensas
    for i in range(n_base):
        folha.valor[pos_codon] = np.array([i])     #para mais de um codon, complica-se a escolha
        no_recalc.cria_L_condicional_vetor()
        raiz.L_condicional = np.multiply(raiz.filhos[0].P_trs[pos_codon], 
                                         raiz.filhos[1].P_trs[pos_codon])
        P_soma += np.dot(raiz.L_condicional, raiz.priori)
    P_cond = P_conj/P_soma
    return(P_cond)

### teste
#print(prob_condicional(S0, S1, 0))
#print(prob_condicional(S0,  S1, 1))
#print(prob_condicional(S0, S1, 2))
#print(prob_condicional(S0, S1, 3))
#print(prob_condicional(S0, S1, 0) + prob_condicional(S0, S1, 1) + 
# prob_condicional(S0, S1, 2) + prob_condicional(S0, S1, 3))

S0 = Especie(None, np.full(2,None), minha_priori = priori)
S6 = Especie(S0, np.full(2,None), meu_tempo = 20)
S1 = Especie(S6, np.array([1,2]), meu_tempo = 0.2)
S2 = Especie(S6, np.array([0,2]), meu_tempo = 0.7)
S8 = Especie(S0, np.full(2,None), meu_tempo = 30)
S3 = Especie(S8, np.array([0,3]), meu_tempo = 0.5)
S7 = Especie(S8, np.full(2,None), meu_tempo = 20)
S4 = Especie(S7, np.array([None,2]), meu_tempo = 0.25)
S5 = Especie(S7, np.array([1,3]), meu_tempo = 0.12)
print(prob_condicional(S0, S4, 0,0))
print(prob_condicional(S0, S4, 1,0))
print(prob_condicional(S0, S4, 2,0))
print(prob_condicional(S0, S4, 3,0))




#fazendo para outra arvore definida pelo exemplo
X7 = Especie(None, np.full(13, None), minha_priori = priori)
X5 = Especie(X7, np.full(13,None), meu_tempo = 10)
X6 = Especie(X7, np.full(13, None), meu_tempo = 12)
X1 = Especie(X5, np.array([1, None, 3, 1, 0, 2, 3, 1 , 2 ,3 ,1 ,2 ,0]), meu_tempo= 0.25)
X2 = Especie(X5, np.array([0, 1, 3, 0, 2, 1, 2, 3, 0, 2, 3, 1, 2]), meu_tempo = 12)
X3 = Especie(X6, np.array([0, 3, 2, 0, 1, 2, 3, 1, 0, 1 ,2, 3, 1]), meu_tempo = 2)
X4 = Especie(X6, np.array([3, 1, 1, 2, 0, 0, 1, 2, 3, 1, 0, 2, 1]), meu_tempo = 3)

### teste com X1 desconhecido
print(prob_condicional(X7, X1, 0,1))
print(prob_condicional(X7, X1, 1,1))
print(prob_condicional(X7, X1, 2,1))
print(prob_condicional(X7, X1, 3,1))
print(prob_condicional(X7, X1, 0,0) + prob_condicional(X7, X1, 1,0)+ 
      prob_condicional(X7, X1, 2,0) +
prob_condicional(X7, X1, 3,0))



class No:
    def __init__(self, indice, vizinhos = None, atributos = None, priori = None):
        self.indice = indice
        self.valor = atributos
        self.num_codon = self.valor.size
        self.L_condicional = np.full((self.num_codon,n_base), None)
        self.tempo = None
        self.P_trs = np.full((self.num_codon,n_base), None)
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
            print('No {} e {} não são vizinhos'.format(self.indice, indice_viz))
            
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
            
        elif(raiz.filhos is not None):     #eh interno
            for filho in raiz.filhos:
                no = No(filho, None, filho.valor, filho.priori)   #adicionando valores da arvore no n�
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
      if np.max(erros) > e:
        return(self.maxim_L_vetor(e,it))
      else:
        return([erros,it])
      
    def maxim_L(self,especie_1,especie_aresta):
      v = self.no(especie_1).retorna_peso(especie_aresta)
      q = math.exp(-v)
      p = 1 - q
      K = self.get_num_aresta()
      self.L_condicional_g_vetor(especie_1,especie_aresta)
      A = np.dot(np.multiply(self.no(especie_1).L_condicional,
                             self.no(especie_aresta).L_condicional)
      ,self.no(especie_1).priori)
      B = np.multiply(np.dot(self.no(especie_1).L_condicional,self.no(especie_1).priori),
      np.dot(self.no(especie_aresta).L_condicional,self.no(especie_aresta).priori))
      v_atual = -(np.log(1 - ((1/K)*(np.sum((p*B)/((A*q) + (B*p)))))))
      self.muda_peso(especie_1,especie_aresta,v_atual)
      erro = abs(v_atual - v)
      return(erro)
      
    
    def L_condicional_g_vetor(self,especie_1, especie_aresta, especie_2 = None):
        self.L_condicional_g(especie_1, especie_aresta, especie_2)
    
    def L_condicional_g(self,especie_1, especie_aresta,especie_2):  
        if especie_2 == None:  #especie_1 � a que vai, especie_2 � a que vem, especie_aresta 
            for especie in self.no(especie_1).indices_vizinhos():
                self.L_condicional_g(especie, especie_aresta,especie_1)
                self.no(especie).P_trs = np.dot(self.no(especie).L_condicional,self.no(especie).trs)
            viz = np.array(self.no(especie_1).indices_vizinhos())
            viz_temp = viz[viz == especie_aresta]     #aresta de interesse isolada
            viz = viz[viz != especie_aresta]
            if not (viz.size):
              self.no(especie_1).L_condicional = especie_1.L_condicional
            else:
              self.no(especie_1).L_condicional = np.multiply(self.no(viz[0]).P_trs,self.no(viz[1]).P_trs)
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
                    self.no(especie).P_trs = np.dot(self.no(especie).L_condicional,self.no(especie).trs)
                viz = viz[viz != especie_2]
                self.no(especie_1).L_condicional = np.multiply(self.no(viz[0]).P_trs,self.no(viz[1]).P_trs)


#g = Transforma_arv(S0)
#g.transforma_grafo(S0)
#g.L_condicional_g_vetor(S8,S7)
#print(g.no(S1).L_condicional)
#print(g.no(S7).L_condicional)
#print(g.no(S8).L_condicional)
#print(g.maxim_L_vetor((10**(-6))))
#g.no(S5).retorna_peso(S7)

### posteriormente criar classe do tipo "arvore"
### predição só use valores temporarios (val.temp)
### resolver problema na predição (slice e recalc na raiz) e testar para caso com 2 e 3 codons
### foco em apresentar
