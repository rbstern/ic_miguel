import numpy as np
import ML_TFK81 as fels
import matplotlib.pyplot as plt
import scipy.stats as scp

# 0 para T
# 1 para C
# 2 para A
# 3 para G

np.random.seed(0)
n_base = 4
priori = np.full((n_base, 1), 1/n_base)
taxa = 0.5

def value_gen(n_site):
    valor =  np.zeros(n_site)
    cumsum = np.cumsum(priori)
    for i in range(valor.size):
        u = np.random.uniform()
        if cumsum[0] <= u < cumsum[1]:
            valor[i] = 1
        elif cumsum[1] <= u < cumsum[2]:
            valor[i] = 2
        elif cumsum[2] <= u < cumsum[3]:
            valor[i] = 3
    return valor            

x = np.linspace(scp.expon.ppf(0.01, scale = 1),
                scp.expon.ppf(0.99, scale = 1), 100)
plt.plot(x, scp.expon.pdf(x, scale = 1),
         'r-', lw = 5, alpha = 0.6, label = 'expond_pdf')
plt.xlabel('x')
plt.ylabel('density')
plt.title('exponential(0.625) density')
plt.show()
        

class Especie_sim:
    def __init__(self, meu_pai, minha_priori, valor, meu_tempo = 0, taxa_subst = 1):
        self.filhos = []
        self.pai = meu_pai
        self.valor = valor
        self.taxa = taxa_subst
        if(meu_pai):
            meu_pai.filhos.append(self)
            self.priori = meu_pai.priori
            self.tempo = meu_tempo
            self.trs = self.cria_transicao()
            self.num_site = self.pai.num_site
        else:
            self.num_site = self.valor.size
            self.priori = minha_priori
    
    def cria_transicao(self):
        P = np.identity(n_base)*(np.exp(-self.taxa*self.tempo))
        for i in range(P[0,].size):
            P[i,] += ((1 - np.exp(-self.taxa*self.tempo))*self.priori[i,0])
        return P.transpose()
    
    def gera_val(self, site):
        pai_val = int(self.pai.valor[site])
        trs_sort = self.trs[pai_val]
        cumsum = np.cumsum(trs_sort)
        u = np.random.uniform()
        if cumsum[0] <= u < cumsum[1]:
            self.valor[site] = 1
        elif cumsum[1] <= u < cumsum[2]:
            self.valor[site] = 2
        elif cumsum[2] <= u < cumsum[3]:
            self.valor[site] = 3
    
    def simula_arv_starting(self):
        self.simula_arv
    
    def simula_arv(self, dici_base = {}):
        if(not(self.filhos)):
            self.valor = np.zeros(self.num_site)
            for i in range(self.num_site):
                self.gera_val(i)
            dici_base[self] = self.valor
        
        elif(not(self.pai)):
            for filho in self.filhos:
                filho.simula_arv()
            return dici_base
                
        else:
            self.valor = np.zeros(self.num_site)
            for i in range(self.num_site):
                self.gera_val(i)
            for filho in self.filhos:
                filho.simula_arv()

# gerando as bases e topologia original a ser predita
n_site = 2500
val_inicial = value_gen(n_site)
S0 = Especie_sim(meu_pai = None, minha_priori = priori, valor = val_inicial, taxa_subst = taxa)
S6 = Especie_sim(meu_pai = S0, minha_priori = priori, valor = None, 
                 meu_tempo = np.random.exponential(scale = 1), taxa_subst = taxa)
S8 = Especie_sim(meu_pai = S0, minha_priori = priori, valor = None,
                 meu_tempo = np.random.exponential(scale = 1), taxa_subst = taxa)
S1 = Especie_sim(meu_pai = S6, minha_priori = priori, valor = None,
                 meu_tempo = np.random.exponential(scale = 1), taxa_subst = taxa)
S2 = Especie_sim(meu_pai = S6, minha_priori = priori, valor = None,
                 meu_tempo = np.random.exponential(scale = 1), taxa_subst = taxa)
S3 = Especie_sim(meu_pai = S8, minha_priori = priori, valor = None,
                 meu_tempo = np.random.exponential(scale = 1), taxa_subst = taxa)
S4 = Especie_sim(meu_pai = S8, minha_priori = priori, valor = None,
                 meu_tempo = np.random.exponential(scale = 1), taxa_subst = taxa)
dici_base = S0.simula_arv()

num_0 = 0
num_1 = 0
num_2 = 0
num_3 = 0

for key in dici_base:
    num_0 += (dici_base[key] == 0).sum()
    num_1 += (dici_base[key] == 1).sum()
    num_2 += (dici_base[key] == 2).sum()
    num_3 += (dici_base[key] == 3).sum()

print(num_0/(2500*4))
print(num_1/(2500*4))
print(num_2/(2500*4))
print(num_3/(2500*4))



x = np.linspace(scp.expon.ppf(0.01, scale = 3),
                scp.expon.ppf(0.99, scale = 3), 100)
plt.plot(x, scp.expon.pdf(x, scale = 3),
         'r-', lw = 5, alpha = 0.6, label = 'expond_pdf')
plt.xlabel('x')
plt.ylabel('density')
plt.title('exponential(1/3) density')
plt.show()

# criando varias topologias na mao
def otimin_arv(raiz):
    g = fels.Grafo()
    g.transforma_grafo(raiz)
    g.maxim_L_vetor((10**(-4)))
    g.muda_tempo_arv()
    raiz.modifica_L_condicional()
    veros = np.sum(np.log(raiz.L_arvore))
    return veros

def display_veros(lista_raiz):
    for i in range(len(lista_raiz)):
        veros = otimin_arv(lista_raiz[i])
        print("Topology{}:".format(i + 1), veros)

lista_raiz = []

# topologia 1
S0_top1 = fels.Especie(None, np.full(n_site, None), priori)
S6_top1 = fels.Especie(S0_top1, np.full(n_site, None), 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S8_top1 = fels.Especie(S0_top1, np.full(n_site, None), 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S1_top1 = fels.Especie(S6_top1, dici_base[S1], 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S2_top1 = fels.Especie(S6_top1, dici_base[S2],
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S3_top1 = fels.Especie(S8_top1, dici_base[S3],
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S4_top1 = fels.Especie(S8_top1, dici_base[S4],
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)

# calculando verossimilhança
S0_top1.cria_L_condicional_vetor()
print(np.sum(np.log(S0_top1.L_arvore)))
lista_raiz.append(S0_top1)


# topologia 2
S0_top2 = fels.Especie(None, np.full(n_site, None), priori)
S2_top2 = fels.Especie(S0_top2, dici_base[S2], 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S8_top2 = fels.Especie(S0_top2, np.full(n_site, None), 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S1_top2 = fels.Especie(S8_top2, dici_base[S1], 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S7_top2 = fels.Especie(S8_top2, np.full(n_site, None), 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S3_top2 = fels.Especie(S7_top2, dici_base[S3],
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S4_top2 = fels.Especie(S7_top2, dici_base[S4],
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)

# calculando verossimilhança
S0_top2.cria_L_condicional_vetor()
print(np.sum(np.log(S0_top2.L_arvore)))
lista_raiz.append(S0_top2)



# topologia 3
S0_top3 = fels.Especie(None, np.full(n_site, None), priori, taxa_subst = taxa)
S1_top3 = fels.Especie(S0_top3, dici_base[S1], 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S8_top3 = fels.Especie(S0_top3, np.full(n_site, None), 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S2_top3 = fels.Especie(S8_top3, dici_base[S2], 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S7_top3 = fels.Especie(S8_top3, np.full(n_site, None), 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S3_top3 = fels.Especie(S7_top3, dici_base[S3],
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S4_top3 = fels.Especie(S7_top3, dici_base[S4],
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)

# calculando verossimilhança
S0_top3.cria_L_condicional_vetor()
print(np.sum(np.log(S0_top3.L_arvore)))
lista_raiz.append(S0_top3)


# topologia 4
S0_top4 = fels.Especie(None, np.full(n_site, None), priori, taxa_subst = taxa)
S4_top4 = fels.Especie(S0_top4, dici_base[S4], 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S8_top4 = fels.Especie(S0_top4, np.full(n_site, None), 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S3_top4 = fels.Especie(S8_top4, dici_base[S3], 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S7_top4 = fels.Especie(S8_top4, np.full(n_site, None), 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S1_top4 = fels.Especie(S7_top4, dici_base[S1],
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S2_top4 = fels.Especie(S7_top4, dici_base[S2],
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)

# calculando verossimilhança
S0_top4.cria_L_condicional_vetor()
print(np.sum(np.log(S0_top4.L_arvore)))
lista_raiz.append(S0_top4)


# topologia 5
S0_top5 = fels.Especie(None, np.full(n_site, None), priori)
S3_top5 = fels.Especie(S0_top5, dici_base[S3], 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S8_top5 = fels.Especie(S0_top5, np.full(n_site, None), 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S4_top5 = fels.Especie(S8_top5, dici_base[S4], 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S7_top5 = fels.Especie(S8_top5, np.full(n_site, None), 
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S1_top5 = fels.Especie(S7_top5, dici_base[S1],
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
S2_top5 = fels.Especie(S7_top5, dici_base[S2],
                       meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)

# calculando verossimilhança
S0_top5.cria_L_condicional_vetor()
print(np.sum(np.log(S0_top5.L_arvore)))
lista_raiz.append(S0_top5)

display_veros(lista_raiz)








        




            
                
                
                
                
                
                
                
            
            