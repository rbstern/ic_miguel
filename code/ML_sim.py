import numpy as np
import ML_TFK81 as fels
import matplotlib.pyplot as plt
import scipy.stats as scp
from itertools import combinations
import time
from Bio import Phylo
from io import StringIO

# 0 para T
# 1 para C
# 2 para A
# 3 para G

# graficos
x = np.linspace(scp.expon.ppf(0.01, scale = 1),
                scp.expon.ppf(0.99, scale = 1), 100)
plt.plot(x, scp.expon.pdf(x, scale = 1),
         'r-', lw = 5, alpha = 0.6, label = 'expond_pdf')
plt.xlabel('x')
plt.ylabel('density')
plt.title('exponential(1) density')
plt.show()

x = np.linspace(scp.expon.ppf(0.01, scale = 2),
                scp.expon.ppf(0.99, scale = 2), 100)
plt.plot(x, scp.expon.pdf(x, scale = 2),
         'r-', lw = 5, alpha = 0.6, label = 'expond_pdf')
plt.xlabel('x')
plt.ylabel('density')
plt.title('exponential(1/2) density')
plt.show()

np.random.seed(0)
n_base = 4
priori = np.full((n_base, 1), 1/n_base)
taxa = 1  

# all functions 

def list_comb(folhas_atuais):
    comb = list(combinations(np.arange(len(folhas_atuais)), 2))
    new_comb = []
    for i in range(len(comb)):
        new_comb.append(list(comb[i]))
    return new_comb

def gera_topologias(folhas_atuais):
    if len(folhas_atuais) == 1:
        tops.append(folhas_atuais)
        return None
    escolhas = list_comb(folhas_atuais)
    for i in range(len(escolhas)):
        pares = escolhas[i]
        pares_aux = np.arange(len(folhas_atuais))
        mask = np.ones(len(folhas_atuais), bool)
        mask[pares] = False
        pares_aux = list(pares_aux[mask])
        if list(np.array(folhas_atuais)[pares_aux]) == []:
            folhas_atuais_n = ["(" + str(folhas_atuais[pares[0]]) + "," + str(folhas_atuais[pares[1]])+ ")"]
            gera_topologias(folhas_atuais_n)
        else:
            folhas_atuais_n = list(np.array(folhas_atuais)[pares_aux])
            string = "(" + str(folhas_atuais[pares[0]]) + "," + str(folhas_atuais[pares[1]]) + ")"
            folhas_atuais_n.append(string)
            gera_topologias(folhas_atuais_n)

def modif_type(tops):
    for i in range(len(tops)):
        tops[i][0] = eval(tops[i][0])
    return tops

def read_top(top, especie_pai = None):
    if not(especie_pai):
        S0 = fels.Especie(None, np.full(n_site, None), priori, taxa_subst = taxa)
        top1 = top[0]
        top2 = top[1]
        read_top(top1, S0)
        read_top(top2, S0)
        return S0
        
    elif isinstance(top, tuple):
        S = fels.Especie(especie_pai, np.full(n_site, None), 
                         meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
        top1 = top[0]
        top2 = top[1]
        read_top(top1, S)
        read_top(top2, S)
    
    else:
        valor = dici_base[list(dici_base)[top]]
        S = fels.Especie(especie_pai, meu_valor = valor, 
                         meu_tempo = np.random.exponential(scale = 2), taxa_subst = taxa)
                 
def read_all_tops(tops):
    if isinstance(tops[0][0], str):
        tops = modif_type(tops)
    fels_roots = []
    for i in range(len(tops)):
        top = tops[i][0]
        S = read_top(top)
        fels_roots.append(S)
    return fels_roots

def change_time(especie, exp_par):
    if not(especie.pai):
        especie.filhos[0].tempo = np.random.exponential(scale = exp_par)
        especie.filhos[1].tempo = np.random.exponential(scale = exp_par)
        especie.filhos[0].trs =  especie.filhos[0].cria_transicao()
        especie.filhos[1].trs =  especie.filhos[1].cria_transicao()
        change_time(especie.filhos[0], exp_par)
        change_time(especie.filhos[1], exp_par)
        return especie
    
    elif(not(especie.filhos)):
        especie.tempo = np.random.exponential(scale = exp_par)
        especie.trs = especie.cria_transicao()
        
    else:
        especie.filhos[0].tempo = np.random.exponential(scale = exp_par)
        especie.filhos[1].tempo = np.random.exponential(scale = exp_par)
        especie.filhos[0].trs =  especie.filhos[0].cria_transicao()
        especie.filhos[1].trs =  especie.filhos[1].cria_transicao()
        change_time(especie.filhos[0], exp_par)
        change_time(especie.filhos[1], exp_par)
        

def change_all_time(roots, exp_par):
    for i in range(len(roots)):
        roots[i] = change_time(roots[i], exp_par)
    return roots
        
        

def permute(top):
    new_top = (top[1], top[0])
    return new_top
    
def top_pop(tops):
    i = 0
    while i < len(tops):
        if (isinstance(tops[i][0][0], tuple)) and (isinstance(tops[i][0][1], tuple)):
            permutation = [permute(tops[i][0])]
            if permutation in tops:
                tops.remove(permutation)
        i += 1  
    return tops

                    
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

        

class Especie_sim:
    def __init__(self, meu_pai, minha_priori, valor, taxa_subst = 1):
        self.filhos = []
        self.pai = meu_pai
        self.valor = valor
        self.taxa = taxa_subst
        if(meu_pai):
            meu_pai.filhos.append(self)
            self.priori = meu_pai.priori
            self.gera_tempo()
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
    
    def gera_tempo(self, scale = 1):
        self.tempo = 1

    def simula_arv(self, dici_filhos = {}):    
        if(not(self.pai)):
            for filho in self.filhos:
                filho.simula_arv()
            return dici_filhos
                
        self.valor = np.zeros(self.num_site)
        for i in range(self.num_site):
            self.gera_val(i)
        for filho in self.filhos:
            filho.simula_arv()
        if(not(self.filhos)):
            dici_filhos[self] = self.valor

# gerando as bases e topologia original a ser predita
n_site = 1500
val_inicial = value_gen(n_site)
S0 = Especie_sim(meu_pai = None, minha_priori = priori, valor = val_inicial, taxa_subst = taxa)
S6 = Especie_sim(meu_pai = S0, minha_priori = priori, valor = None, taxa_subst = taxa)
S8 = Especie_sim(meu_pai = S0, minha_priori = priori, valor = None, taxa_subst = taxa)
S1 = Especie_sim(meu_pai = S6, minha_priori = priori, valor = None, taxa_subst = taxa)
S2 = Especie_sim(meu_pai = S6, minha_priori = priori, valor = None, taxa_subst = taxa)
S3 = Especie_sim(meu_pai = S8, minha_priori = priori, valor = None, taxa_subst = taxa)
S4 = Especie_sim(meu_pai = S8, minha_priori = priori, valor = None, taxa_subst = taxa)
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

print(num_0/(n_site*4))
print(num_1/(n_site*4))
print(num_2/(n_site*4))
print(num_3/(n_site*4))


# criando e avaliando varias topologias
def display_veros(lista_raiz):
    veros_array = np.zeros(len(lista_raiz))
    for i in range(len(lista_raiz)):
        # before optimization
        lista_raiz[i].cria_L_condicional_vetor()
        veros = fels.otimin_arv(lista_raiz[i])
        veros_array[i] = veros
    return veros_array

# função para plotar todas as arvores e armazena-las na pasta certa
def plot_all_trees(tops, indices):
    for i in range(len(indices)):
        top_str = str(tops[indices[i]][0])
        topologia = Phylo.read(StringIO(top_str), 'newick')
        Phylo.draw(topologia)
        plt.title("Topologia {}".format(indices[i] + 1))
        plt.savefig("topologia{}".format(indices[i] + 1))


# gerando topologias
n_filhos = 4
tops = []
indices  = [0, 1, 2, 13, 14]
folhas = list(np.arange(n_filhos))
gera_topologias(folhas)
tops = modif_type(tops)
tops = top_pop(tops)          


# plotando as variaveis que mais apareceram nos graficos
plot_all_trees(tops, indices)

# fazendo varias simulações de acordo com o nsim
inicio = time.time()
nsim = 150
num_max = 5
max_data = np.zeros((nsim, num_max))
start = True

for i in range(nsim):
    if start == True:
        roots = read_all_tops(tops)
        veros_array = display_veros(roots)
        max_tops = veros_array.argsort()[-num_max:][::-1]
        max_tops += 1
        max_data[i] = np.array(max_tops)
        start = False
    else:
        roots = change_all_time(roots, 2)
        veros_array = display_veros(roots)
        max_tops = veros_array.argsort()[-num_max:][::-1]
        max_tops += 1
        max_data[i] = np.array(max_tops)
    print(i + 1)
fim = time.time()
tempo = (fim - inicio)/60
print("O tempo de execução do programa para", nsim, "simulações foi de", tempo, "minutos")


# plotando o grafico 
fig ,axs = plt.subplots(2,3)
fig.suptitle("Porcentagem de topologias em cada posição dos 5 maiores valores",
             fontsize=10)

for i in range(num_max):
    if 0 <= i <= 2:
        (unique, counts) = np.unique(max_data[:, i], return_counts=True)
        frequencies = np.asarray((unique, counts))
        # convertendo para porcentagem
        counts = counts/nsim
        percentages = np.asarray((unique, counts))
        axs[0, i].bar(percentages[0], height = percentages[1], color = 'lightblue', 
           edgecolor='blue')
        axs[0, i].title.set_text("Posição {}".format(i + 1))
        axs[0, i].set_xlabel("Topologias")
        axs[0, i].set_ylabel("Porcentagem")
        axs[0, i].set_ylim([0, 0.275])
    else:
        (unique, counts) = np.unique(max_data[:, i], return_counts=True)
        frequencies = np.asarray((unique, counts))
        # convertendo para porcentagem
        counts = counts/nsim
        percentages = np.asarray((unique, counts))
        axs[1, i - 3].bar(percentages[0], height = percentages[1], color = 'lightblue',
           edgecolor='blue')
        axs[1, i - 3].title.set_text("Posição {}".format(i + 1))
        axs[1, i - 3].set_xlabel("Topologias")
        axs[1, i - 3].set_ylabel("Porcentagem")
        axs[1, i - 3].set_ylim([0, 0.275])
axs[1,2].axis('off')
plt.show()

# plotando as arvores







# forma mais automatizada para simulações
# =============================================================================
# 
# def all_sim(n_site, nsim, num_max):
#   times = np.zeros(len(n_site))
#   resultados = []
#   for i in range(len(n_site)):
#     val_inicial = value_gen(n_site[i])
#     S0 = Especie_sim(meu_pai = None, minha_priori = priori, valor = val_inicial, taxa_subst = taxa)
#     S6 = Especie_sim(meu_pai = S0, minha_priori = priori, valor = None, taxa_subst = taxa)
#     S8 = Especie_sim(meu_pai = S0, minha_priori = priori, valor = None, taxa_subst = taxa)
#     S1 = Especie_sim(meu_pai = S6, minha_priori = priori, valor = None, taxa_subst = taxa)
#     S2 = Especie_sim(meu_pai = S6, minha_priori = priori, valor = None, taxa_subst = taxa)
#     S3 = Especie_sim(meu_pai = S8, minha_priori = priori, valor = None, taxa_subst = taxa)
#     S4 = Especie_sim(meu_pai = S8, minha_priori = priori, valor = None, taxa_subst = taxa)
#     new_base = S0.simula_arv()
#     # fazendo varias simulações de acordo com o nsim
#     inicio = time.time()  
#     max_data = np.zeros((nsim, num_max))
#     start = True
#     for j in range(nsim):
#       if start == True:
#         print(new_base[S1])
#         roots = read_all_tops(tops, new_base, n_site)
#         print(roots[0].filhos[0].filhos[0].valor)
#         veros_array = display_veros(roots)
#         max_tops = veros_array.argsort()[-num_max:][::-1]
#         max_tops += 1
#         max_data[i] = np.array(max_tops)
#         start = False
#       else:
#         roots = change_all_time(roots, 2)
#         veros_array = display_veros(roots)
#         max_tops = veros_array.argsort()[-num_max:][::-1]
#         max_tops += 1
#         max_data[j] = np.array(max_tops)
#     resultados.append(max_data)
#     fim = time.time()
#     tempo = (fim - inicio)/60
#     times[i] = tempo
#   print(times)
#   return resultados
#   
# def plot_all(results, nsim, num_max):
#   for i in range(len(results)):
#     fig ,axs = plt.subplots(2,3, figsize = (18 ,12))
#     fig.suptitle("Porcentagem de topologias em cada posição dos 5 maiores valores",
#              fontsize=10)
#     plt.setp(axs, xticks = range(1, len(tops) + 1, 2))
#     for j in range(num_max):
#       if 0 <= j <= 2:
#         (unique, counts) = np.unique(results[i][:, j], return_counts=True)
#         frequencies = np.asarray((unique, counts))
#         # convertendo para porcentagem
#         counts = counts/nsim
#         percentages = np.asarray((unique, counts))
#         axs[0, j].bar(percentages[0], height = percentages[1], color = 'lightblue', 
#            edgecolor='blue')
#         axs[0, j].title.set_text("Posição {}".format(j + 1))
#         axs[0, j].set_xlabel("Topologias")
#         axs[0, j].set_ylabel("Porcentagem")
#         axs[0, j].set_ylim([0, 0.295])
#       else:
#         (unique, counts) = np.unique(results[i][:, j], return_counts=True)
#         frequencies = np.asarray((unique, counts))
#         # convertendo para porcentagem
#         counts = counts/nsim
#         percentages = np.asarray((unique, counts))
#         axs[1, j - 3].bar(percentages[0], height = percentages[1], color = 'lightblue',
#            edgecolor='blue')
#         axs[1, j - 3].title.set_text("Posição {}".format(j + 1))
#         axs[1, j - 3].set_xlabel("Topologias")
#         axs[1, j - 3].set_ylabel("Porcentagem")
#         axs[1, j - 3].set_ylim([0, 0.295])
#     axs[1,2].axis('off')
#     plt.show()
# 
# =============================================================================
