import numpy as np
import math
import matplotlib as mat
import matplotlib.pyplot as plt
import scipy
import random
from random import randrange
import pandas as pd
import pylab as pl
import seaborn as sb

adj = np.matrix([[0,2,2,1,2,2,1,1,2,2,2,2,2,2,1,1,1,2,2,1],
                [2,0,2,2,1,1,2,1,1,1,1,1,1,2,1,1,1,1,2,2],
                [2,2,0,1,2,2,2,2,1,1,2,1,2,2,2,1,1,3,1,2],
                [1,2,1,0,2,2,1,1,1,2,2,2,3,2,2,2,2,3,1,1],
                [2,1,2,2,0,3,3,1,2,2,2,3,3,1,2,1,2,1,1,2],
                [2,1,2,2,3,0,1,2,1,2,1,1,2,3,1,2,2,2,2,2],
                [1,2,2,1,3,1,0,1,2,2,2,1,2,3,2,2,2,2,2,1],
                [1,1,2,1,1,2,1,0,2,2,2,2,2,2,2,1,2,1,2,1],
                [2,1,1,1,2,1,2,2,0,2,1,2,3,2,1,2,2,3,1,2],
                [2,1,1,2,2,2,2,2,2,0,1,1,1,1,2,1,1,3,2,1],
                [2,1,2,2,2,1,2,2,1,1,0,2,1,1,1,1,2,1,2,1],
                [2,1,1,2,3,1,1,2,2,1,2,0,1,3,2,2,1,2,2,2],
                [2,1,2,3,3,2,2,2,3,1,1,1,0,2,2,2,1,2,3,1],
                [2,2,2,2,1,3,3,2,2,1,1,3,2,0,2,1,2,2,1,1],
                [1,1,2,2,2,1,2,2,1,2,1,2,2,2,0,1,1,2,2,2],
                [1,1,1,2,1,2,2,1,2,1,1,2,2,1,1,0,1,1,1,2],
                [1,1,1,2,2,2,2,2,2,1,2,1,1,2,1,1,0,2,2,2],
                [2,1,3,3,1,2,2,1,3,3,1,2,2,2,2,1,2,0,2,2],
                [2,2,1,1,1,2,2,2,1,2,2,2,3,1,2,1,2,2,0,2],
                [1,2,2,1,2,2,1,1,2,1,1,2,1,1,2,2,2,2,2,0]])
n_amn = adj[1].size
#Analisando cada array
#Combinações realizadas pelo computador
def Comparing_l (l,pr1,pr2):
    sum_list = []
    n1 = pr1.size
    n2 = pr2.size
    for i in range (n1):
        if l+i <= n1:
            npr1 = pr1[i:l+i-1]
            for j in range (n2):
                if l+j <= n2:
                    npr2 = pr2[j:l+j-1]
                    npr = adj[npr2,npr1]
                    nsum = npr.sum()
                    sum_list.append(nsum)
    sum_array = np.array(sum_list)
    return sum_array

def P_calc (l,values_arrays,pr1,pr2):
    n = values_arrays.size
    n_prob = np.zeros(n)
    sum_array = Comparing_l(l,pr1,pr2)
    p_tot = sum_array.size
    for i in range (p_tot):
        value = sum_array[i]
        n_prob[value] += 1
    n_prob = n_prob * (1/p_tot)
    return n_prob

#Distribuição multinomial e P_tot
from scipy.stats import multinomial

def Ps (pr1,pr2):
    count0 = 0
    count1 = 0
    count2 = 0
    count3 = 0
    for i in range(n_amn):
        for j in range(n_amn):
            if adj[i,j] == 0:
                prop0 = ((pr1==i).sum()/(pr1.size))*((pr2==j).sum()/(pr2.size))
                count0 += prop0
            elif adj[i,j] == 1:
                prop1 = ((pr1==i).sum()/(pr1.size))*((pr2==j).sum()/(pr2.size))
                count1 += prop1
            elif adj[i,j] == 2:
                prop2 = ((pr1==i).sum()/(pr1.size))*((pr2==j).sum()/(pr2.size))
                count2 += prop2
            elif adj[i,j] == 3:
                prop3 = ((pr1==i).sum()/(pr1.size))*((pr2==j).sum()/(pr2.size))
                count3 += prop3
    return count0,count1,count2,count3


def P_tot(l,n_array,pr1,pr2):
    p0,p1,p2,p3 = Ps(pr1,pr2)
    va = multinomial(l,[p0,p1,p2,p3])
    lista_soma = []
    for h in range(0,n_array.size):
        n = n_array[h]
        soma = 0
        for k in range(0,int(n/3)+1):
            for j in range(0,int(n/2)+1):
                i = 0
                while i+2*j+3*k <= n:
                    if i + 2*j + 3*k == n:
                        soma += va.pmf([l-(i+j+k),i,j,k])
                    i += 1
        lista_soma.append(soma)
    n_soma = np.array(lista_soma)
    return n_soma


def P_rat(l,n_array,pr1,pr2):
    lista = []
    prob_total = P_tot(l,n_array,pr1,pr2)
    prob_calc = P_calc(l,n_array,pr1,pr2)
    for i in range (n_array.size):
        if prob_calc[i] != 0:
            lista.append(prob_calc[i]/prob_total[i])
        else:
            lista.append(0)
    rat_array = np.array(lista)
    return rat_array
        
    
                    
print(adj[[1,2],[1,3]])   
#
#pr1 = np.array([12,13,1,2,10,19,12,14,13,12,0])
#pr2 = np.array([13,14,15,16,17,18,19,19,12,3,15,5,7])

#tabela de frequencias e histograma usando pandas
#tabela de frequencia
#tab_freq = pd.crosstab(index = Comparing_l(5,pr1,pr2), columns ="counts")
#tab_freq = pd.DataFrame(tab_freq)
#print(tab_freq)

#histograma
#df = (Comparing_l(10,pr1,pr2))
#hist = df.hist(bins = 20 )

#histograma comparando P_calc com P_tot
#x = np.arange(0,30,1)
#fig, ax = plt.subplots(1,1)
#n, bins , patches = ax.hist( Comparing_l(5,pr1,pr2), bins=20,
#    weights = np.zeros_like(Comparing_l(5,pr1,pr2))+1/Comparing_l(5,pr1,pr2).size)
#ax.plot(x, P_tot(5,x,pr1,pr2),'--')
#ax.set_xlabel("Mudanças de bases")
#ax.set_ylabel("Frequência")
#ax.set_xlim([0,15])
#ax.set_title("Histograma de frequencia/Comparação com densidade")
#plt.show()


#grafico de linhas comparando P_calc com P_tot
#x = np.arange(0,100)
#plt.plot(x,P_calc(10,x,pr1,pr2),alpha=0.5,lw = 2, color = 'r',linestyle='--',marker='o')
#plt.plot(x,P_tot(10,x,pr1,pr2),alpha = 0.5, lw = 2,
#         label = "Comparação densidade e observados", color = 'b',linestyle='--',marker='o')
#plt.xlabel("Base diferences")
#plt.ylabel("Probability")



#print(P_tot(10,np.array([8]),pr1,pr2))


#plotando
x = np.arange(0,100,1)
#h = plt.plot(x,P_rat(10,x,pr1,pr2),lw = 2, alpha = 0.6, label='comparação',linestyle = '--'
#             , color = 'b', marker = 'o')
#plt.show()
pr1_l = []
pr2_l = []
for i in range(0,100):
    pr1_l.append(randrange(19))
    pr2_l.append(randrange(19))
pr1 = np.array(pr1_l)
pr2 = np.array(pr2_l)

#plt.plot(range(10),linestyle = '--', marker = 'o', color = 'b')
x = np.arange(0,100,1)
h = plt.plot(x,P_rat(30,x,pr1,pr2),lw = 2, alpha = 0.6, label='comparação',linestyle = '--'
             , color = 'b', marker = 'o')
plt.yscale(value = "log")
plt.show()
