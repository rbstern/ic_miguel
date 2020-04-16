import numpy as np
import math
import time
from scipy.optimize import minimize

# Matrizes do numero de transições e transversoes observadas entre especies
# 0 - RATO 
# 1 - BOI
# 2 - GIBAO
# 3 - ORANGOTANGO
# 4 - GORILA
# 5 - CHIMPANZE
# 6 - HOMEM

# Nucleotideos de classe 1 --> 3ª posiçao no codon
# Nucleotideos de classe 2 --> 1ª e 2ª posiçoes no codon
num_especie = 7

# Quantidades amostrais
# Numero de transversoes da classe 1
V1 = np.array([[0, 82, 83, 85, 77, 79, 77],
               [82, 0, 71, 65, 67, 67, 67],
               [83, 71, 0, 34, 26, 26, 26],
               [85, 65, 34, 0, 18, 18, 20],
               [77, 67, 26, 18, 0, 4, 4],
               [79, 67, 26, 18, 4, 0, 2],
               [77, 67, 26, 20, 4, 2, 0]])

# Numero de transversoes de classe 2
V2 = np.array([[0, 91, 83, 90, 85, 86, 89],
               [91, 0, 69, 65, 72, 71, 70],
               [83, 69, 0, 18, 19, 18, 19],
               [90, 65, 18, 0, 15, 16, 15],
               [85, 72, 19, 15, 0, 5, 4],
               [86, 71, 18, 16, 5, 0, 3],
               [89, 70, 19, 15, 4, 3, 0]])

# Numero de transiçoes de classe 1
S1 = np.array([[0, 39, 53, 48, 46, 50, 51],
               [39, 0, 42, 44, 52, 61, 57],
               [53, 42, 0, 59, 59, 64, 58],
               [48, 44, 59, 0, 52, 60, 53],
               [46, 52, 59, 52, 0, 58, 52],
               [50, 61, 64, 60, 58, 0, 50],
               [51, 57, 58, 53, 52, 50, 0]])

# Numero de transiçoes de classe 2
S2 = np.array([[0, 68, 81, 81, 87, 79, 79],
               [68, 0, 80, 81, 93, 85, 86],
               [81, 80, 0, 57, 65, 61, 59],
               [81, 81, 57, 0, 64, 59, 55],
               [87, 93, 65, 64, 0, 28, 32],
               [79, 85, 61, 59, 28, 0, 24],
               [79, 86, 59, 55, 32, 24, 0]])

# Medidas estacionarias de cada base
# T = 0, C = 1, A = 2, G = 3
pi_1 = np.array([0.169, 0.429, 0.364, 0.038])
pi_2 = np.array([0.297, 0.267, 0.310, 0.126])
r_1 = 232
r_2 = 667
t_2 = 65

# Matriz de transição
def P_trans(t,beta, alpha, classe):
    if classe == 1:
        beta1 = beta
        alpha1 = alpha
        pi_y = pi_1[0] + pi_1[1]
        pi_r = pi_1[2] + pi_1[3]
        part_1 = np.dot(np.array([[1,1,1,1]]).transpose(), np.array([[pi_1
                              [0], pi_1[1], pi_1[2], pi_1[3]]]))
        part_2 = np.exp(-(beta1* t))*np.dot(np.array([[1/pi_y, 
                           1/pi_y, -1/pi_r, -1/pi_r]]).transpose(), np.array([[pi_r*pi_1[0],
                                            pi_r*pi_1[1], -pi_y*pi_1[2], -pi_y*pi_1[3]]]))
        part_3 = np.exp(-t*(pi_y*beta1 + pi_r*alpha1))*np.dot(np.array([[0, 0, 
                       pi_1[3]/pi_r, -pi_1[2]/pi_r]]).transpose(), np.array([[0, 0, 1, -1]]))
        part_4 = np.exp(-t*(pi_y*alpha1 + pi_r*beta1))*np.dot(np.array([[pi_1[1]/pi_y,
                          -pi_1[0]/pi_y, 0, 0]]).transpose(), np.array([[1, -1, 0, 0]]))
        mat_trans = part_1 + part_2 + part_3 + part_4
    else:
        beta2 = beta
        alpha2 = alpha
        pi_y = pi_2[0] + pi_2[1]
        pi_r = pi_2[2] + pi_2[3]
        part_1 = np.dot(np.array([[1,1,1,1]]).transpose(), np.array([[pi_2
                              [0], pi_2[1], pi_2[2], pi_2[3]]]))
        part_2 = np.exp(-(beta2* t))*np.dot(np.array([[1/pi_y, 
                           1/pi_y, -1/pi_r, -1/pi_r]]).transpose(), np.array([[pi_r*pi_2[0],
                                        pi_r*pi_2[1], -pi_y*pi_2[2], -pi_y*pi_2[3]]]))
        part_3 = np.exp(-t*(pi_y*beta2 + pi_r*alpha2))*np.dot(np.array([[0, 0, 
                           pi_2[3]/pi_r, -pi_2[2]/pi_r]]).transpose(), np.array([[0, 0, 1, -1]]))
        part_4 = np.exp(-t*(pi_y*alpha2 + pi_r*beta2))*np.dot(np.array([[pi_2[1]/pi_y,
                           -pi_2[0]/pi_y, 0, 0]]).transpose(), np.array([[1, -1, 0, 0]]))
        mat_trans = part_1 + part_2 + part_3 + part_4     
    return mat_trans



# Medias e variancias
def V_barra (f, beta, especie, classe, t_esp):
    if classe == 1:
        pi_y = pi_1[0] + pi_1[1]
        pi_r = pi_1[2] + pi_1[3]   
        media = 2*f*r_1*(pi_y*pi_r)*(1 - np.exp(-2*beta*t_esp))
    else:
        pi_y = pi_2[0] + pi_2[1]
        pi_r = pi_2[2] + pi_2[3]
        if especie == 2:
            media = 2*f*r_2*(pi_y*pi_r)*(1 - np.exp(-2*beta*t_2))
        else:
            media = 2*f*r_2*(pi_y*pi_r)*(1 - np.exp(-2*beta*t_esp))
    return media


def S_barra (f, beta, alpha, especie, classe, t_esp):
    if classe == 1:
        pi_y = pi_1[0] + pi_1[1]
        pi_r = pi_1[2] + pi_1[3]
        if especie == 2:
            media = 2*f*r_1*(((pi_1[0]*pi_1[1])+(pi_1[2]*pi_1[3])) + 
                       (((pi_1[0]*pi_1[1]*pi_r)/pi_y + (pi_1[2]*pi_1[3]*pi_y)/pi_r)*
                       np.exp(-2*beta*t_2)) -
                       ((pi_1[0]*pi_1[1]/pi_y) * np.exp(-(2*t_2)*(alpha*pi_y + beta*pi_r))) -
                       ((pi_1[2]*pi_1[3]/pi_r) * np.exp(-(2*t_2)*(alpha*pi_r + beta*pi_y))))
        else:
            media = 2*f*r_1*(((pi_1[0]*pi_1[1])+(pi_1[2]*pi_1[3])) + 
                       (((pi_1[0]*pi_1[1]*pi_r)/pi_y + (pi_1[2]*pi_1[3]*pi_y)/pi_r)*
                       np.exp(-2*beta*t_esp)) -
                       ((pi_1[0]*pi_1[1]/pi_y) * np.exp(-(2*t_esp)*(alpha*pi_y + beta*pi_r))) -
                       ((pi_1[2]*pi_1[3]/pi_r) * np.exp(-(2*t_esp)*(alpha*pi_r + beta*pi_y))))
    else:
        pi_y = pi_2[0] + pi_2[1]
        pi_r = pi_2[2] + pi_2[3]
        if especie == 2:
            media = 2*f*r_2*(((pi_2[0]*pi_2[1])+(pi_2[2]*pi_2[3])) + 
                       (((pi_2[0]*pi_2[1]*pi_r)/pi_y + (pi_2[2]*pi_2[3]*pi_y)/pi_r)*
                       np.exp(-2*beta*t_2)) -
                       ((pi_2[0]*pi_2[1]/pi_y) * np.exp(-(2*t_2)*(alpha*pi_y + beta*pi_r))) -
                       ((pi_2[2]*pi_2[3]/pi_r) * np.exp(-(2*t_2)*(alpha*pi_r + beta*pi_y))))
        else:
            media = 2*f*r_2*(((pi_2[0]*pi_2[1])+(pi_2[2]*pi_2[3])) + 
                       (((pi_2[0]*pi_2[1]*pi_r)/pi_y + (pi_2[2]*pi_2[3]*pi_y)/pi_r)*
                       np.exp(-2*beta*t_esp)) -
                       ((pi_2[0]*pi_2[1]/pi_y) * np.exp(-(2*t_esp)*(alpha*pi_y + beta*pi_r))) -
                       ((pi_2[2]*pi_2[3]/pi_r) * np.exp(-(2*t_esp)*(alpha*pi_r + beta*pi_y))))
    return media



def V_var (f, beta, especie, classe, t_esp):
    if classe == 1:
        var = V_barra(f, beta, especie, classe, t_esp)*(1 - 
                   (V_barra(f, beta, especie, classe, t_esp)/r_1))
    else:
        var = V_barra(f, beta, especie, classe, t_esp)*(1 - 
                   (V_barra(f, beta, especie, classe, t_esp)/r_2))
    return var



def S_var (f, beta, alpha, especie, classe, t_esp):
    if classe == 1:
        var = S_barra(f, beta, alpha, especie, classe, t_esp)*(1 - 
                     (S_barra(f, beta, alpha, especie, classe, t_esp)/r_1))
        return var
    else:
        var = S_barra(f, beta, alpha, especie, classe, t_esp)*(1 - 
                     (S_barra(f, beta, alpha, especie, classe, t_esp)/r_2))
        return var
    
def cov_V_S(V_barra, S_barra, classe):
    if classe == 1:
       cov = -(V_barra*S_barra)/r_1
    else:
       cov = -(V_barra*S_barra)/r_2
    return cov
    

# Definindo transversoes e transiçoes amostrais
def V_k(classe, especie):
    soma = 0
    if classe == 1:
        j = especie + 1
        while j <= 6:
            soma += V1[especie,j]
            j += 1
        v_k = soma/(7 - especie)
        return v_k
    else:
        j = especie + 1
        while j <= 6:
            soma += V2[especie, j]
            j += 1
        v_k = soma/(7 - especie)
        return v_k

def S_k(classe, especie):
    soma = 0
    if classe == 1:
        j = especie + 1
        while j <= 6:
            soma += S1[especie,j]
            j += 1
        s_k = soma/(7 - especie)
        return s_k
    else:
        j = especie + 1
        while j <= 6:
            soma += S2[especie, j]
            j += 1
        s_k = soma/(7 - especie)
        return s_k


# Definindo o vetor de observações
D_tilde = np.array([[V_k(1,1), V_k(1,2),V_k(1,3), V_k(1,4), V_k(1,5),
V_k(1,6), S_k(1,1), S_k(1,2), S_k(1,3), S_k(1,4), S_k(1,5),
S_k(1,6), V_k(2,1), V_k(2,2), V_k(2,3), V_k(2,4), V_k(2,5), 
V_k(2,6), S_k(2,1), S_k(2,2), S_k(2,3), S_k(2,4), S_k(2,5),
S_k(2,6)]])
    


    
# Definindo a covarincia e variancia de um V_k^(i) e S_k(i)
# se tipo == 0 --> V_k^(i)
# se tipo == 1 --> S_k^(i)
def q_1(f, t, alpha, beta, i, j, especie_1, classe):
    if classe == 1:
        if especie_1 != 2:
            q = f*pi_1[i]*P_trans(t, beta, alpha, classe)[i, j]
            return q
        else:
            q = f*pi_1[i]*P_trans(t_2, beta, alpha, classe)[i, j]
            return q
    else:
        if especie_1 != 2:
            q = f*pi_2[i]*P_trans(t, beta, alpha, classe)[i, j]
            return q
        else:
            q = f*pi_2[i]*P_trans(t_2, beta, alpha, classe)[i, j]
            return q
        
def q_2(t_I, t_J, t_K, t_L, f , alpha, beta, i, j, k, l, 
        especie_1, especie_2, especie_3, especie_4, classe):
    if classe == 1:
        if especie_1 < especie_2 < especie_3 < especie_4:
            soma = 0
            for x in range(0, pi_1.size):
                for y in range(0, pi_1.size):
                    soma += pi_1[x]*P_trans(2*t_I - t_J, beta, alpha,classe)[x, i]*P_trans(t_J,
                                beta, alpha, classe)[x, j]*P_trans(t_J - t_K,
                                                   beta, alpha, classe)[x, y]*P_trans(t_K,
                                       beta, alpha, classe)[y, k]*P_trans(t_K, 
                                                          beta, alpha, classe)[y, l]
            q = f*soma
            return q
        elif especie_2 > especie_3 and especie_2 < especie_4:
            soma = 0
            for x in range(0, pi_1.size):
                for y in range(0, pi_1.size):
                    soma += pi_1[x]*P_trans(2*t_I - t_K, beta, alpha, classe)[x, i]*P_trans(t_K,
                                beta, alpha, classe)[x, k]*P_trans(t_K - t_J, 
                                                   beta, alpha, classe)[x, y]*P_trans(t_J,
                                       beta, alpha, classe)[y, j]*P_trans(t_J, 
                                                          beta, alpha, classe)[y, l]
            q = f*soma
            return q
        elif especie_2 == especie_3:
            soma = 0
            if j == k:
                for x in range(0, pi_1.size):
                    soma += pi_1[x]*P_trans(2*t_I - t_J, beta, alpha, classe)[x, i]*P_trans(t_J, 
                                beta, alpha, classe)[x, j]*P_trans(t_J, beta, alpha,
                                                   classe)[x, l]
            q = f*soma
            return q
        elif especie_2 == especie_4:
            soma = 0
            if j == l:
                for x in range(0, pi_1.size):
                    soma += pi_1[x]*P_trans(2*t_I - t_K, beta, alpha, classe)[x, i]*P_trans(t_K,
                                beta, alpha, classe)[x, k]*P_trans(t_K, beta, alpha,
                                                   classe)[x, l]
            q = f*soma
            return q
        elif especie_2 > especie_4:
            soma =  0
            for x in range(0, pi_1.size):
                for y in range(0, pi_1.size):
                    soma += pi_1[x]*P_trans(2*t_I - t_K, beta, alpha, classe)[x, i]*P_trans(t_L,
                                beta, alpha, classe)[x, l]*P_trans(t_K - t_L, beta, alpha,
                                                   classe)[x, y]*P_trans(t_L,
                                       beta, alpha, classe)[y, l]*P_trans(t_L, beta, alpha, 
                                                          classe)[y, j]
            q = f*soma
            return q
        elif especie_1 == especie_3 and especie_2 < especie_4:
            soma = 0
            if i == k:
                for x in range(0, pi_1.size):
                    soma += pi_1[x]*P_trans(2*t_I - t_J, beta, alpha, classe)[x, i]*P_trans(t_J, 
                                beta, alpha, classe)[x, j]*P_trans(t_J, beta, alpha, 
                                                   classe)[x, l]
            q = f*soma
            return q
                    
    else:
        if especie_1 < especie_2 < especie_3 < especie_4:
            soma = 0
            for x in range(0, pi_2.size):
                for y in range(0, pi_2.size):
                    soma += pi_2[x]*P_trans(2*t_I - t_J, beta, alpha, classe)[x, i]*P_trans(t_J,
                                beta, alpha, classe)[x, j]*P_trans(t_J - t_K, beta, alpha,
                                                   classe)[x, y]*P_trans(t_K, beta, alpha,
                                       classe)[y, k]*P_trans(t_K, beta, alpha, classe)[y, l]
            q = f*soma
            return q
        elif especie_2 > especie_3 and especie_2 < especie_4:
            soma = 0
            for x in range(0, pi_2.size):
                for y in range(0, pi_2.size):
                    soma += pi_2[x]*P_trans(2*t_I - t_K, beta, alpha, classe)[x, i]*P_trans(t_K,
                                beta, 
                                alpha, classe)[x, k]*P_trans(t_K - t_J, beta, alpha, 
                                             classe)[x, y]*P_trans(t_J,
                                       beta, alpha, classe)[y, j]*P_trans(t_J, beta, alpha,
                                                          classe)[y, l]
            q = f*soma
            return q
        elif especie_2 == especie_3:
            soma = 0
            if j == k:
                for x in range(0, pi_2.size):
                    soma += pi_2[x]*P_trans(2*t_I - t_J, beta, alpha, classe)[x, i]*P_trans(t_J, 
                                beta, alpha, classe)[x, j]*P_trans(t_J, beta, alpha, 
                                                   classe)[x, l]
            q = f*soma
            return q
        elif especie_2 == especie_4:
            soma = 0
            if j == l:
                for x in range(0, pi_2.size):
                    soma += pi_2[x]*P_trans(2*t_I - t_K, beta, alpha, classe)[x, i]*P_trans(t_K,
                                beta, alpha, classe)[x, k]*P_trans(t_K, 
                                                   beta, alpha, classe)[x, l]
            q = f*soma
            return q
        elif especie_2 > especie_4:
            soma =  0
            for x in range(0, pi_2.size):
                for y in range(0, pi_2.size):
                    soma += pi_2[x]*P_trans(2*t_I - t_K, beta, alpha, classe)[x, i]*P_trans(t_L,
                                beta, alpha, classe)[x, l]*P_trans(t_K - t_L, 
                                                   beta, alpha, classe)[x, y]*P_trans(t_L,
                                       beta, alpha, classe)[y, l]*P_trans(t_L, 
                                                          beta, alpha, classe)[y, j]
            q = f*soma
            return q
        elif especie_1 == especie_3 and especie_2 < especie_4:
            soma  = 0
            if i == k:
                for x in range(0, pi_2.size):
                    soma += pi_2[x]*P_trans(2*t_I - t_J, beta, alpha, classe)[x, i]*P_trans(t_J, 
                                beta, alpha, classe)[x, j]*P_trans(t_J, beta, alpha, 
                                                   classe)[x, l]
            q = f*soma
            return q


# q_1(f, t, alpha, beta, i, j, especie_1, classe)
# q_2(t_I, t_J, t_K, t_L, f , alpha, beta, i, j, k, l, 
# especie_1, especie_2, especie_3, especie_4, classe)
            
def cov(t_I, t_J, t_K, t_L, 
        especie_1, especie_2, especie_3, especie_4, classe, tipo1, tipo2,
        f, alpha, beta):
    comb_transv = [[0,2],[0,3],[1,2],[1,3],
                   [2,0],[2,1],[3,0],[3,1]]
    comb_transi = [[0,1],[1,0],[2,3],[3,2]]
    if tipo1 == tipo2 == 0:
        soma = 0
        if classe == 1:
            for x in comb_transv:
                i, j = x[0], x[1]
                for y in comb_transv:
                    k, l = y[0], y[1]
                    soma += (q_2(t_I, t_J, t_K, t_L, f, alpha, beta, 
                                 i, j, k, l, especie_1, especie_2, 
                                 especie_3, especie_4, classe) -(q_1(f, t_I, alpha, beta, 
                            i, j, especie_1, classe)*q_1(f, t_K, alpha, beta, 
                                                   k, l, especie_3, classe)))
            cov = r_1*soma
            return cov
        else:
            for x in comb_transv:
                i, j = x[0], x[1]
                for y in comb_transv:
                    k, l = y[0], y[1]
                    soma +=  (q_2(t_I, t_J, t_K, t_L, f, alpha, beta, 
                                 i, j, k, l, especie_1, especie_2, 
                                 especie_3, especie_4, classe) -(q_1(f, t_I, alpha, beta, 
                            i, j, especie_1, classe)*q_1(f, t_K, alpha, beta, 
                                                   k, l, especie_3, classe)))
            cov = r_2*soma
            return cov
    elif tipo1 == tipo2 == 1:
        soma = 0
        if classe == 1:
            for x in comb_transi:
                i, j = x[0], x[1]
                for y in comb_transi:
                    k, l = y[0], y[1]
                    soma +=  (q_2(t_I, t_J, t_K, t_L, f, alpha, beta, 
                                 i, j, k, l, especie_1, especie_2, 
                                 especie_3, especie_4, classe) -(q_1(f, t_I, alpha, beta, 
                            i, j, especie_1, classe)*q_1(f, t_K, alpha, beta, 
                                                   k, l, especie_3, classe)))
            cov = r_1*soma
            return cov
        else:
            for x in comb_transi:
                i, j = x[0], x[1]
                for y in comb_transi:
                    k, l = y[0], y[1]
                    soma +=  (q_2(t_I, t_J, t_K, t_L, f, alpha, beta, 
                                 i, j, k, l, especie_1, especie_2, 
                                 especie_3, especie_4, classe) -(q_1(f, t_I, alpha, beta, 
                            i, j, especie_1, classe)*q_1(f, t_K, alpha, beta, 
                                                   k, l, especie_3, classe)))
            cov = r_2*soma
            return cov
    elif tipo1 == 1 and tipo2 == 0:
        soma = 0
        if classe == 1:
            for x in comb_transi:
                i, j = x[0], x[1]
                for y in comb_transv:
                    k, l = y[0], y[1]
                    soma +=  (q_2(t_I, t_J, t_K, t_L, f, alpha, beta, 
                                 i, j, k, l, especie_1, especie_2, 
                                 especie_3, especie_4, classe) -(q_1(f, t_I, alpha, beta, 
                            i, j, especie_1, classe)*q_1(f, t_K, alpha, beta, 
                                                   k, l, especie_3, classe)))
            cov = r_1*soma
            return cov
        else:
            for x in comb_transi:
                i, j = x[0], x[1]
                for y in comb_transv:
                    k, l = y[0], y[1]
                    soma +=  (q_2(t_I, t_J, t_K, t_L, f, alpha, beta, 
                                 i, j, k, l, especie_1, especie_2, 
                                 especie_3, especie_4, classe) -(q_1(f, t_I, alpha, beta, 
                            i, j, especie_1, classe)*q_1(f, t_K, alpha, beta, 
                                                   k, l, especie_3, classe)))
            cov = r_2*soma
            return cov
    elif tipo1 == 0 and tipo2 == 1:
        soma = 0
        if classe == 1:
            for x in comb_transv:
                i, j = x[0], x[1]
                for y in comb_transi:
                    k, l = y[0], y[1]
                    soma +=  (q_2(t_I, t_J, t_K, t_L, f, alpha, beta, 
                                 i, j, k, l, especie_1, especie_2, 
                                 especie_3, especie_4, classe) -(q_1(f, t_I, alpha, beta, 
                            i, j, especie_1, classe)*q_1(f, t_K, alpha, beta, 
                                                   k, l, especie_3, classe)))
            cov = r_1*soma
            return cov
        else:
            for x in comb_transv:
                i, j = x[0], x[1]
                for y in comb_transi:
                    k, l = y[0], y[1]
                    soma +=  (q_2(t_I, t_J, t_K, t_L, f, alpha, beta, 
                                 i, j, k, l, especie_1, especie_2, 
                                 especie_3, especie_4, classe) -(q_1(f, t_I, alpha, beta, 
                            i, j, especie_1, classe)*q_1(f, t_K, alpha, beta, 
                                                   k, l, especie_3, classe)))
            cov = r_2*soma
            return cov

# vetor par =  [t1, t3, t4, t5, t6, f1, alpha1, beta1, f2, alpha2, beta2]       
# V_var (f, beta, especie, classe, t_esp)
# S_var (f, beta, alpha, especie, classe, t_esp)
# se tipo == 0 --> V_k^(i)
# se tipo == 1 --> S_k^(i)
            
def cov_k (par, especie_1, especie_2, classe1, classe2, tipo1, tipo2):
    if classe1 == classe2:
        if classe1 == 1:
            f, alpha, beta = par[5], par[6], par[7]
        elif classe1 == 2:
            f, alpha, beta = par[8], par[9], par[10]
        if especie_1 == especie_2:
            soma = 0
            coef = (1/(num_especie - especie_1)**2)
            i , k = especie_1, especie_1
            if especie_1 == 1:
                t_I, t_K = par[especie_1 - 1], par[especie_1 - 1]
            elif especie_1 == 2:
                t_I, t_K = t_2, t_2
            elif especie_1 > 2:
                t_I, t_K = par[especie_1 - 2], par[especie_1 - 2]
            for j in range(i, num_especie):
                for l in range(j + 1, num_especie):
                    if j + 1 == 2:
                        t_J = t_2
                    elif j + 1 > 2 and j + 1 < 7:
                        t_J = par[(j + 1) - 2]
                    elif j + 1 == 7:
                        t_J = None
                    if l + 1 == 2:
                        t_L = t_2
                    elif l + 1 > 2 and l + 1 < 7:
                        t_L = par[(l + 1) - 2]
                    elif l + 1 == 7:
                        t_L = None
                    if classe1 == 1:
                        f, alpha, beta = par[5], par[6], par[7]
                    elif classe1 == 2:
                        f, alpha, beta = par[8], par[9], par[10]
                    soma += (2*cov(t_I, t_J, t_K, t_L, i, j + 1, k, l + 1, classe1, tipo1, tipo2,
                                f, alpha, beta))
            if tipo1 == tipo2 and tipo1 == 0:
                total = coef*(((7-i)*(V_var(f, beta, especie_1, classe1, t_I))) + soma)
                return total
            elif tipo1 == tipo2 and tipo1 == 1 :
                total = coef*(((7-i)*(S_var(f, beta, alpha, especie_1, classe1, t_I))) + soma)
                return total
            elif tipo1 != tipo2:
                total = coef*(((7-i)*cov_V_S(V_barra(f, beta, especie_1, classe1, t_I),
                                S_barra(f, beta, alpha, especie_1, classe1, t_I), classe1)) 
                        + soma)
                return total
        elif especie_1 != especie_2:
            soma = 0
            coef = (1/(num_especie - especie_1))*(1/(num_especie - especie_2))
            i, k = especie_1, especie_2
            if especie_1 == 1:
                t_I = par[especie_1 - 1]
            elif especie_1 == 2:
                t_I = t_2
            elif especie_1 > 2:
                t_I = par[especie_1 - 2]
            if especie_2 == 2:
                t_K = t_2
            elif especie_2 > 2 and especie_2 < 7:
                t_K = par[especie_2 - 2]
            for j in range(i, num_especie):
                for l in range(k, num_especie):
                    if j + 1 == 2:
                        t_J = t_2
                    elif j + 1 > 2 and j + 1 < 7:
                        t_J = par[(j + 1) - 2]
                    elif j + 1 == 7:
                        t_J = None
                    if l + 1 == 2:
                        t_L = t_2
                    elif l + 1 > 2 and l + 1 < 7:
                        t_L = par[(l + 1) - 2]
                    elif l + 1 == 7:
                        t_L = None
                    soma += cov(t_I, t_J, t_K, t_L, i, j + 1, k, l + 1, classe1, tipo1, tipo2,
                                f, alpha, beta)
            total = coef*soma
            return total
    elif classe1 != classe2:
        return 0
    
# Testando 
par_teste = np.array([80, 50, 30, 20, 10, 1, 0.6, 0.07, 0.5, 0.5, 0.4])
cov_k(par_teste, 2, 5, 2, 2, 1, 0)

a, b = 1, 1

# teste do minimize
def f(par):
    fun = (par[0] - par[1])**2 + (par[1])**2
    return fun

x0 = np.array([0.5,0.2])
res = minimize(f, x0)
print(res.x)

# cov_k (par, especie_1, especie_2, classe1, classe2, tipo1, tipo2)
# Definindo R
def R(par):
    nvar = 24
    omega = np.zeros((nvar,nvar))
    # Montando a matriz omega
    for i in range(0, nvar):
        for j in range(i, nvar):
            if 0 <= i <= 5 and 0 <= j <= 5:
                omega[i,j] = cov_k(par, i + 1, j + 1, 1, 1, 0, 0)
                if i != j:
                    omega[j, i] = omega[i, j]
            elif 0 <= i <= 5 and 6 <= j <= 11 :
                if (j - 5) <  (i + 1):
                    omega[i,j] = cov_k(par, j - 5, i + 1, 1, 1, 1, 0)
                    if i != j:
                        omega[j, i] = omega[i, j]
                elif (j - 5) >= (i + 1):
                    omega[i, j] = cov_k(par, i + 1, j - 5, 1, 1, 0, 1)
                    if i != j:
                        omega[j, i] = omega[i, j]
            elif 0 <= i <= 5 and j >= 12:
                omega[i, j] = 0
                omega[j, i] = 0
            elif 6 <= i <= 11 and 6 <= j <= 11:
                omega[i ,j] = cov_k(par, i - 5, j - 5, 1, 1, 1, 1)
                if i != j:
                    omega[j, i] = omega[i, j]
            elif 6 <= i <= 11 and j >= 12:
                omega[i, j] = 0
                omega[j, i] = 0
            elif 12 <= i <= 17 and 12 <= j <= 17:
                omega[i, j] = cov_k(par, i - 11, j - 11, 2, 2, 0, 0)
                if i != j:
                    omega[j, i] = omega[i, j]
            elif 12 <= i <= 17 and 18 <= j <= 23:
                if (j - 17) < (i - 11):
                    omega[i, j] = cov_k(par, j - 17, i - 11,2, 2, 1, 0)
                    if i != j:
                        omega[j, i] = omega[i, j]
                elif (j - 17) >= (i - 11):
                    omega[i, j] = cov_k(par, i - 11, j - 17, 2, 2, 0, 1)
                    if i != j:
                        omega[j, i] = omega[i, j]
            elif 18 <= i <= 23 and 18 <= j <= 23:
                omega[i, j] = cov_k(par, i - 17, j - 17, 2, 2, 1, 1)
                if i != j:
                    omega[j, i]  = omega[i, j]
    # invertendo omega
    omega
    omega_inv = np.linalg.inv(omega)
    print(np.dot(omega, omega_inv))
    # definindo o vetor de médias
    return omega_inv

inicio = time.time()
R(par_teste)
fim = time.time()
tempo = fim - inicio

