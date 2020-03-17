import numpy as np
import math
import sympy as sm
import time

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
pi_1 = np.array([0.169, 0.429, 0.364, 0.038])
pi_2 = np.array([0.297, 0.267, 0.310, 0.126])
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
pi_1 = sm.Matrix([[0.169, 0.429, 0.364, 0.038]])
pi_2 = sm.Matrix([[0.297, 0.267, 0.310, 0.126]])
r_1 = 232
r_2 = 667
t_2 = 65
t_1, t_3, t_4, t_5, t_6, f_1, alpha_1, beta_1, beta_2, f_2, alpha_2, beta_2 = sm.symbols(
        r't_1 t_3 t_4 t_5 t_6 f_1 \alpha_1 \beta_1 \beta_2 f_2 \alpha_2 \beta_2')

# Matriz de transição
def P_trans(t, classe):
    if classe == 1:
        beta, alpha = sm.symbols('beta_1, alpha_1')
        pi_y = pi_1[0] + pi_1[1]
        pi_r = pi_1[2] + pi_1[3]
        part_1 = sm.Matrix([[1,1,1,1]]).transpose()*sm.Matrix([[pi_1
                          [0], pi_1[1], pi_1[2], pi_1[3]]])
        part_2 = sm.exp(-(beta* t))*sm.Matrix([[1/pi_y, 
                       1/pi_y, -1/pi_r, -1/pi_r]]).transpose()*sm.Matrix([[pi_r*pi_1[0],
                                                pi_r*pi_1[1], -pi_y*pi_1[2], -pi_y*pi_1[3]]])
        part_3 = sm.exp(-t*(pi_y*beta + pi_r*alpha))*sm.Matrix([[0, 0, 
                       pi_1[3]/pi_r, -pi_1[2]/pi_r]]).transpose()*sm.Matrix([[0, 0, 1, -1]])
        part_4 = sm.exp(-t*(pi_y*alpha + pi_r*beta))*sm.Matrix([[pi_1[1]/pi_y,
                       -pi_1[0]/pi_y, 0, 0]]).transpose()*sm.Matrix([[1, -1, 0, 0]])
        mat_trans = part_1 + part_2 + part_3 + part_4
        return mat_trans
    else:
        beta, alpha = sm.symbols('beta_2, alpha_2')
        pi_y = pi_2[0] + pi_2[1]
        pi_r = pi_2[2] + pi_2[3]
        part_1 = sm.Matrix([[1,1,1,1]]).transpose()*sm.Matrix([[pi_2
                          [0], pi_2[1], pi_2[2], pi_2[3]]])
        part_2 = sm.exp(-(beta* t))*sm.Matrix([[1/pi_y, 
                       1/pi_y, -1/pi_r, -1/pi_r]]).transpose()*sm.Matrix([[pi_r*pi_2[0],
                                                pi_r*pi_2[1], -pi_y*pi_2[2], -pi_y*pi_2[3]]])
        part_3 = sm.exp(-t*(pi_y*beta + pi_r*alpha))*sm.Matrix([[0, 0, 
                       pi_2[3]/pi_r, -pi_2[2]/pi_r]]).transpose()*sm.Matrix([[0, 0, 1, -1]])
        part_4 = sm.exp(-t*(pi_y*alpha + pi_r*beta))*sm.Matrix([[pi_2[1]/pi_y,
                       -pi_2[0]/pi_y, 0, 0]]).transpose()*sm.Matrix([[1, -1, 0, 0]])
        mat_trans = part_1 + part_2 + part_3 + part_4
        return mat_trans

# Testando
P_trans(t_1, 1)[0,0]
P_trans(t_1, 2)[0,0]

# Medias e variancias
def V_barra (classe, especie):
    if classe == 1:
        pi_y = pi_1[0] + pi_1[1]
        pi_r = pi_1[2] + pi_1[3]
        f, beta = sm.symbols(r'f_1 \beta_1')
        if especie == 2:
            media = 2*f*r_1*(pi_y*pi_r)*(1 - sm.exp(-2*beta*t_2))
        else:
            t_esp = sm.Symbol('t_{}'.format(especie))
            media = 2*f*r_1*(pi_y*pi_r)*(1 - sm.exp(-2*beta*t_esp))
        return media
    
    else:
        pi_y = pi_2[0] + pi_2[1]
        pi_r = pi_2[2] + pi_2[3]
        f, beta = sm.symbols(r'f_1 \beta_1')
        if especie == 2:
            media = 2*f*r_2*(pi_y*pi_r)*(1 - sm.exp(-2*beta*t_2))
        else:
            t_esp = sm.Symbol('t_{}'.format(especie))
            media = 2*f*r_2*(pi_y*pi_r)*(1 - sm.exp(-2*beta*t_esp))
        return media

def S_barra (classe, especie):
    if classe == 1:
        pi_y = pi_1[0] + pi_1[1]
        pi_r = pi_1[2] + pi_1[3]
        f, beta, alpha = sm.symbols(r'f_1 \beta_1 \alpha_1')
        if especie == 2:
            media = 2*f*r_1*(((pi_1[0]*pi_1[1])+(pi_1[2]*pi_1[3])) + 
                       (((pi_1[0]*pi_1[1]*pi_r)/pi_y + (pi_1[2]*pi_1[3]*pi_y)/pi_r)*
                       sm.exp(-2*beta*t_2)) -
                       ((pi_1[0]*pi_1[1]/pi_y) * sm.exp(-(2*t_2)*(alpha*pi_y + beta*pi_r))) -
                       ((pi_1[2]*pi_1[3]/pi_r) * sm.exp(-(2*t_2)*(alpha*pi_r + beta*pi_y))))
        else:
            t = sm.Symbol('t_{}'.format(especie))
            media = 2*f*r_1*(((pi_1[0]*pi_1[1])+(pi_1[2]*pi_1[3])) + 
                       (((pi_1[0]*pi_1[1]*pi_r)/pi_y + (pi_1[2]*pi_1[3]*pi_y)/pi_r)*
                       sm.exp(-2*beta*t)) -
                       ((pi_1[0]*pi_1[1]/pi_y) * sm.exp(-(2*t)*(alpha*pi_y + beta*pi_r))) -
                       ((pi_1[2]*pi_1[3]/pi_r) * sm.exp(-(2*t)*(alpha*pi_r + beta*pi_y))))
        return media
    else:
        pi_y = pi_2[0] + pi_2[1]
        pi_r = pi_2[2] + pi_2[3]
        f, beta, alpha = sm.symbols(r'f_2 \beta_2 \alpha_2')
        if especie == 2:
            media = 2*f*r_2*(((pi_2[0]*pi_2[1])+(pi_2[2]*pi_2[3])) + 
                       (((pi_2[0]*pi_2[1]*pi_r)/pi_y + (pi_2[2]*pi_2[3]*pi_y)/pi_r)*
                       sm.exp(-2*beta*t_2)) -
                       ((pi_2[0]*pi_2[1]/pi_y) * sm.exp(-(2*t_2)*(alpha*pi_y + beta*pi_r))) -
                       ((pi_2[2]*pi_2[3]/pi_r) * sm.exp(-(2*t_2)*(alpha*pi_r + beta*pi_y))))
        else:
            t = sm.Symbol('t_{}'.format(especie))
            media = 2*f*r_2*(((pi_2[0]*pi_2[1])+(pi_2[2]*pi_2[3])) + 
                       (((pi_2[0]*pi_2[1]*pi_r)/pi_y + (pi_2[2]*pi_2[3]*pi_y)/pi_r)*
                       sm.exp(-2*beta*t)) -
                       ((pi_2[0]*pi_2[1]/pi_y) * sm.exp(-(2*t)*(alpha*pi_y + beta*pi_r))) -
                       ((pi_2[2]*pi_2[3]/pi_r) * sm.exp(-(2*t)*(alpha*pi_r + beta*pi_y))))
        return media

def V_var (classe, especie):
    if classe == 1:
        var = V_barra(classe, especie)*(1 - (V_barra(classe, especie)/r_1))
        return var
    else:
        var = V_barra(classe, especie)*(1 - (V_barra(classe, especie)/r_2))
        return var

def S_var (classe, especie):
    if classe == 1:
        var = S_barra(classe, especie)*(1 - (S_barra(classe, especie)/r_1))
        return var
    else:
        var = S_barra(classe, especie)*(1 - (S_barra(classe, especie)/r_2))
        return var

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

# Definindo o vetor d_barra de medias
D_barra = sm.Matrix([[V_barra(1,1), V_barra(1,2),V_barra(1,3), V_barra(1,4), V_barra(1,5),
V_barra(1,6), S_barra(1,1), S_barra(1,2), S_barra(1,3), S_barra(1,4), S_barra(1,5),
S_barra(1,6), V_barra(2,1), V_barra(2,2), V_barra(2,3), V_barra(2,4), V_barra(2,5), 
V_barra(2,6), S_barra(2,1), S_barra(2,2), S_barra(2,3), S_barra(2,4), S_barra(2,5),
S_barra(2,6)]]).transpose()

# Definindo o vetor de observações
D_tilde = sm.Matrix([[V_k(1,1), V_k(1,2),V_k(1,3), V_k(1,4), V_k(1,5),
V_k(1,6), S_k(1,1), S_k(1,2), S_k(1,3), S_k(1,4), S_k(1,5),
S_k(1,6), V_k(2,1), V_k(2,2), V_k(2,3), V_k(2,4), V_k(2,5), 
V_k(2,6), S_k(2,1), S_k(2,2), S_k(2,3), S_k(2,4), S_k(2,5),
S_k(2,6)]]).transpose()

    
# Definindo a covarincia e variancia de um V_k^(i) e S_k(i)
# se tipo == 0 --> V_k^(i)
# se tipo == 1 --> S_k^(i)
def q_1(i,j, especie_1, especie_2,classe):
    if classe == 1:
        f = sm.Symbol('f_1')
        if especie_1 != 2:
            t = sm.Symbol('t_{}'.format(especie_1))
            q = f*pi_1[i]*P_trans(t, classe)[i, j]
            return q
        else:
            q = f*pi_1[i]*P_trans(t_2, classe)[i, j]
            return q
    else:
        f = sm.Symbol('f_2')
        if especie_1 != 2:
            t = sm.Symbol('t_{}'.format(especie_1))
            q = f*pi_2[i]*P_trans(t, classe)[i, j]
            return q
        else:
            q = f*pi_2[i]*P_trans(t_2, classe)[i, j]
            return q
        
def q_2(i,j,k,l, especie_1, especie_2, especie_3, especie_4, classe):
    if classe == 1:
        f = sm.Symbol('f_1')
        if especie_1 < especie_2 < especie_3 < especie_4:
            soma = 0
            t_I = sm.Symbol('t_{}'.format(especie_1))
            t_J = sm.Symbol('t_{}'.format(especie_2))
            t_K = sm.Symbol('t_{}'.format(especie_3))
            for x in range(0, pi_1.shape[1]):
                for y in range(0, pi_1.shape[1]):
                    soma += pi_1[x]*P_trans(2*t_I - t_J, classe)[x, i]*P_trans(t_J,
                                classe)[x, j]*P_trans(t_J - t_K, classe)[x, y]*P_trans(t_K,
                                       classe)[y, k]*P_trans(t_K, classe)[y, l]
            q = f*soma
            return q
        elif especie_2 > especie_3 and especie_2 < especie_4:
            soma = 0
            t_I = sm.Symbol('t_{}'.format(especie_1))
            t_J = sm.Symbol('t_{}'.format(especie_2))
            t_K = sm.Symbol('t_{}'.format(especie_3))
            for x in range(0, pi_1.shape[1]):
                for y in range(0, pi_1.shape[1]):
                    soma += pi_1[x]*P_trans(2*t_I - t_K, classe)[x, i]*P_trans(t_K,
                                classe)[x, k]*P_trans(t_K - t_J, classe)[x, y]*P_trans(t_J,
                                       classe)[y, j]*P_trans(t_J, classe)[y, l]
            q = f*soma
            return q
        elif especie_2 == especie_3:
            soma = 0
            t_I = sm.Symbol('t_{}'.format(especie_1))
            t_J = sm.Symbol('t_{}'.format(especie_2))
            if j == k:
                for x in range(0, pi_1.shape[1]):
                    soma += pi_1[x]*P_trans(2*t_I - t_J, classe)[x, i]*P_trans(t_J, 
                                classe)[x, j]*P_trans(t_J, classe)[x, l]
            q = f*soma
            return q
        elif especie_2 == especie_4:
            soma = 0
            t_I = sm.Symbol('t_{}'.format(especie_1))
            t_J = sm.Symbol('t_{}'.format(especie_2))
            t_K = sm.Symbol('t_{}'.format(especie_3))
            if j == l:
                for x in range(0, pi_1.shape[1]):
                    soma += pi_1[x]*P_trans(2*t_I - t_K, classe)[x, i]*P_trans(t_K,
                                classe)[x, k]*P_trans(t_K, classe)[x, l]
            q = f*soma
            return q
        elif especie_2 > especie_4:
            soma =  0
            t_I = sm.Symbol('t_{}'.format(especie_1))
            t_J = sm.Symbol('t_{}'.format(especie_2))
            t_K = sm.Symbol('t_{}'.format(especie_3))
            t_L = sm.Symbol('t_{}'.format(especie_4))
            for x in range(0, pi_1.shape[1]):
                for y in range(0, pi_1.shape[1]):
                    soma += pi_1[x]*P_trans(2*t_I - t_K, classe)[x, i]*P_trans(t_L,
                                classe)[x, l]*P_trans(t_K - t_L, classe)[x, y]*P_trans(t_L,
                                       classe)[y, l]*P_trans(t_L, classe)[y, j]
            q = f*soma
            return q
    else:
        f = sm.Symbol('f_2')
        if especie_1 < especie_2 < especie_3 < especie_4:
            soma = 0
            t_I = sm.Symbol('t_{}'.format(especie_1))
            t_J = sm.Symbol('t_{}'.format(especie_2))
            t_K = sm.Symbol('t_{}'.format(especie_3))
            for x in range(0, pi_2.shape[1]):
                for y in range(0, pi_1.shape[1]):
                    soma += pi_2[x]*P_trans(2*t_I - t_J, classe)[x, i]*P_trans(t_J,
                                classe)[x, j]*P_trans(t_J - t_K, classe)[x, y]*P_trans(t_K,
                                       classe)[y, k]*P_trans(t_K, classe)[y, l]
            q = f*soma
            return q
        elif especie_2 > especie_3 and especie_2 < especie_4:
            soma = 0
            t_I = sm.Symbol('t_{}'.format(especie_1))
            t_J = sm.Symbol('t_{}'.format(especie_2))
            t_K = sm.Symbol('t_{}'.format(especie_3))
            for x in range(0, pi_2.shape[1]):
                for y in range(0, pi_1.shape[1]):
                    soma += pi_2[x]*P_trans(2*t_I - t_K, classe)[x, i]*P_trans(t_K,
                                classe)[x, k]*P_trans(t_K - t_J, classe)[x, y]*P_trans(t_J,
                                       classe)[y, j]*P_trans(t_J, classe)[y, l]
            q = f*soma
            return q
        elif especie_2 == especie_3:
            soma = 0
            t_I = sm.Symbol('t_{}'.format(especie_1))
            t_J = sm.Symbol('t_{}'.format(especie_2))
            if j == k:
                for x in range(0, pi_2.shape[1]):
                    soma += pi_2[x]*P_trans(2*t_I - t_J, classe)[x, i]*P_trans(t_J, 
                                classe)[x, j]*P_trans(t_J, classe)[x, l]
            q = f*soma
            return q
        elif especie_2 == especie_4:
            soma = 0
            t_I = sm.Symbol('t_{}'.format(especie_1))
            t_J = sm.Symbol('t_{}'.format(especie_2))
            t_K = sm.Symbol('t_{}'.format(especie_3))
            if j == l:
                for x in range(0, pi_2.shape[1]):
                    soma += pi_2[x]*P_trans(2*t_I - t_K, classe)[x, i]*P_trans(t_K,
                                classe)[x, k]*P_trans(t_K, classe)[x, l]
            q = f*soma
            return q
        elif especie_2 > especie_4:
            soma =  0
            t_I = sm.Symbol('t_{}'.format(especie_1))
            t_J = sm.Symbol('t_{}'.format(especie_2))
            t_K = sm.Symbol('t_{}'.format(especie_3))
            t_L = sm.Symbol('t_{}'.format(especie_4))
            for x in range(0, pi_2.shape[1]):
                for y in range(0, pi_2.shape[1]):
                    soma += pi_2[x]*P_trans(2*t_I - t_K, classe)[x, i]*P_trans(t_L,
                                classe)[x, l]*P_trans(t_K - t_L, classe)[x, y]*P_trans(t_L,
                                       classe)[y, l]*P_trans(t_L, classe)[y, j]
            q = f*soma
            return q
        
def cov(especie_1, especie_2, especie_3, especie_4, classe, tipo1, tipo2):
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
                    soma += (q_2(i, j, k, l, especie_1, especie_2, especie_3, especie_4, classe) -(
q_1(i, j, especie_1, especie_2, classe)*q_1(k, l, especie_3, especie_4, classe)))
            cov = r_1*soma
            return cov
        else:
            for x in comb_transv:
                i, j = x[0], x[1]
                for y in comb_transv:
                    k, l = y[0], y[1]
                    soma += (q_2(i, j, k, l, especie_1, especie_2, especie_3, especie_4, classe) -(
                            q_1(i, j, especie_1, especie_2, classe)*q_1(k, l, especie_3,
                               especie_4, classe)))
            cov = r_2*soma
            return cov
    elif tipo1 == tipo2 == 1:
        soma = 0
        if classe == 1:
            for x in comb_transi:
                i, j = x[0], x[1]
                for y in comb_transi:
                    k, l = y[0], y[1]
                    soma += (q_2(i, j, k, l, especie_1, especie_2, especie_3, especie_4, classe) -(
                            q_1(i, j, especie_1, especie_2, classe)*q_1(k, l, especie_3, 
                               especie_4, classe)))
            cov = r_1*soma
            return cov
        else:
            for x in comb_transi:
                i, j = x[0], x[1]
                for y in comb_transi:
                    k, l = y[0], y[1]
                    soma += (q_2(i, j, k, l, especie_1, especie_2, especie_3, especie_4, classe) -(
                            q_1(i, j, especie_1, especie_2, classe)*q_1(k, l, 
                               especie_3, especie_4, classe)))
            cov = r_2*soma
            return cov
    elif tipo1 == 1 and tipo2 == 0:
        soma = 0
        if classe == 1:
            for x in comb_transi:
                i, j = x[0], x[1]
                for y in comb_transv:
                    k, l = y[0], y[1]
                    soma += (q_2(i, j, k, l, especie_1, especie_2, especie_3, especie_4, classe) -(
                            q_1(i, j, especie_1, especie_2, classe)*q_1(k, l,
                               especie_3, especie_4, classe)))
            cov = r_1*soma
            return cov
        else:
            for x in comb_transi:
                i, j = x[0], x[1]
                for y in comb_transv:
                    k, l = y[0], y[1]
                    soma += (q_2(i, j, k, l, especie_1, especie_2, especie_3, especie_4, classe) -(
                            q_1(i, j, especie_1, especie_2, classe)*q_1(k, l,
                               especie_3, especie_4, classe)))
            cov = r_2*soma
            return cov
    elif tipo1 == 0 and tipo2 == 1:
        soma = 0
        if classe == 1:
            for x in comb_transv:
                i, j = x[0], x[1]
                for y in comb_transi:
                    k, l = y[0], y[1]
                    soma += (q_2(i, j, k, l, especie_1, especie_2, especie_3, especie_4 ,classe) -(
                            q_1(i, j, especie_1, especie_2, classe)*q_1(k, l, 
                               especie_3, especie_4, classe)))
            cov = r_1*soma
            return cov
        else:
            for x in comb_transv:
                i, j = x[0], x[1]
                for y in comb_transi:
                    k, l = y[0], y[1]
                    soma += (q_2(i, j, k, l, especie_1, especie_2, especie_3, especie_4, classe) -(
                            q_1(i, j, especie_1, especie_2, classe)*q_1(k, l, 
                               especie_3, especie_4, classe)))
            cov = r_2*soma
            return cov
        
            
def cov_k (especie_1, especie_2, classe1, classe2, tipo1, tipo2):
    if classe1 == classe2:
        if especie_1 == especie_2:
            soma = 0
            coef = (1/(num_especie - especie_1)**2)
            i , k = especie_1, especie_1
            for j in range(i,num_especie):
                for l in range(j + 1, num_especie):
                    soma += cov(i,j + 1, k,l + 1, classe1, tipo1, tipo2)
            total = coef*(((6-i)*(V_var(classe1, especie_1))) + soma)
            return total
        elif especie_1 != especie_2:
            soma = 0
            coef = (1/(num_especie - especie_1))*(1/(num_especie - especie_2))
            i, k = especie_1, especie_2
            for j in range(i, num_especie):
                for l in range(k, num_especie):
                    soma += cov(i, j+1, k, l+1, classe1, tipo1, tipo2)
            total = coef*soma
            return total
    elif classe1 != classe2:
        return 0

# testando as covariancias 
# Para V_1^(1) e V_1^(2)
inicio = time.time()
h = cov_k(1, 2, classe1 = 1, classe2 = 1, tipo1 = 0, tipo2 = 0)
final  = time.time()
final - inicio
# Demora para o  calculo
# Apesar da demora, testa-se a elaboração da matriz omega
inicio = time.time()
omega = sm.zeros(24)
for i in range(0,omega.shape[1]):
    for j in range(i, omega.shape[1]):
        if j <= 5:
            omega[i,j] = cov_k(i+1, j+1, classe1 = 1, classe2 = 1, tipo1 = 0, tipo2 = 0)
            if i != j:
                omega[j, i] = omega[i, j]
        elif 6 <= j <= 11 and i <= 5:
            omega[i,j] = cov_k(i+1, j - 5, classe1 = 1, classe2 = 1, tipo1 = 0, tipo2 = 1)
            omega[j, i] = omega[i, j]
        elif 12 <= j <= 23 and i <= 5:
            omega[i, j] = 0
            omega[j, i] = omega[i, j]
        elif 6 <= j <= 11 and 6 <= i <= 11:
            omega[i, j] = cov_k(i - 5, j - 5, classe1 = 1, classe2 = 1, tipo1 = 1, tipo2 = 1)
            if i!= j:
                omega[j, i] = omega[i, j]
        elif 12 <= j <= 23 and 6 <= i <= 11:
            omega[i, j]= 0
            omega[j, i] = omega[i, j]
        elif 12 <= j <= 17 and 12 <= i <= 17:
            omega[i, j] = cov_k(i - 11, j - 11, classe1 = 2, classe2 = 2, tipo1 = 0, tipo2 = 0)
            if i != j:
                omega[j, i] = omega[i, j]
        elif 18 <= j <= 23 and 12 <= i <= 17:
            omega[i, j] = cov_k(i - 11, j - 17, classe1 = 2, classe2 = 2, tipo1 = 0, tipo2 = 1)
            omega[j, i] = omega[i, j]
        elif 18 <= j <= 23 and 18 <= i <= 23:
            omega[i, j] = cov_k(i - 17, j - 17, classe1 = 2, classe2 = 2, tipo1 = 1, tipo2 = 1)
            if i != j:
                omega[j, i] = omega[i, j]
final = time.time()
tempo = final - inicio
final - inicio
determin = omega.det()
    
