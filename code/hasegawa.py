# Matrizes do numero de transições e transversoes observadas entre especies
# 0 - RATO 
# 1 - BOI
# 2 - GIBAO
# 3 - ORANGOTANGO
# 4 - GORILA
# 5 - CHIMPANZE
# 6 - HOMEM
import numpy as np
import math

# Nucleotideos de classe 1 --> 3ª posiçao no codon
# Nucleotideos de classe 2 --> 1ª e 2ª posiçoes no codon

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
# T = 0, C = 1, A = 3, G = 4
pi_1 = np.array([0.169, 0.429, 0.364, 0.038])
pi_2 = np.array([0.297, 0.267, 0.310, 0.126])

# Medias e variancias
def V_barra (f, r, beta, t, classe):
    if classe == 1:
        pi_y = pi_1[0] + pi_1[1]
        pi_r = pi_1[2] + pi_1[3]
        media = 2*f*r*(pi_y*pi_r)*(1 - math.exp(-2*beta*t))
        return media
    else:
        pi_y = pi_2[0] + pi_2[1]
        pi_r = pi_2[2] + pi_2[3]
        media = 2*f*r*(pi_y*pi_r)*(1 - math.exp(-2*beta*t))
        return media

def S_barra (f, r, beta, alpha, t, classe):
    if classe == 1:
        pi_y = pi_1[0] + pi_1[1]
        pi_r = pi_1[2] + pi_1[3]
        media = 2*f*r*(((pi_1[0]*pi_1[1])+(pi_1[2]*pi_1[3])) + 
                       (((pi_1[0]*pi_1[1]*pi_r)/pi_y + (pi_1[2]*pi_1[3]*pi_y)/pi_r)*
                       math.exp(-2*beta*t)) -
                       ((pi_1[0]*pi_1[1]/pi_y) * math.exp(-2*t*(alpha*pi_y + beta*pi_r))) -
                       ((pi_1[2]*pi_1[3]/pi_r) * math.exp(-2*t*(alpha*pi_r + beta*pi_y))))
        return media
    else:
        pi_y = pi_2[0] + pi_2[1]
        pi_r = pi_2[2] + pi_2[3]
        media = 2*f*r*(((pi_2[0]*pi_2[1])+(pi_2[2]*pi_2[3])) + 
                       (((pi_2[0]*pi_2[1]*pi_r)/pi_y + (pi_2[2]*pi_2[3]*pi_y)/pi_r)*
                       math.exp(-2*beta*t)) -
                       ((pi_2[0]*pi_2[1]/pi_y) * math.exp(-2*t*(alpha*pi_y + beta*pi_r))) -
                       ((pi_2[2]*pi_2[3]/pi_r) * math.exp(-2*t*(alpha*pi_r + beta*pi_y))))
        return media

def V_var (f, r, beta, t, classe):
    if classe == 1:
        var = V_barra(f, r, beta, t, 1)*(1 - (V_barra(f, r, beta, t, 1)/r))
        return var
    else:
        var = V_barra(f, r, beta, t, 2)*(1 - (V_barra(f, r, beta, t, 2)/r))
        return var

# Definindo transversoes e transiçoes amostrais
def V_k(especie, classe):
    soma = 0
    if classe == 1:
        j = especie + 1
        while j <= 6:
            soma += V1[especie,j]
            j += 1
        v_k = soma/(7 - (especie + 1))
        return v_k
    else:
        j = especie + 1
        while j <= 6:
            soma += V2[especie, j]
            j += 1
        v_k = soma/(7 - (especie + 1))
        return v_k

def S_k(especie, classe):
    soma = 0
    if classe == 1:
        j = especie + 1
        while j <= 6:
            soma += S1[especie,j]
            j += 1
        s_k = soma/(7 - (especie + 1))
        return s_k
    else:
        j = especie + 1
        while j <= 6:
            soma += S2[especie, j]
            j += 1
        s_k = soma/(7 - (especie + 1))
        return s_k

    
