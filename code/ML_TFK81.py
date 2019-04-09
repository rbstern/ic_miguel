class Especie:
    def __init__(self, meu_pai, meu_valor):
        self.filhos = []
        self.pai = meu_pai
        self.valor = meu_valor
        if(meu_pai):
            meu_pai.filhos.append(self)

def priori(a):
    return 0.5

def transicao(a, b):
    if(a == b):
        return 0.75
    return 0.25

def L_condicional(esta_especie, valor_ancestral):
    if(not(esta_especie.filhos)):
        return transicao(esta_especie.valor, valor_ancestral)
    elif(not(esta_especie.pai)):
        return(priori(0)*
               L_condicional(esta_especie.filhos[0], 0)*
                L_condicional(esta_especie.filhos[1], 0))+(priori(1)*
                             L_condicional(esta_especie.filhos[0], 1)*
                             L_condicional(esta_especie.filhos[1], 1))
    else: 
        return (transicao(0, valor_ancestral)*
                L_condicional(esta_especie.filhos[0], 0)*
                L_condicional(esta_especie.filhos[1], 0))+(transicao(1, valor_ancestral)*
                             L_condicional(esta_especie.filhos[0], 1)*
                             L_condicional(esta_especie.filhos[1], 1))

S0 = Especie(None, None)
S6 = Especie(S0, None)
S1 = Especie(S6, 1)
S2 = Especie(S6, 0)
S8 = Especie(S0, None)
S3 = Especie(S8, 0)
S7 = Especie(S8, None)
S4 = Especie(S7, 0)
S5 = Especie(S7, 1)
likelihood = L_condicional(S0, None)
print(likelihood)

S0 = Especie(None, None)
S6 = Especie(S0, None)
S1 = Especie(S6, 1)
S2 = Especie(S6, 0)
S8 = Especie(S0, None)
S5 = Especie(S8, 1)
S7 = Especie(S8, None)
S3 = Especie(S7, 0)
S4 = Especie(S7, 0)
likelihood = L_condicional(S0, None)
print(likelihood)

# Achar a arvore de maxima verossimilhança e
# ler The base substitution probabilities
# Bonus: adaptar o que fizemos para o caso em que
# valor pode ser "A", "T", "C" ou "G".
