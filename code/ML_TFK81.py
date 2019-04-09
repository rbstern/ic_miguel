class Especie:
  def __init__(self, meus_filhos, meu_valor):
    self.filhos = meus_filhos
    self.valor = meu_valor

S1 = Especie(None, 1)
S2 = Especie(None, 0)
S3 = Especie(None, 0)
S4 = Especie(None, 0)
S5 = Especie(None, 1)
S6 = Especie([S1, S2], None)
S7 = Especie([S4, S5], None)
S8 = Especie([S3, S7], None)
S0 = Especie([S6, S8], None)

folhas = [S1, S2, S3, S4, S5]

def transicao(a, b):
  if(a == b):
    return 0.75
  return 0.25

def L_condicional(esta_especie, valor_ancestral):
  if(not(esta_especie.filhos)):
    return transicao(esta_especie.valor, valor_ancestral)
  else:
    return (transicao(0, valor_ancestral)*
L_condicional(esta_especie.filhos[0], 0)*
L_condicional(esta_especie.filhos[1], 0))+(transicao(1, valor_ancestral)*
L_condicional(esta_especie.filhos[0], 1)*
L_condicional(esta_especie.filhos[1], 1))

testeS0=L_condicional(S0,1)
testeS8=L_condicional(S8,1)
testeS6=L_condicional(S6,1)
testeS7=L_condicional(S7,1)
testeS3=L_condicional(S3,1)
print(testeS3)

