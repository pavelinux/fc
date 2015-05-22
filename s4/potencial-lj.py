from numpy import linspace

# Parametros de LJ
# Caso de Ar
E_0 = 0.997 #(kJ/mol)
S = 3.41 #(A)
N = 1000
r = linspace(0.00001, 10.0, N)
U = 4 * E_0 * ((S/r)**12 - (S/r)**6)
for i in range(N):
    print r[i], " ", U[i]
