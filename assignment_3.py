import copy
from math import *

# Funcoes utilizadas
def tr(A):
	n = len(A)
	s = 0.
	for i in range(n):
		s += A[i][i]
	
	return s

def matmatmul(A, B):
	n = len(A)
	R = [[0.] * n for i in range(n)]

	for i in range(n):
		for j in range(n):
			for k in range(n):
				R[i][j] += A[i][k] * B[k][j]
	
	return R

def matvetmul(A, x):
	n = len(A)
	v = [0.] * n
	for i in range(n):
		for j in range(n):
			v[i] += A[i][j] * x[j]
	
	return v

def scvetmul(scalar, v):
	return [ scalar * v[i] for  i in range(len(v)) ]

def vetvetdiv(w, v):
	return [ w[i] / v[i] for  i in range(len(v)) ]

def vetvetsub(w, v):
	return [ w[i] - v[i] for  i in range(len(v)) ]

def norm(v):
	n = len(v)
	max_v = v[0]
	for i in range(n):
		if (v[i] > max_v):
			max_v = v[i]
	return max_v

def SolveSup(A, b):
	n = len(A)
	x = [0.] * n
	for i in range(n -1, -1, -1):
		s = 0.
		for j in range(i+1, n):
			s+= A[i][j] * x[j]

		if(A[i][i] == 0):
			x[i] = 0
		else:
			x[i] = (b[i] - s) / A[i][i]
	
	return x

def SolveInf(A, b):
	n = len(A)
	x = [0.] * n
	for i in range(n):
		s = 0.
		for j in range(i):
			s+= A[i][j] * x[j]

		if(A[i][i] == 0):
			x[i] = 0
		else:
			x[i] = (b[i] - s) / A[i][i]
	
	return x

def decompLU(A):
	n = len(A)
	L = [ [0.] *n for i in range (n)]
	U = [ [0.] *n for i in range (n)]
	
	for i in range(n):
		L[i][i] = 1.
		for j in range(n):
			if(i <= j):
				s = 0.
				for k in range(i):
					s += L[i][k] * U[k][j]
				U[i][j] = (A[i][j] - s)
			
			else:
				s = 0.
				for k in range(j):
					s += L[i][k] * U[k][j]
				if(U[j][j] == 0):
					L[i][j] = 0
				else:
					L[i][j] = (A[i][j] - s)/U[j][j]
	return L, U

def transpose(A):
	n = len(A)
	At = [ [0.] *n for i in range (n)]
	for i in range(n):
		for j in range(n):
			for k in range(n):
				At[i][j] = A[j][i]
	
	return At

# Questao 1
def PowerMethod(A, epsilon):
	yk = [1. ] * len(A)
	lambdak = [1.]

	while(True):
		zk1 = matvetmul(A, yk) # Passo 1
		tmp = vetvetdiv(zk1, yk) # Passo 2
		
		if (norm(vetvetsub(tmp, lambdak)) / norm(lambdak) < epsilon):
			break

		lambdak = tmp
		alpha = max( [ abs(zk1[i]) for i in range(len(zk1)) ] ) # Passo 3
		yk = scvetmul(1./alpha, zk1)

	return lambdak[0]

A = [[1., -1., 3.], [-1., 1., 3.], [3., -3., 9.]]
I = [[1.,0.,0.], [0.,1.,0.], [0.,0.,1.]]
lambdak = PowerMethod(A, 1e-6)

print("Questao 1:")
print("Maior autovalor:", lambdak)

def InversePowerMethod(A, epsilon):
	n = len(A)
	I = [[0.] * n for i in range(n)]

	for i in range(n):
		I[i][i] = 1.

	AinvT = []

	L, U = decompLU(A)
	for i in range(n):
		y = SolveInf(L, I[i])
		x = SolveSup(U, y)

		AinvT.append(x)

	Ainv = transpose(AinvT)

	yk = [1. ] * len(A)
	lambdak = [1.]

	while(True):
		zk1 = matvetmul(Ainv, yk) # Passo 1
		tmp = vetvetdiv(zk1, yk) # Passo 2
		
		if (norm(vetvetsub(tmp, lambdak)) / norm(lambdak) < epsilon):
			break

		lambdak = tmp
		alpha = max( [ abs(zk1[i]) for i in range(len(zk1)) ] ) # Passo 3
		yk = scvetmul(1./alpha, zk1)

	return lambdak[0]

# Questão 2
def Leverrier(A):
	n = len(A)
	Ak = copy.deepcopy(A)
	P = [0.] *(n+1)
	S = [0.] *(n+1)

	for k in range(1, n+1):
		S[k] = tr(Ak)
		s = 0.
		for i in range(1, k):
			s += P[i] * S[k-i]
			
		P[k] = (S[k] - s)/k
		Ak = matmatmul(Ak, A)
	
	return P[1:n+1]

A = [[1., -3., 3.], [3., -5., 3.], [4., -6., 4.]]
coefs = Leverrier(A)

print('')
print('Questao 2:')
print('coeficientes: ', coefs)
print(f"polinomio caracteristico: -x^3 + {coefs[1]}x + {coefs[2]}")
print('')