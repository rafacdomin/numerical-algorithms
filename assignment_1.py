import sys
import copy

# Funcoes usadas:

def printExtMat(M, b):
  n = len(M)
  for i in range(n):
    for j in range(n):
      sys.stdout.write("%.2f\t" %M[i][j])
    sys.stdout.write("| %.2f\t" %b[i])
    print("")

def norm(v):
  n = len(v)
  max_v = v[0]
  for i in range(n):
    if (v[i] > max_v):
      max_v = v[i]
  return max_v

def SubVec(v1, v2):
  n = len(v1)
  r = []
  
  for i in range(n):
    r.append(v1[i] - v2[i])
  
  return r 

def AddVec(v1, v2):
  n = len(v1)
  r = []
  
  for i in range(n):
    r.append(v1[i] + v2[i])
  
  return r

def transpose(A):
  n = len(A)
  At = [ [0.] *n for i in range (n)]
  for i in range(n):
    for j in range(n):
      for k in range(n):
        At[i][j] = A[j][i]
  
  return At

def MulMatVec(A, x):
  n = len(A)
  v = [0.] * n
  for i in range(n):
    for j in range(n):
      v[i] += A[i][j] * x[j]
  
  return v

def SolveSup(A, b):
  n = len(A)
  x = [0.] * n
  for i in range(n -1, -1, -1):
    s = 0.
    for j in range(i+1, n):
      s+= A[i][j] * x[j]

    x[i] = (b[i] - s) / A[i][i]
  
  return x

def SolveInf(A, b):
  n = len(A)
  x = [0.] * n
  for i in range(n):
    s = 0.
    for j in range(i):
      s+= A[i][j] * x[j]

    x[i] = (b[i] - s) / A[i][i]
  
  return x

def GaussElimPartialPivot(A, b):
  n = len(A)
  
  for k in range(0, n-1):
    max_v = abs(A[k][k])
    max_idx = k

    for i in range(k+1, n):
      if(abs(A[i][k]) > max_v):
        max_v = abs(A[i][k])
        max_idx = i

    if(max_idx != k):
      for j in range(n):
        tmp = A[k][j]
        A[k][j] = A[max_idx][j]
        A[max_idx][j] = tmp
      
      tmp = b[k]
      b[k] = b[max_idx]
      b[max_idx] = tmp


    for i in range(k+1, n):
      p = A[i][k] / A[k][k]
      b[i] = b[i] - p*b[k]
      
      for j in range(k, n):
        A[i][j] = A[i][j] - p*A[k][j]

  return A, b

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
        L[i][j] = (A[i][j] - s)/U[j][j]
  return L, U



# Questao 1:
def SolutionRefinePivot(A, b, err):
  Acopy = copy.deepcopy(A)
  bcopy = copy.deepcopy(b)

  GaussElimPartialPivot(A, b)
  x = SolveSup(A, b)
  r = SubVec(bcopy, MulMatVec(Acopy, x))
  print(f"Residuo: {r}, norm(r): {norm(r)}")

  while(norm(r) > err):
    A_ = copy.deepcopy(Acopy)
    GaussElimPartialPivot(A_, r)
    y = SolveSup(A_, r)
    x = AddVec(x, y)
    r = SubVec(bcopy, MulMatVec(Acopy, x))
    print(f"Residuo: {r}, norm:{norm(r)}")
  
  return x

print("Questao 1:")
print("\nSistema Linear:")
A = [[1.,4.,1.], [3.,1.,-1.], [-5., 13., -22.]]
b = [7., 3., 48.]
printExtMat(A, b)

print("\nResiduos refinamento com pivoteamento parcial:")
x = SolutionRefinePivot(A, b, 10e-16)

print(f"\nSolucao Refinada: {x}")

# Questao 2
print("\n\nQuestao 2:")
print("\nSistema Linear:")
A = [[1.,4.,1.], [3.,1.,-1.], [-5., 13., -22.]]
b = [7., 3., 48.]
printExtMat(A, b)

## A) Numero de condicionamento
print("\nLetra A:")

def InvMat(A):
  n = len(A)
  AinvT = []
  L, U = decompLU(A)
  I = [[1.,0,0], [0,1.,0], [0,0,1.]]

  for i in range(n):
    y = SolveInf(L, I[i])
    x = SolveSup(U, y)

    AinvT.append(x)

  Ainv = transpose(AinvT)
  return Ainv

def normMatInfinite(M):
  n = len(M)
  
  max_v = 0
  for i in range(n):
    s = 0.
    for j in range(n):
      s += abs(M[i][j])
    
    if(max_v < s):
      max_v = s
  
  return max_v

A = [[1.,4.,1.], [3.,1.,-1.], [-5., 13., -22.]]
b = [7., 3., 48.]

Ainv = InvMat(A)

normA = normMatInfinite(A)
normAinv = normMatInfinite(Ainv)

condA = normA * normAinv
print(f"Numero de Condicionamento: {condA}")

## B) Solucao
print("\nLetra B:")
A = [[1.,4.,1.], [3.,1.,-1.], [-5., 13., -22.]]
b = [7., 3., 48.]

A_, b_ = GaussElimPartialPivot(A, b)
print("Eliminacao Gaussiana com Pivoteamento Parcial:")
printExtMat(A_, b_)

x = SolveSup(A_, b_)
print(f"\nSolucao: {x}")