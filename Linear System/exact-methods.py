import sys
import copy

def printMat(M):
  n = len(A)
  for i in range(n):
    for j in range(n):
      sys.stdout.write("%.2f\t" %M[i][j])
    print("")

def printExtMat(M, b):
  n = len(M)
  for i in range(n):
    for j in range(n):
      sys.stdout.write("%.2f\t" %M[i][j])
    sys.stdout.write("| %.2f\t" %b[i])
    print("")

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

def MulMatVec(A, x):
  n = len(A)
  v = [0.] * n
  for i in range(n):
    for j in range(n):
      v[i] += A[i][j] * x[j]
  
  return v

def mulMat(A,B):
  n = len(A)
  C = [ [0.] *n for i in range (n)]
  for i in range(n):
    for j in range(n):
      for k in range(n):
        C[i][j] += A[i][k] * B[k][j]
  
  return C 

def transpose(A):
  n = len(A)
  At = [ [0.] *n for i in range (n)]
  for i in range(n):
    for j in range(n):
      for k in range(n):
        At[i][j] = A[j][i]
  
  return At

def UpdateVectorB(A, b):
  n = len(b)
  for k in range(n):
    for i in range(k+1, n):
      pivot = A[i][k]
      b[i] = b[i] - pivot * b[k]

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

def GaussElim(A, b):
  # k = 0,..., n-2
  n = len(A)
  
  for k in range(0, n-1):
    # i = k+1,..., n-1
    for i in range(k+1, n):
      p = -A[k][k] / A[i][k]
      b[i] = p*b[i] + b[k] # Matriz estendida
      
      # j = k,..., n-1
      for j in range(k, n):
        A[i][j] = p*A[i][j] + A[k][j]

  return A, b

def GaussElimPartialPivot(A, b):
  n = len(A)
  
  for k in range(0, n-1): # k = 0,..., n-2
    # trocar as linhas de acordo com o modulo
    # Maior modulo sera o elemento da diagonal
    max_v = abs(A[k][k])
    max_idx = k

    for i in range(k+1, n):
      if(abs(A[i][k]) > max_v):
        max_v = abs(A[i][k])
        max_idx = i

    # troca a linha k pela linha max_idx
    if(max_idx != k):
      for j in range(n):
        tmp = A[k][j]
        A[k][j] = A[max_idx][j]
        A[max_idx][j] = tmp
      
      tmp = b[k]
      b[k] = b[max_idx]
      b[max_idx] = tmp


    for i in range(k+1, n): # i = k+1,..., n-1
      p = A[i][k] / A[k][k]
      b[i] = b[i] - p*b[k] # Matriz estendida ###########
      
      for j in range(k, n): # j = k,..., n-1
        A[i][j] = A[i][j] - p*A[k][j]

  return A, b

def GaussElimPartialPivotCompact(A, b):
  n = len(A)
  
  for k in range(0, n-1): # k = 0,..., n-2
    # trocar as linhas de acordo com o modulo
    # Maior modulo sera o elemento da diagonal
    max_v = abs(A[k][k])
    max_idx = k

    for i in range(k+1, n):
      if(abs(A[i][k]) > max_v):
        max_v = abs(A[i][k])
        max_idx = i

    # troca a linha k pela linha max_idx
    if(max_idx != k):
      for j in range(n):
        tmp = A[k][j]
        A[k][j] = A[max_idx][j]
        A[max_idx][j] = tmp
      
      tmp = b[k]
      b[k] = b[max_idx]
      b[max_idx] = tmp


    for i in range(k+1, n): # i = k+1,..., n-1
      p = A[i][k] / A[k][k]
      
      for j in range(k, n): # j = k,..., n-1
        A[i][j] = A[i][j] - p*A[k][j]
      A[i][k] = p

  return A, b

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

def norm(v):
  n = len(v)
  max_v = v[0]
  for i in range(n):
    if (v[i] > max_v):
      max_v = v[i]
  return max_v

def SolutionRefine(A, b, err):
  Acopy = copy.deepcopy(A)
  bcopy = copy.deepcopy(b)
  GaussElimPartialPivotCompact(A, b)
  UpdateVectorB(A, b)
  x = SolveSup(A, b)
  r = SubVec(bcopy, MulMatVec(Acopy, x))
  print(f"Residuo: {r}")

  while(norm(r) > err):
    UpdateVectorB(A, r)
    y = SolveSup(A, r)
    x = AddVec(x, y)
    r = SubVec(bcopy, MulMatVec(Acopy, x))
    print(f"Residuo: {r}")
  
  return x

"""
# Exemplo 1 Refinamento
Acopy = [[16., 5.], [3., 2.5]]
A = [[16., 5.], [3., 2.5]]
b = [21.,5.5]

# Exemplo 2 Refinamento
Acopy = [[1.,3.,4.], [3.,2.,1.], [2.,4.,3.]]
bcopy = [-5.,8.,4.]

A = [[1.,3.,4.], [3.,2.,1.], [2.,4.,3.]]
b = [-5.,8.,4.]

print (SolutionRefine(A, b, 10e-2))
"""

"""
# Matriz Inversa
n = 3
A = [[1., -1., 3.], [-2., -2., -1.], [4., 3., 8.]]
b = [[1.,0,0], [0,1.,0], [0,0,1.]]

AinvT = []

L, U = decompLU(A)
for i in range(n):
  y = SolveInf(L, b[i])
  x = SolveSup(U, y)

  AinvT.append(x)

Ainv = transpose(AinvT)
printMat (mulMat(A, Ainv))
"""