import scipy
from scipy.sparse import dok_matrix
import numpy as np
from numpy.linalg import lstsq
import tensorly as tl
from tensorly.decomposition import parafac

def e(i,j,n):
    M = dok_matrix((n,n))
    M[i,j] = 1
    return M 

#return matrix(QQ,n,n,{(i,j):1})

def basis(i,n):
    assert i < n**2-1
    i,j = i//n, i%n
    if i==j:
        return e(i,i,n) - e(i+1,i+1,n)
    else:
        return e(i,j,n)
'''
# e_ij
def Tsln(n):
    T = [{} for i in range(n**2-1)]
    B = matrix(QQ,[ list(basis(i,n)) for i in range(n**2-1)])
    for i in range(n**2-1):
        a = basis(i,n)
        for j in range(n**2-1):
            b = basis(j,n)
            c = a*b - b*a
            c = np.linalg.lstsq(B,list(c))
            #c = B.solve_left(vector(QQ,c.list()))
            for k in c.nonzero_positions():
                T[i][(j,k)] = c[k]
    T = [matrix(QQ,n**2-1,n**2-1,m) for m in T]
    return T
'''

D0 = {(1,1):2, (2,2):1, (3,3):-2, (5,5):-1, (6,6):-1, (7,7):1}
D1 = {(0,1):-2, (3,0):1, (4,1):1, (5,2):1, (6,7):-1}
D2 = {(0,2):-1, (3,5):-1, (4,2):-1, (6,0):1, (6,4):1, (7,1):1}
D3 = {(0,3):2, (1,0):-1, (2,5):1, (4,3):-1, (7,6):-1}
D4 = {(1,1):-1, (2,2):1, (3,3):1, (5,5):2, (6,6):-1, (7,7):-2}
D5 = {(0,5):1, (1,2):-1, (4,5):-2, (6,3):1, (7,4):1}
D6 = {(0,6):1, (1,7):1, (2,0):-1, (2,4):-1, (4,6):1, (5,3):-1}
D7 = {(0,7):-1, (2,1):-1, (3,6):1, (4,7):2, (5,4):-1}
Dict = [D0,D1,D2,D3,D4,D5,D6,D7]
