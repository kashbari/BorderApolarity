import numpy as np
import math
from scipy.sparse import csr_matrix,identity

#X is m^2 by m^2-r, with so write in basis for A*B*C* (m^3), so want m^3 by m*(m^2-r)
def csr2matsage(X):
	I,J,K = scipy.sparse.find(X)
	dat = {(i,j):v for i,j,v in zip(I,J,K)}
	X1 = matrix(X.shape[0],X.shape[1],dat,sparse=True)
	return X1

def Sep1(D,k,m):
	for d in D:
		k0 = int(math.floor(d[0]/k))
		k1 = d[0]%k
		if k0 != 0
			D[(k*m*k0+k1,d[1])] = D[d]
			del D[d]
		return D

def Shft1(D,k,m):


#Want X in matrix format
def E110(X,m):
	if type(X) == csr_matrix:
		X = csr2matsage(X)
	else:
		X = X
	X1 = kron(identity(m),X)
	return X1


def E101(X,m):



def E011(X,m):
