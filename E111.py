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
		print(d)
		if k0 != 0:
			D[(k*m*k0+k1,d[1])] = D[d]
			print((k*m*k0+k1,d[1]))
			D.pop(d)
			print(len(D))
	return D

def Shft1(D,k,m):
	return D

#Want X in matrix format
def E110(X,m):
	if type(X) == csr_matrix:
		X = csr2matsage(X)
	X1 = kron(identity(m),X)
	return X1


def E101(X,m):
	if type(X) == csr_matrix:
		X = csr2matsage(X)
	D = X.dict()
	D = Sep1(D,m,m)
	D = Shft1(D,m,m)
	X1 = matrix(m**3,m*X.ncols(),D,sparse=True)
	return X1

def E011(X,m):
	if type(X) == csr_matrix:
		X = csr2matsage(X)
	D = Sep1(D,1,m)
        D = Shft1(D,1,m)
        X1 = matrix(m**3,m*X.ncols(),D,sparse=True)
	return X1

def E111(X,Y,Z,m):
	X1 = E110(X,m)
	Y1 = E101(Y,m)
	Z1 = E011(Z,m)
	DX = X1.dict()
	DY = Y1.dict()
	DZ = Z1.dict()
	ycol = X1.ncols()
	zcol = X1.ncols()+Y1.ycols()
	for d in DY:
		DY[(d[0],d[1]+ycol)] = DY[d]
		del DY[d]
	for d in DZ:
		DZ[(d[0],d[1]+zcol)] = DZ[d]
		del DZ[d] 
	D = DX|DY|DZ
	cols = X1.ncols()+Y1.ncols()+Z1.ncols()
	rows = X1.nrows()+Y1.nrows()+Z1.nrows()
	M = matrix(rows,cols,D,sparse=True)
	return M


	
	

