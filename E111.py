import numpy as np
import scipy.sparse
import math
from scipy.sparse import csr_matrix,identity,kron
import scipy

import sys
from sage.all import *

#X is m^2 by m^2-r, with so write in basis for A*B*C* (m^3), so want m^3 by m*(m^2-r)
def csr2matsage(X):
	I,J,K = scipy.sparse.find(X)
	dat = {(i,j):v for i,j,v in zip(I,J,K)}
	X1 = matrix(QQ,X.shape[0],X.shape[1],dat,sparse=True)
	return X1

def Sep1(D,k,l):
	E = {}
	for d in D:
		q0 = int(math.floor(d[0]/k))
		q1 = d[0]%k
		E[(q0*l+q1,d[1])] = D[d]
	return E

def Shft1(D,k,m,q):
	E = {}
	for d in D:
		for t in range(m):
			E[(d[0]+k*t,d[1]+q*t)] = D[d]
	return E

#Want X in matrix format
def E110(X,m):
	if type(X) == csr_matrix:
		X = csr2matsage(X)
	q = X.ncols()
	D = X.dict()
	D = Sep1(D,m**2,m**2)
	D = Shft1(D,m**2,m,q)
	X1 = matrix(m**3,m*X.ncols(),D,sparse=True)
	#X1 = kron(identity(m),X)
	return X1


def E101(X,m):
	if type(X) == csr_matrix:
		X = csr2matsage(X)
	q = X.ncols()
	D = X.dict()
	D = Sep1(D,m,m**2)
	D = Shft1(D,m,m,q)
	X1 = matrix(m**3,m*X.ncols(),D,sparse=True)
	return X1

def E011(X,m):
	if type(X) == csr_matrix:
		X = csr2matsage(X)
	q = X.ncols()
	D = X.dict()
	D = Sep1(D,1,m)
        D = Shft1(D,1,m,q)
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
	zcol = X1.ncols()+Y1.ncols()
	for d in DY:
		DX[(d[0],d[1]+ycol)] = DY[d]
	for d in DZ:
		DX[(d[0],d[1]+zcol)] = DZ[d]
	cols = X1.ncols()+Y1.ncols()+Z1.ncols()
	rows = X1.nrows()
	M = matrix(rows,cols,DX,sparse=True)
	return M


	
	

