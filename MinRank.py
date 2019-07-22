'''
Sage code for creating matrix with polynomial entries

test if rank \leq k
'''
import sys
from sage.all import *

#Test Example
R = PolynomialRing(QQ,3,'x,y,z')
R.inject_variables()
D = {(27,0):2*x+z,(47,0):2,(20,1):y,(8,1):3*z-y}
M = matrix(R,81,2,D,sparse=True)


# Trim down sparse Matrix
def nonzerorows(M):
	I = []
	for i in range(M.nrows()):
		if M.rows()[i] != (0):
			I.append(i)
	return I

def nonzerocols(M):
	Mt = M.transpose()
	J = nonzerorows(Mt)
	return J

def Trim(M):
	I = nonzerorows(M)
	J = nonzerocols(M)
	tM = M[I,:][:,J]
	return tM


# Echelon Form
def PivotM(M):
	NZ = []
	for i in range(M.nrows()):
		for j in range(M.ncols()):
			if M[i][j] in QQ and M[i][j] != 0:
				NZ.append((i,j))
	if len(NZ) == 0:
		k = None
	else:
		k = NZ[0]
	return k
def swap_cols(N,i,j):
	tN = N.transpose()
	tN.swap_rows(i,j)
	tN.transpose()
	return tN

def RowRedOverR(M):
	

def EchelonOverR(M):
	P = PivotM(M)
	if P == None:
		return M
	elif P == (0,0) and M[0,1:] == M[1:,0] == 0:
		N = EchelonOverR(M[1:,1:])
	else:
		
	return N

def Echelon1(M):
	if M[0][0] == 1 and M
# Minors

