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
def Pivots(M,l):
	NZ = []
	for i in range(l,M.nrows()):
		for j in range(l,M.ncols()):
			if M[i][j] in QQ and M[i][j] != 0:
				NZ.append((i,j))
	return NZ


def SwapRow(N,i,j):
	k = N.nrows()
	S = matrix.identity(QQ,k)
	S[i,j] = 1
	S[j,i] = 1
	S[i,i] = 0
	S[j,j] = 0
	N1 = S*N
	return N1

def SwapCol(N,i,j):
	k = N.ncols()
	S = matrix.identity(QQ,k)
	S[i,j] = 1
	S[j,i] = 1
	S[i,i] = 0
	S[j,j] = 0
	N1 = N*S
	return N1

#Multiply row j by a and add it to i
def RowOp(ring,T,i,j,a):
	k = T.nrows()
	S = matrix.identity(ring,k)
	S[i,j] = a
	T1 = S*T
	return T1

#Multiply col j by a and add it to i
def ColOp(ring,T,i,j,a):
	T = T.transpose()
	T1 = RowOp(ring,T,i,j,a)
	T1 = T1.transpose()
	return T1

def MultRow(T,i,a):
	k = T.nrows()
	S = matrix.identity(QQ,k)
	S[i,i] = a
	T1 = S*T
	return T1


def PSmithForm(M,ring):
	P = copy(M)
	l = 0
	p = Pivots(M,l)
	while len(p) != 0:
		i,j = p[0]
		if i != l:
			P = SwapRow(P,i,l)
		if j != l:
			P = SwapCol(P,j,l)
		if P[l,l] != 1:
			P = MultRow(P,l,1/P[l,l])
		for k in range(l+1,M.nrows()):
			P = RowOp(ring,P,k,l,-P[k,l])
		for k in range(l+1,M.ncols()):
			P = ColOp(ring,P,k,l,-P[l,k])
		l += 1
		p = Pivots(P,l)
	r = l
	B = P[r:,r:]
	return r,B	

# Minors for minimum rank

def MinRank(M,ring,s,dimS2AB):
	r,B = PSmithForm(M,ring)
	m = min(M.nrows(),M.ncols())
	if r+1 >= s:
		return False
	else:
		s1 = dimS2AB - s
		alpha = s1+1-r
		K = []
		for k in range(1,alpha):
			K.extend(B.minors(k))
		I = ideal(K).groebner_basis()
		if I == [1]:
			return False
		else:
			return True

		



