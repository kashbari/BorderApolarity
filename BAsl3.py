import BorderApolarity
import PosetHWV

import numpy as np
import sympy as sp
import queue,itertools,math
from scipy.sparse import csr_matrix, vstack, hstack
from collections import deque

import sys
from sage.all import *

n = 3
r  = 15
p = r - (n**2-1)
dimS2AB = int( (n**2)*(n**2-1)*(n**2-1)/2)


row = [20]
col = [0]*len(row)
data = [1]
V22 = csr_matrix((data,(row,col)),shape=(n**4,1))


row = [11,19]
col = [0]*len(row)
data = [-1,1]
V30 = csr_matrix((data,(row,col)),shape=(n**4,1))


row = [23,47]
col = [0]*len(row)
data = [-1,1]
V03 = csr_matrix((data,(row,col)),shape=(n**4,1))


row = [0,4,8,12,24,28,36,40,44,52,56,68,72,76,80]
col = [0]*len(row)
data = [2,-1,-1,3,3,3,-1,2,-1,3,3,3,-1,-1,2]
V00 = csr_matrix((data,(row,col)),shape=(n**4,1))

row = [2,14,18,22,26,38,46,74]
col = [0]*len(row)
data = [3,9,3,-6,3,-6,9,3]
V11 = csr_matrix((data,(row,col)),shape=(n**4,1))




LowOps = [BorderApolarity.SparseMatLowOp(i,n) for i in range(n-1)]
C = CartanMatrix(['A',n-1])
Cinv = np.linalg.inv(C)
f = lambda x,y: all(i >= 0 for i in map(int,round(Cinv.dot(np.subtract(y,x)))))


def WeightPoset(V,LowOps,f):
	L = BorderApolarity.LatticeElts(V,LowOps)
	E = L.keys()
	P = Poset((E,f))
	LE = list(P.linear_extension())
	LE.reverse()
	N = [L[s].shape[1] for s in LE]
	return L,LE,P,N

Dict22,LE22,P22,N22 = WeightPoset(V22,LowOps,f)

Dict30,LE30,P30,N30 = WeightPoset(V30,LowOps,f)

Dict03,LE03,P03,N03 = WeightPoset(V03,LowOps,f)

Dict00,LE00,P00,N00 = WeightPoset(V00,LowOps,f)

Dict11,LE11,P11,N11 = WeightPoset(V11,LowOps,f)

#precomputed
DIMKER = {}
DIMKER[(0,)] = 35
DIMKER[(1,)] = 35
DIMKER[(0,1)] = 9

C1 = [len(LE22),len(LE30),len(LE03),len(LE00),len(LE11)]
C = [c for c in BorderApolarity.findComp(p,5) if all( i >= 0 for i in np.subtract(C1,c))]


def HWVGrassmannian(c):
	H22=PosetHWV.dfs(c[0],LE22,N22,P22,DIMKER)
	H30=PosetHWV.dfs(c[1],LE30,N30,P30,DIMKER)
	H03=PosetHWV.dfs(c[2],LE03,N03,P03,DIMKER)
	H00=PosetHWV.dfs(c[3],LE00,N00,P00,DIMKER)
	H11=PosetHWV.dfs(c[4],LE11,N11,P11,DIMKER)
	return H22,H30,H03,H00,H11


def AnnPlane(c,dimS2AB,r):
	H22,H30,H03,H00,H11 = HWVGrassmannian(c)
	if len(H22) == 0:
		H22.append(None)
	if len(H30) == 0:
		H30.append(None)
	if len(H03) == 0:
		H03.append(None)
	if len(H00) == 0:
		H00.append(None)
	if len(H11) == 0:
		H11.append(None)
	for z in itertools.product(H22,H30,H03,H00,H11):
		K = []
		if z[0] != None:
			K.append(np.subtract(N22,z[0]))
		else:
			K.append(None)
		if z[1] != None:
			K.append(np.subtract(N30,z[1]))
		else:
			K.append(None)
		if z[2] != None:
			K.append(np.subtract(N03,z[2]))
		else:
			K.append(None)
		if z[3] != None:
			K.append(np.subtract(N00,z[3]))
		else:
			K.append(None)
		if z[4] != None:
			K.append(np.subtract(N11,z[4]))
		else:
			K.append(None)
		A22,a22,b22 = wvs(K[0],N22,Dict22,LE22)
		A30,a30,b30 = wvs(K[1],N30,Dict30,LE30)
		A03,a03,b03 = wvs(K[2],N03,Dict03,LE03)
		A00,a00,b00 = wvs(K[3],N00,Dict00,LE00)
		A11,a11,b11 = wvs(K[4],N11,Dict11,LE11)
		if a22 == a30 == a03 == a00 == a11 == []:
			A = A22 + A30 + A03 + A00 + A11
			t,rk = wvs0(A,dimS2AB,r)
			if t == True:
				print(c,K,rk)
		else:
			A22 = Convert1(A22)
			A30 = Convert1(A30)
			A03 = Convert1(A03)
			A00 = Convert1(A00)
			A11 = Convert1(A11)
			aa = len(a22) + len(a30) + len(a03) + len(a00) + len(a11)
			a = a22+ a30+a03+a00+a11
			bb = Max(b22+ b30+ b03+ b00 +b11)
			R = PolynomialRing(QQ,['x_%d%d' %(i,j) for i in range(aa) for j in range(bb)])
	return 

		

 
