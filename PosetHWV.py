import numpy as np
import sympy as sp
import queue,itertools,math
import scipy.sparse
from scipy.sparse import csr_matrix,vstack,hstack
import BorderApolarity
import PartialSmithForm
from collections import deque

import sys
from sage.all import *

#import sms


#L = LatticeElts(V,LowOps)
#E = L.keys()
#C = CartanMatrix(['A',n-1])
#Cinv = np.linalg.inv(C)
#f = lambda x,y: all(i >= 0 for i in map(int,round(Cinv.dot(np.subtract(y,x)))))
#P = Poset((E,f)) 

# Poset above

def LinExt(L,P):
	LE = list(P.linear_extension())
	LE.reverse()
	N = [L[s].shape[1] for s in Lin]
	return LE,N

def dimker(LowOps,n):
	if len(LowOps) == 0:
		return None
	else:
		M = []
		for k in LowOps:
			M += BorderApolarity.SparseMatLowOp(k,n)
		M1 = vstack(M)
		M2 = sp.Matrix(M1.todense())
		null = M2.nullspace()
		return len(null)

def subsets(s):
	l = []
	for i in range(1,len(s)+1): 
		l += list(map(tuple,itertools.combinations(s,i)))
	return l

 

#dimker = {Lowering ops intersections as tuple:dimker}
# K = [k_1, \cdots, k_l], change k_i so check conditions for k_i
#Lin is linear extension, N is maximal elts, P is poset and DIMKER is dimensions of all kernels
'''
def conditions(K,LE,N,P,DIMKER):
	Test = True
	while Test != False:
		if all(N[i] >= K[i] for i in range(len(N))):
			for j in range(len(K)):
				if K[j] != 0:
					parents = P.upper_covers(LE[j])
					if len(parents) != 0:
						pi = [LE.index(p) for p in parents]
						lo = [np.subtract(p,LE[j]).tolist().index(2) for p in parents]
						lo.sort()
						for s in subsets(lo):
							sum = DIMKER[s]
							for i in s:
								sum += K[pi[lo.index(i)]]
							if K[j] <= sum:
								Test = Test
							else:
								Test = False
		else:
			Test = False
	return Test
'''

def cond(K,i,LE,P,DIMKER):
	Test = True
	parents = P.upper_covers(LE[i])
	if len(parents) != 0:
		pi = [LE.index(p) for p in parents]
		lo = [np.subtract(p,LE[i]).tolist().index(2) for p in parents]
		lo.sort()
		for s in subsets(lo):
			sum = DIMKER[s]
			for j in s:
				sum += K[pi[lo.index(j)]]
			if K[i] <= sum:
				Test = Test
			else:
				Test = False
	return Test
def SUM(q):
	S= 0
	for k in range(len(q)):
		S += q[k]
	return S

def dfs(p,LE,NN,P,DIMKER):
	Hwvs = []
	q = [1]
	q.extend([0]*(len(LE)-1))
	Stack = deque([q])
	while len(Stack) != 0:
		s = Stack.pop()
		if SUM(s) == p:
			Hwvs.append(s)
		else:
			i = [index for index,item in enumerate(s) if item != 0][-1]
			parents = P.upper_covers(LE[i])
			par1 = []
			for q in parents:
				par1.extend(P.lower_covers(q))
			par1 = list(set(par1))
			cci = [LE.index(q) for q in par1 if LE.index(q) >= i]
			child = P.lower_covers(LE[i])
			ci = [LE.index(c) for c in child]
			ci.extend(cci)
			if len(ci) != 0:
				for k in ci:
					r = s[:]
					r[k] += 1
					if (r[k] <= NN[k]) and (cond(r,k,LE,P,DIMKER) == True):
						Stack.append(r)
	return Hwvs

#auxiliary max fcn
def Max(b):
	m = 0
	for i in range(len(b)):
		if b[i][1]-b[i][0] > m:
			m = b[i][1]-b[i][0]
		else:
			m = m
	return m
#R = PolynomialRing(QQ, (n-k)*k,x)

def GrassCharts(k,n,pivots,ring,q,aa,bb):
	assert len(pivots) == k and sorted(pivots) == list(pivots)
	M = matrix(ring,k,n)
	taken = []
	for (row,pivot) in enumerate(pivots):
		M[row,pivots] = 1
		taken.extend((i,pivot) for i in range(k))
		taken.extend((row,j) for j in range(pivot))
	indet_indices = [(i,j) for i in range(k) for j in range(n) if not (i,j) in taken]
	for (idx,y) in zip(indet_indices,ring.gens()[q*bb:(q+1)*bb]):
		M[idx] = y
	return M


def GrassCharts1(k,n,ring,q,aa,bb):
	A = []
	for pivots in itertools.combinations(range(n),k):
		A.extend(GrassCharts(k,n,pivots,ring,q,aa,bb))
	return A 


def wvs(T,N,L,LE):
	if len(T) == 0:
		A,a,b = [None],[],[]
	else:
		a = []
		b = []
		for i in range(len(T)):
			if T[i] < N[i] and T[i] > 0:
				a.append(i)
				b.append((T[i],N[i]))
		A = []
		for j in range(len(T)):
			if T[j] != 0:
				A.append(L[LE[j]])
			else:
				A.append(None)
	return A,a,b

#use if len(a) == 0
def wvs0(A,dimS2AB,r,n):
	B = hstack([a for a in A if a != None])
	B = S2AB(B,n**2-1)
	rk = sms.rank(B)
	return rk

#use if len(a) != 0
def wvs1(A,dimS2AB,r,n,ring):
	r1 = dimS2AB - r
	B = S2AB1(A,n**2-1,ring)
	t = PartialSmithForm.MinRank(B,ring,s,dimS2AB)
	return t
#Convert from csr_matrix to sage sparse matrix
def Convert(A):
	I,J,K = scipy.sparse.find(A)
	dat = {(i,j):v for i,j,v in zip(I,J,K)}
	A1 = matrix(A.shape[0],A.shape[1],dat,sparse=True)
	A1 = A1.sparse_matrix()
	return A1

def Convert1(A):
	for i in range(len(A)):
		if A[i] != None:
			A[i] = Convert(A[i])
	return A
#Multiply A csr_matrix with B a sage Matrix
def Multiply(A,B):
	I,J,K = scipy.sparse.find(A)
	dat = {(i,j):v for i,j,v in zip(I,J,K)}
	A1 = matrix(A.shape[0],A.shape[1],dat,sparse=True)
	B1 = A1*B
	B1 = B1.sparse_matrix()
	return B1

def shstack(A,ring):
	B = [a for a in A if a != None]
	B1 = [b.dict() for b in B]
	D = {}
	m = 0
	for i in range(len(B1)):
		for (j,k) in B1[i]:
			D[(j,k+m)] = B1[i][(j,k)]
		m += B[i].ncols()
	rows = B[0].nrows()
	A1 = matrix(ring,rows,m,D)
	return A1

def DictS2AB(m):
	Dict = {}
	for k in range(m**2):
		q = int(np.floor(k/m))
		r = k % m
		if q == 0:
			Dict[k] = range(k*m,(k+1)*m)
		else:
			a = [Dict[m*i+r][q] for i in range(q)]
			s = np.array([(m-i)*m for i in range(q)]).sum()
			b = range(s+r*(m-q),s+(r+1)*(m-q))
			Dict[k] = a + b
	return Dict

'''
Given E corresponding hyperplane annihilating it, return A* times E
m = dim sln
E is csr_matrix
'''
def S2AB(E,m):
	p,q = E.shape
	I,J,D = scipy.sparse.find(E)
	I1 = []
	J1 = []
	D1 = []
	for i in list(I):
		for k in range(m):
			I1.extend(Dict[i])
	for j in list(J):
		for k in range(m):
			J1.extend([m*j+k])
	for d in list(D):
		D1.extend([d]*m)
	M = csr_matrix((D1,(I1,J1)),shape=(p*m,q*m))
	return M


def S2AB1(E,m,ring):
	p = E.nrows()
	q = E.ncols()
	D = E.dict()
	D1 = {}
	for k in range(m):
		for (i,j) in D:
			D1[(m*i+k,m*j+k)] = D[(i,j)]
	M = matrix(ring,p*m,q*m,D1,sparse=True)
	return M

#A is csr_matrix, removes zero rows and zero columns
def trim(A):
	B = A[A.getnnz(1)>0][:,A.getnnz(0)>0]
	return B

#A is matrix over a polynomial ring
def minRK(A):
	m = min(A.nrows(),A.ncols())
	k = 1
	while k <= m:
		I = ideal(A.minors(k)).groebner_basis()
		if I != [1]:
			break
		k += 1
	return k-1



def minRK1(A,r):
	I = ideal(A.minors(r+1)).groebner_basis()
	if I != [1]:
		t = False
	else:
		t = True
	return t




