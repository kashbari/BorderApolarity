import numpy as np
import sympy as sp
import queue,itertools,math
import scipy.sparse
from scipy.sparse import csr_matrix,vstack,hstack,kron
from scipy.sparse.linalg import inv
import BorderApolarity
import PartialSmithForm
import E111
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

def dimker(w,Dict,LE,LowOps,n):
	if len(LowOps) == 0:
		return None
	else:
		M = []
		for k in LowOps:
			M += csr_matrix.transpose(BorderApolarity.SparseMatLowOp(k,n))
		M1 = vstack(M)
		M3 = Dict[w]
		M2 = M1*M3
		I,J,K = scipy.sparse.find(M2)
		J = list(dict.fromkeys(J))
		ker = int(Dict[w].shape[1]- len(J))
		return ker 

def DimKer(P,Dict,LE,n):
	DIMKER = {}
	for p in P:
		if P.upper_covers(p) != []:
			pi = [LE.index(q) for q in P.upper_covers(p)]
			lo = [np.subtract(LE[q],p).tolist().index(2) for q in pi]
			S = subsets(lo)
			for s in S:
				DIMKER[(p,s)] = dimker(p,Dict,LE,s,n)
	return DIMKER

def subsets(s):
	l = []
	for i in range(1,len(s)+1): 
		l += list(map(tuple,itertools.combinations(s,i)))
	l1 = []
	for k in range(len(l)):
		l1.append(tuple(sorted(l[k])))
	return l1
 

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
			sum = DIMKER[(LE[i],s)]
			for j in s:
				sum += K[pi[lo.index(j)]]
			if K[i] <= sum:
				Test = Test
			else:
				Test = False
	return Test

def firstnonzero(r):
	I = [index for index,item in enumerate(r) if item != 0][-1]
	if I == 1 :
		return True
	else:
		return False

def SUM(q):
	S= 0
	for k in range(len(q)):
		S += q[k]
	return S

def startdfs(S,q):
	S1 = []
	for s in S:
		q1 = [0]*q
		q1[s] = 1
		S1.extend([q1])
	return S1


def dfs(p,LE,NN,P,DIMKER):
	Hwvs = []
	q = [0]*(len(LE))
	Stack = deque([q])
	while len(Stack) != 0:
		s = Stack.pop()
		if SUM(s) == p:
			Hwvs.append(s)
		else:
			I = [index for index,item in enumerate(s) if item != 0]
			if len(I) == 0:
				for k in range(len(s)):
					r = s[:]
					r[k] += 1
					if (r[k] <= NN[k]) and (cond(r,k,LE,P,DIMKER) == True):
						Stack.append(r)
			else:
				i = I[-1]	
				#parents = P.upper_covers(LE[i])
				#par1 = []
				#for q in parents:
				#	par1.extend(P.lower_covers(q))
				#par1 = list(set(par1))
				#cci = [LE.index(q) for q in par1 if LE.index(q) >= i]
				#child = P.lower_covers(LE[i])
				#ci = [LE.index(c) for c in child]
				#ci.extend(cci)
				ci = range(i,len(LE))
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
		if b[i][0]*(b[i][1]-b[i][0]) > m:
			m = b[i][0]*(b[i][1]-b[i][0])
		else:
			m = m
	return m
#R = PolynomialRing(QQ, (n-k)*k,x)

def GrassCharts(k,n,pivots,ring,q,bb):
	assert len(pivots) == k and sorted(pivots) == list(pivots)
	M = Matrix(ring,k,n)
	taken = []
	for (row,pivot) in enumerate(pivots):
		M[row,pivot] = 1
		taken.extend((i,pivot) for i in range(k))
		taken.extend((row,j) for j in range(pivot))
	indet_indices = [(i,j) for i in range(k) for j in range(n) if not (i,j) in taken]
	for (idx,y) in zip(indet_indices,ring.gens()[q*bb:(q+1)*bb]):
		M[idx] = y
	return M


def GrassCharts1(k,n,ring,q,bb):
	A = []
	for pivots in itertools.combinations(range(n),k):
		A.append(matrix(GrassCharts(k,n,pivots,ring,q,bb)))
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
'''
#use if len(a) == 0
def wvs0(A,dimS2AB,r,n,DS2AB):
	B = hstack([a for a in A if a != None])
	S = COB(n)
	B = S.dot(B)
	print(B.shape)
	B = S2AB(B,n**2-1,DS2AB)
	rk = sms.rank(B)
	return rk

def csr2matsage(X):
	I,J,K = scipy.sparse.find(X)
	dat = {(i,j):v for i,j,v in zip(I,J,K)}
	X1 = matrix(X.shape[0],X.shape[1],dat,sparse=True)
	return X1

'''
def wvs01(A,dimS2AB,r,n,DS2AB):
	B = hstack([a for a in A if a != None])
	S = COB(n)
	B = S.dot(B)
	B = S2AB(B,n**2-1,DS2AB)
	B1 = E111.csr2matsage(B)
	rk = B1.rank()

def wvs01M(A,dimS2AB,r,n,DS2AB):
	B = hstack([a for a in A if a != None])
	S = COB(n)
	B = S.dot(B)
	B = S2AB(B,n**2-1,DS2AB)
	B1 = E111.csr2matsage(B)
	return B1

#use if len(a) != 0
def wvs1(A,dimS2AB,r,n,ring,DS2AB):
	r1 = dimS2AB - r
	B = S2AB1(A,n**2-1,DS2AB,ring)
	print(B.nrows(),B.ncols())
	#r,B1 = PartialSmithForm.PSmithForm(B,ring)
	t = PartialSmithForm.MinRank(B,ring,r,dimS2AB)
	return t

def wvs1M(A,dimS2AB,r,n,ring,DS2AB):
	B = S2AB1(A,n**2-1,DS2AB,ring)
	return B
#COMBINE DICTIONARIES FOR POSET; D is list of Dictionaries
def Combine(Dict1,LE1,N1,ListofDic):
	for D in ListofDic:
		assert type(D) == dict
		for d in D:
			if d in Dict1:
				Dict1[d] = hstack([Dict1[d],D[d]])
			else:
				Dict1[d] = D[d]
	for d in LE1:
		N1[LE1.index(d)] = Dict1[d].shape[1]
	return Dict1,N1 
			
				


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

#M is csr_matrix, 
def COB(n):
	row = range(n**2-1)
	col = range(n**2-1)
	data = [1]*len(row)
	for k in range(0,n-1):
		col.append(k*(n+1))
		row.append((k+1)*(n+1))
		data.append(-1)
	B = csr_matrix((data,(row,col)),shape=(n**2,n**2-1))
	A = kron(B,B)
	tA = A.transpose()
	C = inv(tA.dot(A))
	S = C.dot(tA)
	return S

def COB1(B,n):
	S = COB(n)
	B = Multiply(S,B).sparse_matrix()
	return B



'''	

Given E corresponding hyperplane annihilating it, return A* times E
m = dim sln
E is csr_matrix
'''
def S2AB(E,m,DS2AB):
	p,q = E.shape
	I,J,D = scipy.sparse.find(E)
	I1 = []
	J1 = []
	D1 = []
	for i in list(I):
		I1.extend(DS2AB[i])
	for j in list(J):
		for k in range(m):
			J1.extend([m*j+k])
	for d in list(D):
		D1.extend([d]*m)
	M = csr_matrix((D1,(I1,J1)),shape=(p*(m+1)/2,q*m))
	return M


def S2AB1(E,m,DS2AB,ring):
	p = E.nrows()
	q = E.ncols()
	D = E.dict()
	D1 = {}
	for (i,j) in D:
		s = DS2AB[i]
		for k in s:
			s1 = s.index(k)
			D1[(k,m*j+s1)] = D[(i,j)]
	M = matrix(ring,int(p*(m+1)/2),q*m,D1,sparse=True)
	return M

#A is csr_matrix, removes zero rows and zero columns
def trim(A):
	B = A[A.getnnz(1)>0][:,A.getnnz(0)>0]
	return B
'''
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
'''



