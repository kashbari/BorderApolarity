# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 20:16:51 2019

@author: kashb
"""

import numpy as np
import sympy as sp
import queue
from scipy.sparse import csr_matrix, vstack

#Construct all compositions of n into k parts

def findComp(n,k):
	if n<0 or k<0:
		return 
	elif k == 0:
		if n == 0:
			yield []
		return
	elif k == 1:
		yield [n]
		return
	else:
		for i in range(0,n+1):
			for comp in findComp(n-i,k-1):
				yield [i] + comp

#[c for c in findComp(5,3)]	

def part(n, k):
	def _part(n, k, pre):
		if n <= 0:
			return []
		if k == 1:
			if n <= pre:
				return [[n]]
			return []
		ret = []
		for i in range(min(pre, n), 0, -1):
			ret += [[i] + sub for sub in _part(n-i, k-1, i)]
		return ret
	return _part(n, k, n)

def part1(n,k):
	p = []
	for i in range(1,k+1):
		p += part(n,i)
	return p
'''Indexing for Weight Basis'''
def Ind(x,k):
	a = int(np.floor(x/k))
	b = x % k
	return [a,b]

def IInd(x,k):
	return k*x[0] + x[1]

def IND(x,k):
	[a,b] = Ind(x,k**2)
	return [Ind(a,k),Ind(b,k)]

def IIND(x,k):
	a = [IInd(x[0],k),IInd(x[1],k)]
	return IInd(a,k**2)


'''Lie Bracket for E=[m,n] and a1 = [ind(i),ind(j)]=[[i0,i1],[j0,j1]] input E and a1 = [i,j] from Weight Basis'''
#Compute E=[m,n] lie bracket with a = ind(i)= [i0,i1]
def LB(E,a):
        c1 = None
        c2 = None
        if a[0] == E[1]:
                c1 = [E[0],a[1]]
                if a[1] == E[0]:
                        c2 = [a[0],E[1]]
                else:
                        c2 = None
        else:
                c1 = None
                if a[1] == E[0]:
                        c2 = [a[0],E[1]]
                else:
                        c2 = None
        c = [c1,c2]
        return c

#implicitly Lowering for sln, i parametrizes torus basis and j corresponds to indexing of A \otimes B
def LB1(i,j,n):
	a = IND(j,n)
	E = [i+1,i]
	la0 = LB(E,a[0])
	la1 = LB(E,a[1])
	b = []
	if la0[0] == None:
		b.append(None)
	else:
		b.append([la0[0],a[1]])
	if la0[1] == None:
		b.append(None)
	else:
		b.append([la0[1],a[1]])
	if la1[0] == None:
		b.append(None)
	else:
		b.append([a[0],la1[0]])
	if la1[1] == None:
		b.append(None)
	else:
		b.append([a[0],la1[1]])
	for i in range(len(b)):
		if b[i] != None:
			b[i] = IIND(b[i],n)
	return b


def SparseMatLowOp(i,n):
	row = []
	col = []
	data = []
	for k in range(n**4):
		a = LB1(i,k,n)
		for j in range(len(a)):
			if a[j] != None:
				row.append(a[j])
				col.append(k)
				data.append((-1)**j)
	M = csr_matrix((data,(row,col)),shape=(n**4,n**4))
	return M

#LowOps = [SparseMatLowOp(i,n) for i in range(n-1)]

#V is csr_matrix vector n**4 by 1, mu must be hashable tuple
def Weight(V1):
	(l,l0) = V1.shape
	n = int(l**(1/float(4)))
	a = V1.nonzero()[0]
	if a.size == 0:
		return None
	else:
		c = IND(a[0],n)
		c1 = np.zeros(n)
		for i in range(2):
			for j in range(2):
				k = c[i][j]
				c1[k] += (-1)**j
		FW = np.zeros((n-1,n))
		for i in range(n-1):
			FW[i,i] = 1
			FW[i,i+1] = -1
		wt = tuple(map(int,FW.dot(c1)))
		return wt
			

def LatticeElts(V,LowOps):
	q = queue.Queue()
	q.put(V)
	mu = Weight(V)
	Dict = {mu:V}
	while q.empty() == False:
		P = q.get()
		for i in range(len(LowOps)):
			P1 = LowOps[i].dot(P)
			if P1.nnz != 0:
				wt = Weight(P1)
				if wt not in Dict:
					Dict[wt] = P1
					q.put(P1)
				else:
					P0 = Dict[wt].todense()
					P2 = np.concatenate((P0,P1.todense()),axis=1)
					r = np.linalg.matrix_rank(P2)
					if r == P2.shape[1]:
						Dict[wt] = csr_matrix(P2)
						q.put(P1)
	return Dict


					

'''
n = 3
V = np.zeros((n**4,1))
V[IIND([[0,n-1],[0,n-1]],n)] = 1
V = csr_matrix(V)
LowOps = [SparseMatLowOp(i,n) for i in range(n-1)]
L = LatticeElts(V,LowOps)



Sage Codes for Poset

E = L.keys()
C = CartanMatrix(['A',n-1])
Cinv = np.linalg.inv(C)
f = lambda x,y: all(i >= 0 for i in map(int,round(Cinv.dot(np.subtract(y,x)))))



P = Poset((E,f))


Sage Codes for counting all saturated chains

LinExt = list(P.linear_extension())
LinExt.reverse()
N = [L[s].shape[1] for s in LinExt]


K = [0]*len(L)


def dfs_iterate(P,LinExt):
	N = [L[s].shape[1] for s in LinExt]
	K = [0]*len(L)
	stack = [K]
		
'''














'''
Dimensions for Kernel of intersections of lowering operators

K = list of operators considered

def dimker(K,n):
	if len(K) == 0:
		return None
	else:
		M = []
		for k in K:
			M += SparseMatLowOp(k,n)
		M1 = vstack(M)
		M2 = sp.Matrix(M1.todense())
		null = M2.nullspace()
		return len(null)


s = set of all lowering operators
def subsets(s):
	l = []
	for i in range(1,len(s)):
		l += list(map(set,itertools.combinatios(s,i)))
	return l


l = subsets(s)
l1 = [dimker(r,n) for r in subsets(s)]
'''
