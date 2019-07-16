import numpy as np
import sympy as sp
import queue,itertools,math
import scipy.sparse
from scipy.sparse import csr_matrix,vstack,hstack
import BorderApolarity
from collections import deque

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
			I1.extend([m*i+k])
	for j in list(J):
		for k in range(m):
			J1.extend([m*j+k])
	for d in list(D):
		D1.extend([d]*m)
	M = csr_matrix((D1,(I1,J1)),shape=(p*m,q*m))
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






