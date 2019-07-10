import numpy as np
import sympy as sp
import queue,itertools,math
import scipy.sparse
from scipy.sparse import csr_matrix,vstack,hstack
import BorderApolarity

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
		l += list(map(list,itertools.combinations(s,i)))
	return l

 

#dimker = {Lowering ops intersections as tuple:dimker}
# K = [k_1, \cdots, k_l], change k_i so check conditions for k_i
#Lin is linear extension, N is maximal elts, P is poset and DIMKER is dimensions of all kernels
def conditions(K,LE,N,P,DIMKER):
	Test = True
	while Test != False:
		if all(N[i] >= K[i] for i in range(len(N))):
			for j in range(len(K)):
				if K[j] != 0:
					parents = P.upper_covers(LE[j])
					if len(parents) != 0:
						pi = [LE.index(p) for p in parents]
						lo = [np.subtract(p,LE[i]).tolist().index(2) for p in parents]
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


