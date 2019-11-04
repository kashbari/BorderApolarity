import BorderApolarity
import PosetHWV

import numpy as np
import sympy as sp
import queue,itertools,math
from scipy.sparse import csr_matrix, vstack, hstack
from collections import deque

import sys
from sage.all import *

import time

n = 3
r  = 14
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
Ca = CartanMatrix(['A',n-1])
Cinv = np.linalg.inv(Ca)
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


Dict,NN = PosetHWV.Combine(Dict22,LE22,N22,[Dict30,Dict03,Dict00,Dict11]) 
LE = LE22
P = P22

Dict30,N30 = PosetHWV.Combine(Dict30,LE30,N30,[Dict00,Dict11])
'''

#precomputed
DIMKER = PosetHWV.DimKer(P,Dict,LE,n)

DIMKER = {}
DIMKER[(0,)] = 35
DIMKER[(1,)] = 35
DIMKER[(0,1)] = 9


H = PosetHWV.dfs(p,LE,NN,P,DIMKER)
'''


def flip(h,LE):
	h1 = [0]*len(h)
	for i in range(len(h)):
		if h[i] != 0:
			j = LE.index(tuple([x*(-1) for x in LE[i]]))
			h1[j] = h[i]
	return h1

load('structure.sage')

def chnge(u):
	u1 = [0]*len(LE)
	for i in range(len(u)):
		j = LE.index(u[i][0][0])
		u1[j] = u[i][1]
	return u1

H = []
for u in upsets:
	H.append(chnge(u)) 

DS2AB = PosetHWV.DictS2AB(n**2-1)


def HwithGrassCharts(H):
	H1 = []
	for h in H:
		K = list(np.subtract(NN,flip(h,LE)))
		A,a,b = PosetHWV.wvs(K,NN,Dict,LE)
		if a == []:
			H1.append((K,None))
		else:
			aa = len(a)
			bb = PosetHWV.Max(b)
			RING = PolynomialRing(QQ,['x%d_%d' %(i,j) for i in range(aa) for j in range(bb)])
			G = []
			for k in range(aa):
				G.extend([PosetHWV.GrassCharts1(b[k][0],b[k][1],RING,k,bb)])
			for g in itertools.product(*G):
				H1.append((K,g,RING))
	return H1

H1 = HwithGrassCharts(H)

def AnnPlane1(h):
	q = H1.index(h)
	with open("sl3rk14res{}.txt".format(q),'w') as ff:
                ff.write('h is'+str(h)+'\n')
                K = h[0]
                A,a,b = PosetHWV.wvs(K,NN,Dict,LE)
		print(A)
                print(a)
                print(b)
		g = h[1]
                if g == None:
                        rk = PosetHWV.wvs01(A,dimS2AB,r,n,DS2AB)
                        if dimS2AB-r >= rk:
                                ff.write('CANDIDATE\n')
                                ff.write(str(K)+'\n')
                        else:
                                ff.write('Nope!\n')
                                ff.write(str(K)+'\n')
		else:
			A = PosetHWV.Convert1(A)
                        aa = len(a)
                        bb = PosetHWV.Max(b)
                        print('aa and bb are')
                        print(aa)
                        print(bb)
                        RING = PolynomialRing(QQ,['x_%d%d' %(i,j) for i in range(aa) for j in range(bb)])
                        print(RING)
			print(g)
			B = A[:]
			for j in range(len(a)):
				print(B[a[j]].nrows(),B[a[j]].ncols())
				B[a[j]] = B[a[j]]*(matrix(g[j]).transpose()).sparse_matrix()
			for j in range(len(A)):
				if B[j] != None:
					B[j] = matrix(RING,B[j])
			B = PosetHWV.shstack(B,RING)
			B = PosetHWV.COB1(B,n)
			print(B.nrows(),B.ncols())      
			t = PosetHWV.wvs1(B,dimS2AB,r,n,RING,DS2AB)     
			if t == True:
				ff.write('CANDIDATE with parameters\n')
				ff.write(str(K)+'\n')
				ff.write(str(g)+'\n')
			else:
				ff.write('Nope! w/ p\n')
				ff.write(str(K)+'\n')
				ff.write(str(g)+'\n')
	return




'''
def AnnPlane0(h):
	q = H.index(h)
	with open("sl3rk14/sl3rk14hwvgr{}.txt".format(q),'w') as ff:
		ff.write('h is'+str(h)+'\n')
		K = list(np.subtract(N,h))
		A,a,b = PosetHWV.wvs(K,N,Dict,LE)
		print(A)
		print(a)
		print(b)
		if a == []:
			rk = PosetHWV.wvs01(A,dimS2AB,r,n,DS2AB)
			if dimS2AB-r >= rk:
				ff.write('CANDIDATE\n')
				ff.write(str(K)+'\n')
			else:
				ff.write('Nope!\n')
				ff.write(str(K)+'\n')
		else:
			A = PosetHWV.Convert1(A)
			aa = len(a)
			bb = PosetHWV.Max(b)
			print('aa and bb are')
			print(aa)
			print(bb)
			RING = PolynomialRing(QQ,['x_%d%d' %(i,j) for i in range(aa) for j in range(bb)])
			print(RING)
			G = []
			for k in range(aa):
				G.extend([PosetHWV.GrassCharts1(b[k][0],b[k][1],RING,k,bb)])
			print('The G are')
			print(G)
			print('The products are')
			for g in itertools.product(*G):
				print(g)
				B = A[:]
				for j in range(len(a)):
					print(B[a[j]].nrows(),B[a[j]].ncols())
					print(g[j])
					B[a[j]] = B[a[j]]*(matrix(g[j]).transpose()).sparse_matrix()
				for j in range(len(A)):
					if B[j] != None:
						B[j] = matrix(RING,B[j])
				B = PosetHWV.shstack(B,RING)
				B = PosetHWV.COB1(B,n)
				print(B.nrows(),B.ncols())	
				t = PosetHWV.wvs1(B,dimS2AB,r,n,RING,DS2AB)	
				if t == True:
					ff.write('CANDIDATE with parameters\n')
					ff.write(str(K)+'\n')
					ff.write(str(g)+'\n')
				else:
					ff.write('Nope! w/ p\n')
					ff.write(str(K)+'\n')
					ff.write(str(g)+'\n')
	return 

# Main Code to run- In Series
for c in C:
	AnnPlane(c)
# Main Code to run- Parallel Computing over 8 processors
import multiprocessing as mp

def main():
	pool = mp.Pool(mp.cpu_count())
	result = pool.map(AnnPlane0,H)

if __name__=="__main__":
	main()
# Main Code to run- SLURM!
AnnPlane1(H1[int(sys.argv[1])])

for h in H[int(sys.argv[1])*350:(1+int(sys.argv[1]))*350]:
	AnnPlane1(h)
'''

