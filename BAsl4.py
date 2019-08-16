import BorderApolarity
import PosetHWV

import numpy as np
import sympy as sp
import queue,itertools,math
from scipy.sparse import csr_matrix, vstack, hstack
from collections import deque

import sys
from sage.all import *

n = 4
r  = input('rank to test:')
p = r - (n**2-1)
dimS2AB = int( (n**2)*(n**2-1)*(n**2-1)/2)

row = [51]
col = [0]*len(row)
data = [1]
V202 = csr_matrix((data,(row,col)),shape=(n**4,1))

row = [35,50]
col = [0]*len(row)
data = [-1,1]
V210 = csr_matrix((data,(row,col)),shape=(n**4,1))

row = [55,115]
col = [0]*len(row)
data = [-1,1]
V012 = csr_matrix((data,(row,col)),shape=(n**4,1))

row = [39,54,99,114]
col = [0]*len(row)
data = [1,-1,-1,1]
V020 = csr_matrix((data,(row,col)),shape=(n**4,1))


row = [0,5,10,15,20,40,60,65,80,85,90,95,105,125,130,150,160,165,170,175,190,195,215,235,240,245,250,255]
col = [0]*len(row)
data = [3,-1,-1,-1,4,4,4,4,-1,3,-1,-1,4,4,4,4,-1,-1,3,-1,4,4,4,4,-1,-1,-1,3]
V000 = csr_matrix((data,(row,col)),shape=(n**4,1))

row = [3,23,43,48,53,58,63,83,113,163,178,243]
col = [0]*len(row)
data = [8,16,16,8,-8,-8,8,-8,16,-8,16,8]
V101 = csr_matrix((data,(row,col)),shape=(n**4,1))



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

Dict202,LE202,P202,N202 = WeightPoset(V202,LowOps,f)

Dict210,LE210,P210,N210 = WeightPoset(V210,LowOps,f)

Dict012,LE012,P012,N012 = WeightPoset(V012,LowOps,f)

Dict020,LE020,P020,N020 = WeightPoset(V020,LowOps,f)

Dict000,LE000,P000,N000 = WeightPoset(V000,LowOps,f)

Dict101,LE101,P101,N101 = WeightPoset(V101,LowOps,f)


#precomputed
DIMKER = {}
DIMKER[(0,)] = 126 
DIMKER[(1,)] = 126
DIMKER[(2,)] = 126
DIMKER[(0,1)] = 42
DIMKER[(0,2)] = 60
DIMKER[(1,2)] = 42
DIMKER[(0,1,2)] = 10

C1 = [len(LE202),len(LE210),len(LE012),len(LE020),len(LE000),len(LE101)]
C = [c for c in BorderApolarity.findComp(p,6) if all( i >= 0 for i in np.subtract(C1,c))]


def HWVGrassmannian(c):
	H202=PosetHWV.dfs(c[0],LE202,N202,P202,DIMKER)
	H210=PosetHWV.dfs(c[1],LE210,N210,P210,DIMKER)
	H012=PosetHWV.dfs(c[2],LE012,N012,P012,DIMKER)
	H020=PosetHWV.dfs(c[3],LE020,N020,P020,DIMKER)
	H000=PosetHWV.dfs(c[4],LE000,N000,P000,DIMKER)
	H101=PosetHWV.dfs(c[5],LE101,N101,P101,DIMKER)
	return H202,H210,H012,H020,H000,H101

DS2AB = PosetHWV.DictS2AB(n**2-1)

def AnnPlane(c):
	H202,H210,H012,H020,H000,H101 = HWVGrassmannian(c)
	if len(H202) == 0:
		H202.append(None)
	if len(H210) == 0:
		H210.append(None)
	if len(H012) == 0:
		H012.append(None)
	if len(H020) == 0:
		H020.append(None)
	if len(H000) == 0:
		H000.append(None)
	if len(H101) == 0:
		H101.append(None)
	c1 = map(str,c)
	c1 = ''.join(c1)
	with open("sl3rk{}.txt".format(c1),'w') as ff:
		ff.write('c is'+str(c)+'\n')
		for z in itertools.product(H202,H210,H012,H020,H000,H101):
			K = []
			if z[0] != None:
				K.append(list(np.subtract(N202,z[0])))
			else:
				K.append(N202)
			if z[1] != None:
				K.append(list(np.subtract(N210,z[1])))
			else:
				K.append(N210)
			if z[2] != None:
				K.append(list(np.subtract(N012,z[2])))
			else:
				K.append(N012)
			if z[3] != None:
				K.append(list(np.subtract(N020,z[3])))
			else:
				K.append(N020)
			if z[4] != None:
				K.append(list(np.subtract(N000,z[3])))
			else:
				K.append(N000)
			if z[5] != None:
				K.append(list(np.subtract(N101,z[4])))
			else:
				K.append(N101)
			A202,a202,b202 = PosetHWV.wvs(K[0],N202,Dict202,LE202)
			A210,a210,b210 = PosetHWV.wvs(K[1],N210,Dict210,LE210)
			A012,a012,b012 = PosetHWV.wvs(K[2],N012,Dict012,LE012)
			A020,a020,b020 = PosetHWV.wvs(K[3],N020,Dict020,LE020)
			A000,a000,b000 = PosetHWV.wvs(K[4],N000,Dict000,LE000)
			A101,a101,b101 = PosetHWV.wvs(K[5],N101,Dict101,LE101)
			if a202 == a210 == a012 == a020 == a000 == a101 == []:
				A = A202 + A210 + A012 + A020 + A000 + A101
				rk = PosetHWV.wvs0(A,dimS2AB,r,n,DS2AB)
				if dimS2AB -r >= rk:
					ff.write('CANDIDATE\n')
					ff.write(str(K)+'\n')
				else:
					ff.write('No Candidate\n')
					ff.write(str(K)+'\n')
			else:
				A202 = PosetHWV.Convert1(A202)
				A210 = PosetHWV.Convert1(A210)
				A012 = PosetHWV.Convert1(A012)
				A020 = PosetHWV.Convert1(A020)
				A000 = PosetHWV.Convert1(A000)
				A101 = PosetHWV.Convert1(A101)
				A = A202 + A210 + A012 + A020 + A000 + A101
				for k in range(len(a210)):
					a210[k] = a210[k] + len(A202)
				for k in range(len(a012)):
					a012[k] = a012[k] + len(A202)+len(A210)
				for k in range(len(a020)):
					a020[k] = a020[k] + len(A202)+len(A210)+len(A012)
				for k in range(len(a000)):
					a000[k] = a000[k] + len(A202)+len(A210)+len(A012)+len(A020)
				for k in range(len(a101)):
					a101[k] = a101[k] + len(A202)+len(A210)+len(A012)+len(A020)+len(A00)
				a = a202+a210+a012+a020+a000+a101
				aa = len(a)
				bb = PosetHWV.Max(b202+ b210+ b012+ b020 +b000 +b101)
				RING = PolynomialRing(QQ,['x_%d%d' %(i,j) for i in range(aa) for j in range(bb)])
				print(RING)
				G202 = []
				G210 = []
				G012 = []
				G020 = []
				G000 = []
				G101 = []
				for k in range(len(a202)):
					G202.extend([PosetHWV.GrassCharts1(b202[k][0],b202[k][1],RING,k,aa,bb)])
				for k in range(len(a210)):
					G210.extend([PosetHWV.GrassCharts1(b210[k][0],b210[k][1],RING,k+len(a202),aa,bb)])
				for k in range(len(a012)):
					G012.extend([PosetHWV.GrassCharts1(b012[k][0],b012[k][1],RING,k+len(a202+a210),aa,bb)])
				for k in range(len(a020)):
					G020.extend([PosetHWV.GrassCharts1(b020[k][0],b020[k][1],RING,k+len(a202+a210+a012),aa,bb)])
				for k in range(len(a000)):
					G000.extend([PosetHWV.GrassCharts1(b000[k][0],b000[k][1],RING,k+len(a202+a210+a012+a020),aa,bb)])
				for k in range(len(a101)):
					G101.extend([PosetHWV.GrassCharts1(b101[k][0],b101[k][1],RING,k+len(a202+a210+a012+a020+a000),aa,bb)])
				G = G202 + G210 + G012 + G020 + G000+ G101
				for g in itertools.product(*G):
					B = A[:]
					for j in range(len(a)):
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
						ff.write('No Candidate with parameters\n')
						ff.write(str(K)+'\n')
						ff.write(str(g)+'\n')
	return 


'''
# Main Code to run- Local Parallel Processing (8 cpu)
import multiprocessing as mp

def main():
        pool = mp.Pool(mp.cpu_count())
        result = pool.map(AnnPlane,C)

if __name__=="__main__":
        main()

# Main Code to run- SLURM
AnnPlane(C[int(sys.argv[1])])
'''
