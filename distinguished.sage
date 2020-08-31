# Determines distinguished vectors for each weight space (all elementary lowerings)
# Generates ghwv, uses upsets[k] and Distinguished Poset
import itertools

def hwv(wt,a,M):
	X = block_matrix([[i] for i in a[0]])
	K = (X*M[wt]).right_kernel()
	Mwt = M[wt]*matrix(K.basis()[0]).transpose()
	return (wt,Mwt)


def same(q,S):
	B = False
	for s in S:
		b = block_matrix(1,2,[q[1],s[1]])
		if rank(b) == 1:
			B = True
	return B

#Let V = (wt,vector)
#C = Cartan Matrix
def Lowerings(V,a,C):
	LowOps = a[1]
	queue = [V]
	S = []
	while len(queue) != 0:
		q = queue.pop(0)
		if same(q,S) == False:
			S.append(q)
		for i in range(len(LowOps)):
			l = LowOps[i]
			if (l*q[1]).dict() != {}:
				nwt = tuple(a_i - b_i for a_i,b_i in zip(list(q[0]),[C[i,j] for j in range(C.dimensions()[0])]))
				queue.append((nwt,l*q[1]))
	return S

def DistDict(S):
	D = {}
	for s in S:
		if s[0] not in D:
			D[s[0]] = s[1]
		else:
			D[s[0]] = block_matrix(1,2,[D[s[0]],s[1]])
	return D


def Distinguished(wt,data):
	a = data['a']
	C = data['C']
	M = data['M']
	Vwt = hwv(wt,a,M)
	S = Lowerings(Vwt,a,C)
	return DistDict(S)


#  Generates weight diagram with distinguished vectors
def Poset_Distinguished(wts,data):
	Dicts = []
	for wt in wts:
		Dicts.append(Distinguished(wt,data))
	FDict = {}
	for d in Dicts:
		for k in d.keys():
			if k not in FDict:
				FDict[k] = d[k]
			else:
				FDict[k] = block_matrix(1,2,[FDict[k],d[k]])
	return FDict




def findsubsets(n,k):
	n = list(range(0,n))
	return iter(map(list,itertools.combinations(n,k)))

def SubsetsD(k,upsets):
	S = []
	for u in upsets[k]:
		n = D[u[0][0]].dimensions()[1]
		k = u[1]
		S.append(findsubsets(n,k))
	S = product(*S)
	return S

def dist_hwv(k,upsets,D):
	S = Subsets(k,upsets)
	
