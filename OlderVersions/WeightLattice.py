#from sage.all import *
import numpy as np

'''Generate Cartan matrix for sln'''
def Cartan(n):
	C = np.zeros((n-1,n-1))
	for i in range(0,n-1):
		C[i,i] = 2
	for i in range(0,n-2):	
		C[i,i+1] = -1
		C[i+1,i] = -1
	return C

#mu = list(input("mu = "))
#n = int(len(mu)+1)
#d = int(input("d = "))

'''Coordinates not same'''
def LowOp(v):
	k = 1
	for i in range(0,len(v)-1):
		if v[i] == v[i+1]:
			k = 0
		else:
			k = k
	return k

'''Generate All strings of length k from alphabet v'''
def allstrings(alphabet, length):
    """Find the list of all strings of 'alphabet' of length 'length'"""

    c = [[]]
    for i in range(length):
        c = [[x]+y for x in alphabet for y in c]
    return c


'''Compute LieBracket LB1 applied [i,j] to [[k,l],[m,n], ..] where it is +-+- of weight basis'''
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

def LB1(E,a):
	if isinstance(a[0],list) == False:
		L = LB(E,a)
	else:
		L = []
		for i in range(0,len(a)):
			if a[i] == None:
				L += [None,None]
			else:
				L += LB(E,a[i])
	for k in range(0,len(L)):
		if k % 4 == 2:
			L[k],L[k+1] = L[k+1],L[k]
	return L
#Applies Sequence,s, of lowering operators to a i.e. [0,1,0] applies [[1,0],[2,1],[1,0]] to a

def Shorten(a1):
	K = int(len(a1)/2)
	for k in range(0,K):
		if a1[2*(K-k)-1] == None and a1[2*(K-k)-2] == None:
			del a1[2*(K-k)-1]
			del a1[2*(K-k)-2]
	return a1

def LoweringOps(s,a):
	a1 = a
	for i in range(0,len(s)):
		E = [s[i]+1,s[i]]
		a1 = LB1(E,a1)
	a1 = Shorten(a1)
	return a1
	
#Computes all Sequences of length k from set of n
def Seq(n,k):
	v = range(0,n)
	#	v = [i+1 for i in v]
	W = [item for item in allstrings(v,k)]
	return W

def Seq2(n,k):
	S = []
	for i in range(0,k+1):
		S += Seq(n,i)
	return S


# Get Weight Basis, S is set of lowering ops, and V is hwvi

def NONEA(s):
	if all(i == None for i in s):
		return True
	else:
		return False

def WeightBasis(S,V):
	S1 = []
	for s in S:
		if s == []:
			S1 += [V]
		else:
			S1 += [LoweringOps(s,V)]
	return S1

def IND(S,S1):
	ind = [S.index(x) for x in S if NONEA(S1[S.index(x)]) == False]
	return ind



def WeightBasis2(n,d,V):
	S = Seq2(n,d)
	S1 = WeightBasis(S,V)
	SS = [S[i] for i in IND(S,S1)]
	SS1 = [S1[i] for i in IND(S,S1)]
	return SS,SS1














#Take list [i_1, \cdots, i_n] and reduces using lie bracket, where i_1 = [i_1+1,i_1]
def pre(s):
	s1 = [0]*len(s)
	for i in range(0,len(s)):
		s1[i] = [s[i]+1,s[i]]
	return s1

def rec(v):
	w = reduce((lambda x,y: LieB(x,y)),v)
	return w

def LieB(x,y):
	if x == 0 or y == 0:
		return 0
	else:
		z = LBD(x,y)
		if z[0] == 0:
			z1 = z[1]
		else:
			z1 = z[0]
		return z1

def perm_parity(L):
	exp = -1
	for i in range(0,len(L)):
		for j in range(i+1,len(L)):
			if L[i]<L[j]:
				exp += 1
			else:
				exp = exp
	p = (-1)**exp
	return p

def LBequiv(a,b):
	if len(a) == len(b) and rec(pre(a)) == rec(pre(b)) and perm_parity(a) == perm_parity(b):
		return True
	else:
		return False

#Takes Sequences and mods by equivalence under Lie bracket lowering operations
def Seq1(n,k):
	S = Seq(n,k)
	classes = []
	for x in S:
		for c in classes:
			if LBequiv(x,c[0]):
				c.append(x)
				break
		else:
			classes.append([x])
	return classes

'''Subtract vectors componentwise'''
def sub(v,w):
	u = [i-j for i,j in zip(v,w)]
	return u

def CartanWeight(v,n):
	C = Cartan(n)
	S = np.zeros((1,n-1))
	for i in range(0,len(v)):
		S += C[v[i]]
	return S

'''Given HWV, v, compute all lowering vectors'''
#Find Weight of vector
def WW(v,n):
	z = [0]*n
	for k in range(0,len(v)):
		if (k % 2) == 0:
			z[v[k]] = z[v[k]] + 1
		else:
			z[v[k]] = z[v[k]] - 1
	return z

def Weight(v,n):
	k = 0
	for i in range(len(v)):
		if v[i] == 0:
			k = k
		else:
			k = np.dot((-1)**i,WW(v[i],n)).tolist()
			break
	return k
	
#Apply Lowering operator [i,j,k,...] to hwv V

''' Compare two weights q,p and return True is q \leq p'''
def f1(v):
	if sum(v) <= 0:
		return False
	return True

'''	for i in range(0,len(v)):
		if v[i] < 0:
			return False
'''
def f2(v):
	n = int(len(v)+1)
	l = len(v)
	C = Cartan(n)
	C1 = np.matmul(C,v)
	return C1

f = lambda q,p: f1(f2(sub(p,q))) == True
ff = lambda q,p: f(W[q],W[p])


'''Plot Hasse Diagram of Poset... Sage Code'''
'''
P = Poset((W1,ff))
P.show()
'''



