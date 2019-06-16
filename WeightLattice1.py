#from sage.all import *
import numpy as np
import hwv_3 as h

'''Generate Cartan matrix for sln'''
def Cartan(n):
	C = np.zeros((n-1,n-1))
	for i in range(0,n-1):
		C[i,i] = 2
	for i in range(0,n-2):	
		C[i,i+1] = -1
		C[i+1,i] = -1
	return C

mu = list(input("mu = "))
n = int(len(mu)+1)
d = int(input("d = "))

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


'''Generate y_(i_1) \cdots y_(i_n)'''
def Seq(n,k):
	v = range(0,n)
#	v = [i+1 for i in v]
	W = [item for item in allstrings(v,k) if LowOp(item) == 1]
	return W
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
		z = h.LBD(x,y)
		if z[0] == 0:
			z1 = z[1]
		else:
			z1 = z[0]
		return z1

def LBequiv(a,b):
	if len(a) == len(b) and rec(pre(a)) == rec(pre(b)) and h.perm_parity(a) == h.perm_parity(b):
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

'''Given irred rep with highest weight mu, gives all weights spaces in irred rep, depth is at most 2sum(mu)+1 '''
def PosetElts(mu,depth):
	S = []
	n = len(mu)+1
	C = Cartan(n)
	for k in range(0,depth):
		I = Seq(n-1,k)
		for j in I:
			nu = sub(mu,CartanWeight(j,n))
			S += nu		
	return S


def WLPosetElts(mu,depth):
	S = PosetElts(mu,depth)
	S = {a.tostring(): a for a in S}
	S = S.values()
	return S

W = WLPosetElts(mu,d)
W1 = range(len(W))


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



