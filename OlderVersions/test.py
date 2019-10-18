import numpy as np
import sympy as sp
from scipy.sparse import csr_matrix

#Create Weight Lattice; meant to run in SAGE

w = [1,1]
n = int(len(w) +1)
#C = CartanMatrix(['A',len(w)])


#Convert index to coordinates i.e. 0 to [0,0] and 1 to [0,1]
def ind(x,k):
        a = int(np.floor(x/k))
        b = x % k
        c = [a,b]
        return c

#Inverse of prior
def iind(x,k):
        c = k*x[0] +x[1]
        return c

#Auxiliary to not count too many sequences (
def maxrep(l):
	if l == []:
		return [-1,-1]
	else:
		n = len(l)
		count = 0
		res = l[0]
		for i in range(n):
			cur_count = 1
			for j in range(i+1,n):
				if (l[i] != l[j]):
					break
				cur_count += 1
			if cur_count > count :
				count = cur_count
				res = l[i]
	return [res,count]

#Consruct all Sequences of length d from first n numbers. Corresponds to computing all lowering operators of depth d
def allstrings(a,l):
        c = [[]]
        for i in range(l):
                c = [[x]+y for x in a for y in c]
        return c
#Gives all sequences of length d from the first n digits (starting at 0), let k = w
def maxrep1(l,k):
	K = maxrep(l)
	if k[K[0]] < K[1]:
		return False
	else:
		return True

def Seq(n,d,k):
        v = range(0,n)
        W = [item for item in allstrings(v,d) if maxrep1(item,k)]
	return W
		
#Gives all sequences of length at most d
def Seq2(n,d,k):
        S = []
        for i in range(0,d+1):
                S+= Seq(n,i,k)
        return S

#Computes New Weight given sequence of lowering operators, C is Cartan Matrix, s is string of lowering ops and w is Highest Weight
def LowerWeight(s,w,C):
	Low = w
	for i in range(len(s)):
		Low = np.subtract(Low,C[s[i]])
	return Low


def PosetElts(d,w,C):
	n = len(w)
	S = Seq2(n,d,w)
	SS = []
	for s in S:
		L = LowerWeight(s,w,C)
		SS.append(L)
	SS = map(tuple,SS)
	return S,SS


def fcn(p,q):
	p1 = np.dot(C,np.subtract(np.array(p),np.array(q)))
	if sum(p1) >= 0:
		return True
	else:
		return False




