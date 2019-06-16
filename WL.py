# -*- coding: utf-8 -*-
"""
Created on Tue May 21 04:48:51 2019

@author: kashbari

This code will encode HWV as tensor product tuples [i,j] \in A \otimes B where ind(i) tells us element in A and ind(j) tells us element in B. We create a function to get all lowering operators for the dth level and compute poset from that. 

Perhaps wish to encode vectors as sparse matrices....
"""

import numpy as np
import sympy as sp

'''Translate index to coordinates for computation of Lie Bracket '''
n = 4
d = 2
t = int(n*(n-1)/2)

#Convert index to coordinates i.e. 0 to [0,0] and 1 to [0,1]
def ind(x):
	a = int(np.floor(x/n))
	b = x % n
	c = [a,b]
	return c

def ind2coor1(s):
	s[:] = [ind(x) for x in s];
	f = [x for subslist in s for x in sublist]
	return f

#Inverse of prior
def iind(x):
	c = n*x[0] +x[1]
	return c

def coor2ind1(x):
	k = int(len(x)/2)
	y = []
	for i in range(k):
		x1 = [x[2*i],x[2*i+1]]
		y.append(iind(x1))
	return y

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

def Shorten(L):
	K = int(len(L)/2)
	for k in range(0,K):
		if L[2*(K-k)-1] == None and L[2*(K-k)-2] == None:
			del L[2*(K-k)-1]
			del L[2*(K-k)-2]
	return L

#We have a1 = [i,j] as above
def LB1(E,a1):
	la1 = LB(E,ind(a1[0]))
	la2 = LB(E,ind(a1[1]))
	for i in range(len(la1)):
		if la1[i] == None:
			la1[i] = la1[i]
		else:
			la1[i] = iind(la1[i])
	for i in range(len(la2)):
		if la2[i] == None:
			la2[i] = la2[i]
		else:
			la2[i] = iind(la2[i])
	L = [ [la1[0],a1[1]], [la1[1],a1[1]],[a1[0],la2[0]],[a1[0],la2[1]]]
	for k in range(len(L)):
		if L[k][0] == None or L[k][1] == None:
			L[k] = None
	L = Shorten(L)
	return L

#Apply LB1 for a1 = [[i,j], None, [k,l]] etc.
def LB2(E,a1):
	if len(a1) == 2:
		return LB1(E,a1)
	else:
		a2 = []
		for i in range(len(a1)):
			if a1[i] == None:
				a2 += [[None,None]]
			else:
				a2 += [LB1(E,a1[i])]
		for j in range(len(a2)):
			if (j % 2) == 1 and a2[j] != [None,None]:
				a2[j].reverse() 
		a2 = [item for sublist in a2 for item in sublist]
		a2 = Shorten(a2)
		return a2
#Consruct all Sequences of length d from first n numbers. Corresponds to computing all lowering operators of depth d
def allstrings(a,l):
	c = [[]]
	for i in range(l):
		c = [[x]+y for x in a for y in c]
	return c
#triangular numbers for generating lowering ops
def tri(k,n):
	i = int(n-2-floor(sqrt(-8*k+4*n*(n-1)-7)/2.0-5.0))
	j = int(k+i+1-n*(n-1)/2 + (n-i)*(n-i-1)/2)
	return [i,j]

def Seq(n,d):
	v = range(0,n)
	W = [item for item in allstrings(v,d)]
	return W

def Seq2(n,d):
	S = []
	for i in range(0,d+1):
		S+= Seq(n,i)
	return S


#Get Weight Lattice Elements for Poset, S is set of Lowering ops and V is HWV
def NONEA(s):
	if all(i == None for i in s):
		return True
	else:
		return False

def LoweringOps(s,a):
	a1 = a
	for i in range(len(s)):
		E = [s[i]+1,s[i]]
		a1 = LB2(E,a1)
	a1 = Shorten(a1)
	return a1

def LatElts(S,V):
	S1 = []
	for s in S:
		if s == []:
			S1 +=[V]
		else:
			S1 += [LoweringOps(s,V)]
	return S1

def IND(S,S1):
	I = [S.index(x) for x in S if NONEA(S1[S.index(x)]) == False]
	return I

def LatElts2(n,d,V):
	t = int(n*(n-1)/2)
	S = Seq2(t,d)
	S1 = LatElts(S,V)
	SS = [S[i] for i in IND(S,S1)]
	SS1 = [S1[i] for i in IND(S,S1)]
	return SS,SS1



#Find Weight of Vector given in above form from SS1
#def Weight(V1):
#	val = next((i for i,v in enumerate(V1) if v != None),None)
#	return Wt(V1[val])

SS,SS1 = LatElts2(4,2,[2,2])

def fcn(p,q):
	if len(SS[p]) > len(SS[q]):
		return True
	else:
		return False

Elt = range(len(SS1))

P = Poset((Elt,fcn))

P.order_filter([0])



