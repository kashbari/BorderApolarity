# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 20:16:51 2019

@author: kashb
"""

import numpy as np
import sympy as sp
from scipy.sparse import csr_matrix

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

[c for c in findComp(5,3)]	

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


def SparseMatrixForLoweringOperator(i,n):
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




