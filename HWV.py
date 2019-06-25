# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 20:16:51 2019

@author: kashb
"""

import numpy as np
import sympy as sp
import queue
from scipy.sparse import csr_matrix



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


def LB2(Q,j,n):
        a = IND(j,n)
        la0 = LB(Q,a[0])
        la1 = LB(Q,a[1])
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


def SparseMat(Q,n):
        row = []
        col = []
        data = []
        for k in range(n**4):
                a = LB2(Q,k,n)
                for j in range(len(a)):
                        if a[j] != None:
                                row.append(a[j])
                                col.append(k)
                                data.append((-1)**j)
        M = csr_matrix((data,(row,col)),shape=(n**4,n**4))
        return M


#V is csr_matrix vector n**4 by 1, mu must be hashable tuple
def Weight(V1):
        (l,l0) = V1.shape
        n = int(l**(1/float(4)))
        a = V1.nonzero()[0]
        if a.size == 0:
                return None
        else:
                c = IND(a[0],n)
                c1 = np.zeros(n)
                for i in range(2):
                        for j in range(2):
                                k = c[i][j]
                                c1[k] += (-1)**j
                FW = np.zeros((n-1,n))
                for i in range(n-1):
                        FW[i,i] = 1
                        FW[i,i+1] = -1
                wt = tuple(map(int,FW.dot(c1)))
                return wt

from scipy.sparse import vstack

#Find Elements of Weight Basis which are correct Weight
def FindWtVecs(wt,n):
	Q = []
	for i in range(n**4):
		row = [i]
		col = [0]
		data = [1]
		V = csr_matrix((data,(row,col)),shape=(n**4,1))
		if wt == Weight(V):
			Q.append(i)
	return Q

def SparseMat1(R,v,n):
	row = []
	col = []
	data = []
	for k in range(len(v)):
		a = LB2(R,v[k],n)
		for j in range(len(a)):
			if a[j] != None:
				row.append(a[j])
				col.append(k)
				data.append((-1)**j)
	M = csr_matrix((data,(row,col)),shape=(n**4,len(v)))
	return M


#compute matrix of raising operators and non-weight vectors being 0
def HWV0(wt,n):
	v = FindWtVecs(wt,n)
	Q = []
	for i in range(n-1):
		Q.append(SparseMat1([i,i+1],v,n))
	Q1 = vstack(Q)
	return Q1

def TorusRestrict(wt,n):
	v = FindWtVecs(wt,n)
	v1 = [IND(k,n) for k in v]
	w0 = []
	for i in range(len(v1)):
		if v1[i][0] == [0,0] or v1[i][1] == [0,0]:
			w0.append(v1[i])
	if len(w0) == 0:
		return None
	else:
		w = []
		for i in range(len(w0)):
			if w0[i][0] == [0,0]:
				k = w0[i][1]
				w1 = [0]*len(v)
				for j in range(n):
					w1[v1.index([[j,j],k])] = 1
				w.append(w1)
		for i in range(len(w0)):
			if w0[i][1] == [0,0]:
				k = w0[i][0]
				w1 = [0]*len(v)
				for j in range(n):
					w1[v1.index([k,[j,j]])] = 1
				w.append(w1)
	w = np.array(w)
	return w	




def HWV(wt,n):
	v = FindWtVecs(wt,n)
	Q = HWV0(wt,n).todense()
	w = TorusRestrict(wt,n)
	Q1 = np.concatenate((Q,w),axis=0)
	Q2 = sp.Matrix(Q1)
	null = Q2.nullspace()
	v1 = [IND(k,n) for k in v]
	return v,v1,null

