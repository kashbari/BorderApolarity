# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 20:16:51 2019

@author: kashb
"""

import numpy as np
import sympy as sp
#from rank_nullspace import rank, nullspace
#from basis_koszul_sln import printSeqUtil, printSeq

print("Let A=B=C=sl_n. Obtain HWVs for weight, w, in sl_n \otimes \Wedge^p sl_n")
n = 3
t = (n**2)-1
t2 = t**2
m = int(n*(n-1)/2)
p=1


w = [2,0,-2]

#y = [0, 0, 1, 2, 0]
#z = np.nonzero(y)
#print(z[0])


''' Translating indexing to coordinates for computation of Lie Bracket '''
def ind(x):
    a = int(np.floor(x/n))
    b= x % n
    c = [a,b]
    return c

def ind2coor1(s):
    s[:] =[ind(x) for x in s];
    f = [x for sublist in s for x in sublist]
    return f

''' Inverse of Prior '''
def iind(x):
    c = n*x[0] +x[1]
    return c


def coor2ind1(x):
        k= int(len(x)/2)
        y=[]
        for i in range(k):
                x1=[x[2*i],x[2*i+1]]
                y.append(iind(x1))
        return y

'''Basis for sln \otimes \Wedge^p sln'''

def printSeqUtil(n,k,len1,arr):
    if (len1 == k):
        A = arr[:]
        LI.append(A)
        return;
    i = 1 if (len1 == 0) else (arr[len1-1]+1);
    len1 += 1;
    while (i <= n):
        arr[len1 - 1] = i;
        printSeqUtil(n,k,len1,arr);
        i += 1;
    len1 -= 1;

def printSeq(n,k):
    arr = [0]*k;
    len1 = 0;
    printSeqUtil(n,k,len1,arr);
    return LI

def BasisWedgeP(n,k):
    S = printSeq(n,k)
    S1 = []
    for s in S:
        s[:] = [x-1 for x in s]
    for s in S:
        for i in range(n):
            s1 = [i] + s
            S1.append(s1)
    return S1

LI = []
S1 = BasisWedgeP(t,p)

S2 = [ind2coor1(s) for s in S1]



#''' Gives Upper triangular indexing as coordinates '''
#def tri(k):
#    i = n - 2 - np.floor(((-8*k + 4*n*(n-1)-7)**(0.5))/2.0 - 0.5)
#    j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2
#    E = [int(i),int(j)]
#    return E

'''Weight Basis '''
def W(v):
    z = [0]*n
    l = len(v)
    for k in range(0,l):
        if (k % 2) == 0:
            z[v[k]] = z[v[k]] +1
        else:
            z[v[k]] = z[v[k]] -1
    return z



def WB(w):
	S = []
	for i in S2:
		v = i
		if W(v) == w:
			S = S + [v]
		else:
			S = S
	return S

for s in WB(w):
	print(s)
	print('\n')





