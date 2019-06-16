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
           
#for s in WB(w):
#    print(s,end='\n')


''' Computes Lie bracket E = [m,n] a = [i,j] '''
def LB(E,a):
    c1 = 0
    c2 = 0
    if a[0] ==E[1]:
        c1 = [E[0],a[1]]
        if a[1] == E[0]:
            c2 = [a[0],E[1]]
        else:
            c2 = 0
    else:
        c1 = 0
        if a[1] == E[0]:
            c2 = [a[0],E[1]]
        else:
            c2 = 0
    c = [c1,c2]
    return c

#def LB1(E,a):
#    [c1,c2] = LB(E,a)
#    if c1 == 0 and c2 == 0:
#        d = print(0)    
#    elif c1 == 0 and c2 != 0:
#        d = print("-",c2, sep="")
#    elif c2 == 0 and c1 != 0:
#        d = print(c1)
#    else:
#        d = print(c1,c2,sep="-")
#    return d


def LBD(E,a):
    if a[0] == a[1]:
        c = LB(E,a)
        a1 = [a[0]+1,a[1]+1]
        d = LB(E,a1)
        c = c+ [d[1],d[0]]
    else:
        c=LB(E,a)
    return c

#E= [2,3]
#a = [2,2]
#b = LBD(E,a)
#print(b)


''' Computes Lie bracket E = [m,n] a = [a1,a2,...] = [i,j,k,l,...], returns + - + - '''
def replace(q,v,k):
    q1 = q
    q1[k] = v[0]
    q1[k+1] = v[1]
    return q1

def LBT(E,a):
    b = []
    l = len(a)
    l2 = int(l/2)
    for k in range(0,l2):
        a1 = [a[2*k],a[2*k+1]]
        L = LBD(E,a1)
        for i in range(0,len(L)):
            if L[i] != 0:
                q = replace(a,L[i],2*k)
                b = b + [q]
            else:
                b = b + [0]
    return b

#def LBT1(E,a):
#    b = LBT(E,a)
#    for i in range(0,len(b)):
#        if b[i] != 0:
#            if i == 0:
#                print(b[0],end="")
#            elif i%2 == 1:
#                print('-',b[i],sep="",end="")
#            else:
#                print('+',b[i],sep="",end="")
#        else:
#            pass
#    return 

''' Nice Output for viewing for Weightvectors '''
#def LieBracketEqn(E,w):
#    for s in WB(w):
#        print(WB(w).index(s),s,sep=".",end=' '*25)
#        LBT1(E,s) 
#        print("\n","="*50)
#    return

#E = [2,3]
#for s in WB(w):
#    print(WB(w).index(s),s,sep=".",end=' '*25)
#    LBT1(E,s)
#    print("\n","="*50)
    
''' Get Equations from Lie Algebra Action '''    


def LBasMat(E,w):
    M=[]
    for a in WB(w):
        M = M+ [LBT(E,a)]
    length = len(sorted(M,key=len, reverse=True)[0])
    for i in range(0,len(M)):
        z = [0]*(length-len(M[i]))
        M[i] = M[i] + z 
    return M


E=[2,3]
M = LBasMat(E,w)

''' Equality of elements in WedgeP A'''
def qsort(L):
	if L  == []:
		return []
	else:
		p = L[0]
		l = qsort([x for x in L[1:] if x < p])
		g = qsort([x for x in L[1:] if x >= p])
		return l + [p] + g

def perm_parity(L):
	exp = -1
	for i in range(0,len(L)):
		for j in range(i+1, len(L)):
			if L[i]<L[j]:
				exp += 1
			else:
				exp = exp
	p = (-1)**exp
	return p

def perm(a,b):
	if qsort(a) == qsort(b):
		c = perm_parity(a)*perm_parity(b)
	else:
		c = 0
	return c


'''Partitions set under equivalence relation'''
def partition(a, equiv):
    partitions = [] # Found partitions
    for e in a: # Loop over each element
        found = False # Note it is not yet part of a know partition
        for p in partitions:
            if equiv(e, p[0]): # Found a partition for it!           
                found = True
                break
        if not found: # Make a new partition for it.
            partitions.append([e])
    return partitions



''' Find unique elements in array '''

def UniqInArray(A): 
    n = len(A)
    m = len(A[0])  
    S = []
    for i in range(n): 
        for j in range(m):
            if A[i][j] not in S:
                S = S +[A[i][j]]
            else:
                S = S
    if 0 in S:
        S.remove(0)
    return S 

def UniqInArray1(A):
	A1 = [item for sublist in partition(UniqInArray(A),perm) for item in sublist]
	return A1

U = UniqInArray1(M)
print(U)


''' Getting index of element in array '''

get_indexes = lambda x, xs: [i for (y, i) in zip(xs, range(len(xs))) if x == y]

def IndInArray(n,A):
    S = []
    for i in range(0,len(A)):  
        b = A[i]
        D = [] 
        for x in get_indexes(n,b):
            D = D + [[i,x]] 
        S = S+D
    return S

def EqnfromInd(n,A):
    z = [0]*len(A)
    for s in IndInArray(n,A):
        z[s[0]] = z[s[0]]+(-1)**s[1]
    return z



def EqnsFromMat(w):
    N = []
    for k in range(0,m):
        E = tri(k)
        M = LBasMat(E,w)
        U = UniqInArray(M)
        for u in U:
            eqn = EqnfromInd(u,M)
            N = N +[eqn]
    return N


#def HWV(w):
#    WW = np.array(WB(w))
#    N = EqnsFromMat(w)
#    N= sp.Matrix(N)
#    null = N.nullspace()
#    return WW,N,null
#
#WW = np.array(WB(w))
#N = EqnsFromMat(w)
#N = sp.Matrix(N)
#null = N.nullspace()
#
#print(WW)
#print('='*25)
#print(null)




''' FINAL OUTPUT '''
#w = list(input("w = "))
#
#WW,N,null = HWV(w)
#
#print(WW)
#print('='*25)
#print(null)
