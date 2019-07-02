import numpy as np
import sympy as sp
import queue
from scipy.sparse import csr_matrix, vstack


n=4

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







