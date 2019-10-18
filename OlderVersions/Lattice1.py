import numpy as np

n=3
V = np.zeros((9,9))
V[2,2] = 1

def ConvertToMatrix(v,n):
	V1 = np.zeros((n**2,n**2))
	if len(v) == 2:
		V1[v[0],v[1]] += 1
	else:
		for i in range(len(v)):
			if v[i] != None:
				V1[v[i][0],v[i][1]] += (-1)**i
	return V1



