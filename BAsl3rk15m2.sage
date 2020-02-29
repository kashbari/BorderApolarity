import BorderApolarity
import PosetHWV
import scipy.sparse 

def proj(A,b):
	tA = A.transpose()
	iA = (tA*A).inverse()
	P = A*iA*tA
	if len((P*b).dict()) == 0:
		return 'No proj'
	else:
		return 'Proj!'


D = {(i,i):1 for i in range(0,8)}
D[(4,0)] = -1
D[(8,4)] = -1

M = matrix(9,8,D)
T = M.tensor_product(M)
Rest = T.transpose()

load('BAsl3rk15m.py')

def csr2matsage(X):
	I,J,K = scipy.sparse.find(X)
	dat = {(i,j):v for i,j,v in zip(I,J,K)}
	return matrix(X.shape[0],X.shape[1],dat,sparse=True)


def Proj(Dict,wt,Rest,cand):
	R = Rest*csr2matsage(Dict[wt])
	return proj(R,cand)
	
