import BorderApolarity
import PosetHWV

n = 3

row = []
col = [0]*len(row)
data = []
V22 = csr_matrix((data,(row,col)),shape=(n**4,1))


row = []
col = [0]*len(row)
data = []
V30 = csr_matrix((data,(row,col)),shape=(n**4,1))


row = []
col = [0]*len(row)
data = []
V03 = csr_matrix((data,(row,col)),shape=(n**4,1))


row = []
col = [0]*len(row)
data = []
V00 = csr_matrix((data,(row,col)),shape=(n**4,1))

row = []
col = [0]*len(row)
data = []
V11 = csr_matrix((data,(row,col)),shape=(n**4,1))




LowOps = [SparseMatLowOp(i,n) for i in range(n-1)]
C = CartanMatrix(['A',n-1])
Cinv = np.linalg.inv(C)
f = lambda x,y: all(i >= 0 for i in map(int,round(Cinv.dot(np.subtract(y,x)))))


def WeightPoset(V,LowOps,f):
	L = LatticeElts(V,LowOps)
	E = L.keys()
	P = Poset((E,f))
	return L,P

L22,P22 = WeightPoset(V22,LowOps,f)
LinExt22 = list(P22.linear_extension())
LinExt22.reverse()
N22 = [L22[s].shape[1] for s in LinExt22]


L30,P30 = WeightPoset(V30,LowOps,f)
LinExt30 = list(P30.linear_extension())
LinExt30.reverse()
N30 = [L30[s].shape[1] for s in LinExt30]


L03,P03 = WeightPoset(V03,LowOps,f)
LinExt03 = list(P03.linear_extension())
LinExt03.reverse()
N03 = [L03[s].shape[1] for s in LinExt03]


L00,P00 = WeightPoset(V00,LowOps,f)
LinExt00 = list(P00.linear_extension())
LinExt00.reverse()
N00 = [L00[s].shape[1] for s in LinExt00]


L11,P11 = WeightPoset(V11,LowOps,f)
LinExt11 = list(P11.linear_extension())
LinExt11.reverse()
N11 = [L11[s].shape[1] for s in LinExt11]








