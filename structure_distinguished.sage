load('borderapolarity3.sage')
load('distinguished.sage')

def e(i,j,n):
    return matrix(QQ,n,n,{(i,j):1})

def basis(i,n):
    assert i < n**2-1
    i,j = i//n, i%n
    if i==j:
        return e(i,i,n) - e(i+1,i+1,n)
    else:
        return e(i,j,n)

# e_ij
def Tsln(n):
    T = [{} for i in range(n**2-1)]
    B = matrix(QQ,[ basis(i,n).list() for i in range(n**2-1)])
    for i in range(n**2-1):
        a = basis(i,n)  
        for j in range(n**2-1):
            b = basis(j,n)
            c = a*b - b*a
            c = B.solve_left(vector(QQ,c.list()))
            for k in c.nonzero_positions():
                T[i][(j,k)] = c[k]
    T = [matrix(QQ,n**2-1,n**2-1,m) for m in T] 

    reps = [[-T[i*n+i+1] for i in range(n-1)],
            [-T[(i+1)*n+i] for i in range(n-1)],
            [-T[i*n+i] for i in range(n-1)]]

    C = matrix(QQ,n-1,n-1)
    for i in range(n-1):
        C[i,i] = 2
    for i in range(n-2):
        C[i+1,i] = -1
        C[i,i+1] = -1

    return T,[reps,reps,module_dual(reps)],C

T,reps,C = Tsln(3)
data,em = border_apolarity_110data(T,reps,C)
r = 16
upsets = list(grassmannian_hwvs_upsets(data,em.dimensions()[0]-r))



D = Poset_Distinguished([(2,2),(3,0),(0,3),(1,1),(0,0)],data)

#print len(upsets)
# for v in grassmannian_hwvs_upsets(data,r):
#     print v
#G = grassmannian_hwvs(data,r)
# for v in grassmannian_hwvs(data,r):
#     print v



#SLURM IT UP
#k = int(sys.argv[1])

#H = list(grassmannian_hwvs_for_upset(data,upsets[k],verbose=True))


def Grassmannian_hwvs(k,mdata,verbose=True):
	for hwt in grassmannian_hwvs_for_upset(data,upsets[k],verbose):
		yield hwt


def border_apolarity_110(T,reps,C,r,k):
	with open("RESULTS/sl3rk14res{}.txt".format(k),'w') as ff:
		mdata,em = border_apolarity_110data(T,reps,C)
		admin = len(T)
		cand110 = []
		i = 0
		ff.write(str(len(H))+'\n')
		for ghwv in Grassmannian_hwvs(k,mdata,em.dimensions()[0]-r):
			cand = em*ghwv
			cand = AB_grass_restrict_ok(cand,admin,r)
			if cand is not None:
				cand110.append(cand)
				ff.write(str(i)+'. Candidate is\n')
			else:
				ff.write(str(i)+'. None\n')
			i = i+1
	return 

#border_apolarity_110(T,reps,C,r,k)

# vim: ft=python
