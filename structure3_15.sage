load('borderapolarity3.sage')

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
r = 15
upsets = list(grassmannian_hwvs_upsets(data,em.dimensions()[0]-r))




##### Try 111 add S = {4,21,114} or set of upsets with candidate hwv
#border_apolarity_cycl_inv(T,reps,C,r)

#print len(upsets)
# for v in grassmannian_hwvs_upsets(data,r):
#     print v
G = grassmannian_hwvs(data,r)
# for v in grassmannian_hwvs(data,r):
#     print v


########## SLURM IT UP
#k = int(sys.argv[1])

#H = list(grassmannian_hwvs_for_upset(data,upsets[131],verbose=True))

def Grassmannian_hwvs(k,mdata,verbose=True):
	for hwt in grassmannian_hwvs_for_upset(data,upsets[k],verbose):
		yield hwt

def border_apolarity_110(T,reps,C,r,k):
	with open("RESULTS3_15/sl3rk15res{}_0.txt".format(k),'w') as ff:
		mdata,em = border_apolarity_110data(T,reps,C)
		admin = len(T)
		cand110 = []
		i = 0
		ff.write(str(len(H))+'\n')
		if len(H) != 0:
			G = Grassmannian_hwvs(k,mdata,em.dimensions()[0]-r)
			for ghwv in G:
				cand = em*ghwv
				cand = AB_grass_restrict_ok(cand,admin,r)
				if cand is not None:
					cand110.append(cand)
					ff.write(str(i)+'. Candidate\n')
				else:
					ff.write(str(i)+'. None\n')
				i = i+1
	return

# border_apolarity_110(T,reps,C,r,k)

######### General for 110 and 111 simultaneously
def border_apolarity_110N(T,reps,C,r,S):
	mdata,em = border_apolarity_110data(T,reps,C)
	admin = len(T)
	cand110 = []
	for k in S:
		print(k)
                G = Grassmannian_hwvs(k,mdata,em.dimensions()[0]-r)
		if k != 4:
                        for ghwv in G:
                                cand = em*ghwv
                                cand = AB_grass_restrict_ok(cand,admin,r)
                                if cand is not None:
                                        cand110.append(cand)
		else:
			next(G)
			for ghwv in G:
				cand = em*ghwv
				cand = AB_grass_restrict_ok(cand,admin,r)
				if cand is not None:
					cand110.append(cand)
        return cand110


S = {4,21,114}
cand110 = border_apolarity_110N(T,reps,C,r,S)
print('The number of 110 candidates is:\n')
print(len(cand110))
load('borderapolarity.sage')
load('misc.sage')

cand110 = refine_candidates(cand110)
#Cand110 = [[t] for t in cand110]
P = product(cand110,repeat=3)

def border_apolarity_111N(P):
	cand111 = []
	for xs in P:
		Rf,Rems = adjoin_rings([x.base_ring() for x in xs])
		Q = [x.apply_map(m,Rf) for m,x in zip(Rems,xs)]
		Q1 = [Q[0],-Q[1],-Q[2]]
		W = matrix_to_111(*Q1)
		eqs = matrix_rank_le_eqs(W,W.dimensions()[0]-r)
		if 1 in eqs: continue
		print 'candidate'
		cand111.append([W.change_ring(W.base_ring().quo(eqs))])
	return cand111

cand111 = border_apolarity_111N(P)
print('The number of 111 candidate is:\n')
print(len(cand111))


# vim: ft=python
