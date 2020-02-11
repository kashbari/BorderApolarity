

# This file was *autogenerated* from the file structure3_15.sage
from sage.all_cmdline import *   # import sage library

_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_4 = Integer(4); _sage_const_114 = Integer(114); _sage_const_15 = Integer(15); _sage_const_21 = Integer(21)
load('borderapolarity3.sage')

def e(i,j,n):
    return matrix(QQ,n,n,{(i,j):_sage_const_1 })

def basis(i,n):
    assert i < n**_sage_const_2 -_sage_const_1 
    i,j = i//n, i%n
    if i==j:
        return e(i,i,n) - e(i+_sage_const_1 ,i+_sage_const_1 ,n)
    else:
        return e(i,j,n)

# e_ij
def Tsln(n):
    T = [{} for i in range(n**_sage_const_2 -_sage_const_1 )]
    B = matrix(QQ,[ basis(i,n).list() for i in range(n**_sage_const_2 -_sage_const_1 )])
    for i in range(n**_sage_const_2 -_sage_const_1 ):
        a = basis(i,n)  
        for j in range(n**_sage_const_2 -_sage_const_1 ):
            b = basis(j,n)
            c = a*b - b*a
            c = B.solve_left(vector(QQ,c.list()))
            for k in c.nonzero_positions():
                T[i][(j,k)] = c[k]
    T = [matrix(QQ,n**_sage_const_2 -_sage_const_1 ,n**_sage_const_2 -_sage_const_1 ,m) for m in T] 

    reps = [[-T[i*n+i+_sage_const_1 ] for i in range(n-_sage_const_1 )],
            [-T[(i+_sage_const_1 )*n+i] for i in range(n-_sage_const_1 )],
            [-T[i*n+i] for i in range(n-_sage_const_1 )]]

    C = matrix(QQ,n-_sage_const_1 ,n-_sage_const_1 )
    for i in range(n-_sage_const_1 ):
        C[i,i] = _sage_const_2 
    for i in range(n-_sage_const_2 ):
        C[i+_sage_const_1 ,i] = -_sage_const_1 
        C[i,i+_sage_const_1 ] = -_sage_const_1 

    return T,[reps,reps,module_dual(reps)],C

T,reps,C = Tsln(_sage_const_3 )
data,em = border_apolarity_110data(T,reps,C)
r = _sage_const_15 
upsets = list(grassmannian_hwvs_upsets(data,em.dimensions()[_sage_const_0 ]-r))




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
		i = _sage_const_0 
		ff.write(str(len(H))+'\n')
		if len(H) != _sage_const_0 :
			G = Grassmannian_hwvs(k,mdata,em.dimensions()[_sage_const_0 ]-r)
			for ghwv in G:
				cand = em*ghwv
				cand = AB_grass_restrict_ok(cand,admin,r)
				if cand is not None:
					cand110.append(cand)
					ff.write(str(i)+'. Candidate\n')
				else:
					ff.write(str(i)+'. None\n')
				i = i+_sage_const_1 
	return

# border_apolarity_110(T,reps,C,r,k)

######### General for 110 and 111 simultaneously
def border_apolarity_110N(T,reps,C,r,S):
	mdata,em = border_apolarity_110data(T,reps,C)
	admin = len(T)
	cand110 = []
	for k in S:
		print(k)
                G = Grassmannian_hwvs(k,mdata,em.dimensions()[_sage_const_0 ]-r)
		if k != _sage_const_4 :
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


S = {_sage_const_4 ,_sage_const_21 ,_sage_const_114 }
cand110 = border_apolarity_110N(T,reps,C,r,S)
print('The number of 110 candidates is:\n')
print(len(cand110))



load('borderapolarity.sage')
load('misc.sage')
def border_apolarity_111N(cands):
	cand111 = []
	for xs in product(*map(enumerate,cands)):
		ixs = tuple(i for i,x in xs)
		print ixs,
		sys.stdout.flush()
		xs = tuple(x for _,x in xs)
		Rf,Rems = adjoin_rings([x.base_ring() for x in xs])
		W = matrix_to_111(*[x.apply_map(em,Rf) for em,x in zip(Rems,xs)])
		eqs = matrix_rank_le_eqs(W,W.dimensions()[_sage_const_0 ]-r)
		if _sage_const_1  in eqs: continue
		print 'candidate'
		cand111.append((W.change_ring(W.base_ring().quo(eqs)),ixs))
	return cand111

cand111 = border_apolarity_111N(cand110)
print('The number of 111 candidate is:\n')
print(len(cand111))
# vim: ft=python

