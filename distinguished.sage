# Determines distinguished vectors for each weight space (all elementary lowerings)
# Generates ghwv, uses upsets[k] and Distinguished Poset
import itertools

def hwv(wt,a,M):
	X = block_matrix([[i] for i in a[0]])
	K = (X*M[wt]).right_kernel()
	Mwt = M[wt]*matrix(K.basis()[0]).transpose()
	return (wt,Mwt)


def same(q,S):
	B = False
	for s in S:
		b = block_matrix(1,2,[q[1],s[1]])
		if rank(b) == 1:
			B = True
	return B

#Let V = (wt,vector)
#C = Cartan Matrix
def Lowerings(V,a,C):
	LowOps = a[1]
	queue = [V]
	S = []
	while len(queue) != 0:
		q = queue.pop(0)
		if same(q,S) == False:
			S.append(q)
		for i in range(len(LowOps)):
			l = LowOps[i]
			if (l*q[1]).dict() != {}:
				nwt = tuple(a_i - b_i for a_i,b_i in zip(list(q[0]),[C[i,j] for j in range(C.dimensions()[0])]))
				queue.append((nwt,l*q[1]))
	return S

def DistDict(S):
	D = {}
	for s in S:
		if s[0] not in D:
			D[s[0]] = s[1]
		else:
			D[s[0]] = block_matrix(1,2,[D[s[0]],s[1]],subdivide=False)
	return D


def Distinguished(wt,data):
	a = data['a']
	C = data['C']
	M = data['M']
	Vwt = hwv(wt,a,M)
	S = Lowerings(Vwt,a,C)
	return DistDict(S)


#  Generates weight diagram with distinguished vectors
def Poset_Distinguished(wts,data):
	Dicts = []
	for wt in wts:
		Dicts.append(Distinguished(wt,data))
	FDict = {}
	for d in Dicts:
		for k in d.keys():
			if k not in FDict:
				FDict[k] = d[k]
			else:
				FDict[k] = block_matrix(1,2,[FDict[k],d[k]],subdivide=False)
	return FDict


#########################################

def findsubsets(n,k):
	n = list(range(0,n))
	return iter(map(list,itertools.combinations(n,k)))

def SubsetsD(upset,D):
	S = []
	for u in upset:
		n = D[u[0][0]].dimensions()[1]
		k = u[1]
		S.append(findsubsets(n,k))
	S = product(*S)
	return S

#input s from SubsetsD and upsets[k], output ghwv
def dist_hwv(s,upset,D):
	u = upset
	d = [D[u[i][0][0]][:,s[i]] for i in range(len(s))]
	return block_matrix(1,len(s),d)



def Dtot(D):
	S = []
	for d in D.keys():
		S.append((d,D[d].dimensions()[1]))
	return S


#W is dictionary

def Raise(E,V):
	#E = block_matrix(2,1,a[1])
	Rv = [ E*V[:,i] for i in range(V.dimensions()[1])]
	F = [ block_matrix(1,3,[ V[:,range(i)],Rv[i],V[:,range(i+1,len(Rv))]],subdivide=False) for i in range(len(Rv))]
	if all(rank(F[i]) < rank(V) for i in range(len(F))):
		return True
	else:
		return False
	


def grassmannian_hwvs_for_upset_distinguished(mdata,up,verbose):
    C = mdata['C']
    A = mdata['a']
    M = mdata['M'] 
    ssrank = len(A[0])

    if verbose:
        upd = dict(up)
        missing = [p[0] for p in mdata['tot'] if p[-1] > upd.get(p,0)]
        print 'grass hwv missing',missing

    if all((p[-1]-m)*m == 0 for p,m in up):
        # no parameters needed, so its full in every nontrivial coordinate
        # This must be fine, due to conditions on up, so immediately return
        yield block_matrix([[M[p[0]] for p,_ in up]],subdivide=False)
        return
    #print('block_matrix made')

    # Parameters needed. We must consider the product of grassmannians in each
    # weight space with dimensions dictated by up. In the following, letters
    # m and f will refer m dimensional planes inside f dimensional space.

    # Each grassmannian will be covered by (f choose m) charts, corresponding to
    # the choices of m lexicographically smallest pivot rows in a f x m matrix with
    # column space the plane.

    # The following loops over all tuples of all such choices over all the weight spaces
    for nzsi,nzs in enumerate(product(*[combinations(range(p[-1]),m) \
            for p,m in up])):

        if verbose:
            print nzsi
            sys.stdout.flush()
            sys.stdout.write("\033[F")

        nvars = sum( (nzi+1)*(j-i-1) for ((wt,f),m),nz in zip(up,nzs)
            for nzi,(i,j) in enumerate(zip(nz,nz[1:]+(f,))))
        R = PolynomialRing(QQ,'t',nvars,implementation='singular') if nvars > 0 else QQ

        W = {}
        pi = 0
        for q,nz in zip(up,nzs):
            p, m = q
            wt, f = p

            t = matrix(R,f,m,sparse=True)
            t[nz,:] = identity_matrix(R,m,sparse=True)
            for j,ks in enumerate(zip(nz,nz[1:]+(f,))):
                inc = (ks[1]-ks[0]-1)*(j+1)
                t[ks[0]+1:ks[1],:j+1] = matrix(R,ks[1]-ks[0]-1,j+1,R.gens()[pi:pi+inc])
                pi += inc
	    #print(q,nz)
            W[wt] = (M[p[0]]*t, m,f)
	
	cur = block_matrix([[B for B,m,f in W.values()]],subdivide=False)
	
	#print(cur.dimensions())
        # cur will be the result. However, first we must restrict the parameters
        # of cur to those for which the full grassmannian plane is closed under
        # raising operators. In the below, eqs will be populated with the
        # corresponding list of polynomial conditions on the parameters
'''
        eqs = []
        for wt,v in W.items():
            B,m,f = v
            def raisewt(wt,k):
                return tuple(a+b for a,b in zip(wt,C[:,k].list()))
            for k in range(ssrank):
                rwt = raisewt(wt,k)
                curd = W.get(rwt,None)
                if curd is not None:
                    Braise,mr,fr = curd
                    if mr == fr: continue
                    mm = Braise.augment(A[0][k]*B)
                    eqs.extend(minors_sparse(mm,mr+1))
                elif rwt in M:
                    eqs.extend((A[0][k]*B).coefficients())

        if len(eqs) > 0:
            I = R.ideal(eqs)
            if R.one() in I: continue
            Rbar = R.quo(I)
            cur = cur.apply_map(Rbar,sparse=True,R=Rbar)
        #print(cur.dimensions())
	#print(cur.base_ring().gens())
	yield cur
	#if cur.base_ring() != QQ:
        #	for q in Vari(cur):
	#    		yield q
	#else:
	#	yield cur
    if verbose: print
'''

#Evaluate parameters 
def Vari(A):
	S = A.base_ring()
	for i in range(len(S.gens())):
		Dc = {}
		for j in range(len(S.gens())):
			if j != i:
				Dc[S.gens()[j]] = 0
			else:
				Dc[S.gens()[j]] = 1
		yield (A.subs(Dc)).change_ring(QQ)


#Convert upsets[k] to distinguished values
def upsetsD(upset,D):
	S = []
	for p in upset:
		p0,p1 = p
		S.append( (p0[0],D[p0[0]].dimensions()[1],p1))
	return S



def hwvs_for_dist_upsets(D,upset,data):
	M = data['M']
	c = []
	for i in range(len(upset)):
		q0,q1 = upset[i]
		if q1 == q0[-1]:
			c.append((q0[-1],True))
		else:
			c.append((D[q0[0]].dimensions()[1],False))	
	for p in product(*[combinations(range(c[i][0]),upset[i][-1]) for i in range(len(upset))]):
		W = {}
		for i in range(len(upset)):
			if c[i][1] == False:
				W[upset[i][0][0]] = (D[upset[i][0][0]][:,p[i]],upset[i][1],D[upset[i][0][0]].dimensions()[1])
			else:
				W[upset[i][0][0]] = (M[upset[i][0][0]],upset[i][1],D[upset[i][0][0]].dimensions()[1])
		cur = block_matrix([[ B for B,m,f in W.values()]],subdivide=False)
		if rank(cur) == cur.dimensions()[1]:
			yield cur

			

def EQS(W,data):
		M = data['M']
		A = data['a']
		C = data['C']
		ssrank = len(A[0])

	        eqs = []
	        for wt,v in W.items():
		    print('wt = ',wt)
	            B,m,f = v
	            def raisewt(wt,k):
	                return tuple(a+b for a,b in zip(wt,C[:,k].list()))
	            for k in range(ssrank):
	                rwt = raisewt(wt,k)
			print('rwt = ',rwt)
	                curd = W.get(rwt,None)
			print('curd = ',curd)
	                if curd is not None:
	                    Braise,mr,fr = curd
	                    #if mr == fr: continue
	                    mm = Braise.augment(A[0][k]*B)
	                    eqs.extend(minors_sparse(mm,mr+1))
	                elif rwt in M:
			    print('coeff?')
			    print((A[0][k]*B).coefficients())
	                    eqs.extend((A[0][k]*B).coefficients())
		    print('+++')
		
		return eqs

'''	
	        if len(eqs) > 0:
	            I = R.ideal(eqs)
	            if R.one() in I: continue
	            Rbar = R.quo(I)
	            cur = cur.apply_map(Rbar,sparse=True,R=Rbar)
	        #print(cur.dimensions())
	        #print(cur.base_ring().gens())
'''
		#yield cur
		



