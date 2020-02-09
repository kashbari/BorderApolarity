def border_apolarity(T,r,ssonly=True):
    if ssonly:
        Lg,reps = stabilizer_reps_ss(T)
    else:
        Lg,reps = stabilizer_reps(T)
    return border_apolarity_gen(T,Lg,reps,r)

def border_apolarity_gen(T,Lg,reps,r):
    dims = (len(T),) + T[0].dimensions()
    global cands
    cands = [[] for i in range(3)]
    for missing in range(3):
        mdata,em = border_apolarity_110data(T,Lg,reps,missing)
        for ghwv in grassmannian_hwvs(mdata, em.dimensions()[0]-r):
            cand = em*ghwv
            cand = AB_grass_restrict_ok(cand,dims[(missing+1)%3],r,fast=False)
            if cand is not None:
                cands[missing].append(cand)
                print 'candidates',len(cands[missing])
    cand111 = []
    for xs in product(*map(enumerate,cands)):
        ixs = tuple(i for i,x in xs)
        print ixs,
        sys.stdout.flush()
        xs = tuple(x for _,x in xs)
        Rf,Rems = adjoin_rings([x.base_ring() for x in xs])
        W = matrix_to_111(*[x.apply_map(em,Rf) for em,x in zip(Rems,xs)])
        eqs = matrix_rank_le_eqs(W,W.dimensions()[0]-r)
        if 1 in eqs: continue
        print 'candidate'
        cand111.append((W.change_ring(W.base_ring().quo(eqs)),ixs))
    return cand111,cands

# this is special as it exploits cyclic symmetry of Mn to avoid 
# computing I101 and I011. It also uses the natural representation of the
# symmetry group of Mn
def border_apolarity_mmult(n,r):
    T = matrixmult(n,n,n)
    Lg,reps = stabilizer_reps_mmult(n,n,n)
    return border_apolarity_cycl_inv(T,Lg,reps,r)

# border apolarity specialized to cyclic invariant tensors, so only one of I110,
# I101 and I011 are computed. The symmetry lie algebra is also a parameter, so
# special forms of the symmetry group can be used (as in matrix multiplication)
def border_apolarity_cycl_inv(T,Lg,reps,r):
    mdata,em = border_apolarity_110data(T,Lg,reps)
    adim = len(T)
    global cand110
    cand110 = []
    for ghwv in grassmannian_hwvs(mdata, em.dimensions()[0]-r):
        cand = em*ghwv
        cand = AB_grass_restrict_ok(cand,adim,r,fast=False)
        if cand is not None:
            cand110.append(cand)
            print 'candidates',len(cand110)
    cand111 = []
    for xs in product(*map(enumerate,[cand110]*3)):
        ixs = tuple(i for i,x in xs)
        if ixs != min([ixs,cyclicr(ixs),cyclicr(cyclicr(ixs))]):
            continue
        print ixs,
        sys.stdout.flush()
        xs = tuple(x for _,x in xs)
        Rf,Rems = adjoin_rings([x.base_ring() for x in xs])
        W = matrix_to_111(*[x.apply_map(em,Rf) for em,x in zip(Rems,xs)])
        eqs = matrix_rank_le_eqs(W,W.dimensions()[0]-r)
        if 1 in eqs: continue
        print 'candidate'
        cand111.append((W.change_ring(W.base_ring().quo(eqs)),ixs))
    return cand111,cand110

def border_apolarity_110data(T,Lg,reps,missing=2):
    reps = reps[missing+1:] + reps[:missing]
    print 'Afull..'
    Afull = module_dual(module_product(*reps))

    Tperp = T
    for i in range(missing):
        Tperp = tensor_cycl(Tperp)
    Tperp = matrix(QQ,[m.list() for m in Tperp],sparse=True)
    Tperp = Tperp.right_kernel_matrix().transpose().sparse_matrix()

    M = submodule_from_basis(Afull,Tperp)
    data = weight_decomposition(Lg,M)

    return data,Tperp

# # data describes a module, This function returns a generator yielding a set of
# # highest weight vectors of the grassmannian, each possibly parameterized
# # (having entries in a polynomial ring). All highest weight vectors should
# # appear, but perhaps not uniquely.
def grassmannian_hwvs(mdata,subdim,verbose=True):
    for upi,up in enumerate(grassmannian_hwvs_upsets(mdata,subdim)):
        if verbose: print upi,
        for hwt in grassmannian_hwvs_for_upset(mdata,up,verbose):
            yield hwt

def grassmannian_hwvs_FES(mdata,subdim):
    upsets = tuple(grassmannian_hwvs_upsets(mdata,subdim))
    return FESset(upsets).map(lambda up:
            tuple(grassmannian_hwvs_for_upset(mdata,up,False)))

def grassmannian_hwvs_for_upset(mdata,up,verbose=False):
    P = mdata['P']
    C = mdata['C']
    A = mdata['a']
    # M,S = M_to_number_field(mdata['M'])
    M,S = mdata['M'],QQ
    ssdim = len(A[0])

    if verbose: 
        upd = dict(up)
        missing = [p[0]+(p[-1]-upd.get(p,0),) 
                for p in mdata['tot'] if p[-1] > upd.get(p,0)]
        print 'grass hwv missing',missing

    if all((p[-1]-m)*m == 0 for p,m in up):
        yield block_matrix([[M[p[0]] for p,_ in up]],subdivide=False)
    else:
        for j,nzs in enumerate(product(*[combinations(range(p[-1]),m) \
                for p,m in up])):
            if verbose:
                print j,
                sys.stdout.flush()

            names = [x+str(i) for x,q,nz in zip(alphagen(),up,nzs)
                    for i in range(sum((ks[1]-ks[0]-1)*(j+1) 
                        for j,ks in enumerate(zip(nz,nz[1:]+(q[0][-1],))))) ]
            R = PolynomialRing(S,names,implementation='singular') if len(names) > 0 else S

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

                W[wt] = (M[p[0]]*t, m,f)
            cur = block_matrix([[B for B,m,f in W.values()]],subdivide=False)

            eqs = []
            for wt,v in W.items():
                B,m,f = v
                if m==0: continue # for clarity and future changes, should never proc
                def raisewt(wt,k):
                    return tuple(a+b for a,b in zip(wt,C[k,:].list())) + wt[ssdim:]
                for k in range(ssdim):
                    rwt = raisewt(wt,k)
                    curd = W.get(rwt,None)
                    if curd is not None:
                        Braise,mr,fr = curd
                        if mr == fr: continue
                        # ttt = M[(rwt,hwt)]
                        # ttt = ttt.change_ring(R)
                        # ttt = ttt.augment(A[0][k]*B)
                        # if not all(e==0 for e in
                        #         minors_sparse(ttt,ttt.dimensions()[1])):
                        #     import IPython
                        #     IPython.embed()
                        #     assert False
                        mm = Braise.augment(A[0][k]*B)
                        eqs.extend(minors_sparse(mm,mr+1))
                    elif rwt in M:
                        eqs.extend(minors_sparse(A[0][k]*B,1))

            if len(eqs) > 0:
                I = R.ideal(eqs)
                if R.one() in I: continue
                Rbar = R.quo(I)
                cur = cur.apply_map(Rbar,sparse=True,R=Rbar)
            yield cur
        if verbose: print

# Returns an iterator over a set of choices of subspace dimensions for a
# grassmannian highest weight vector which contains at least all legal choices. 

# This is done by deriving inequalities all legal such choices must satisfy
# and enumerating the integer points in the resulting polytope.
def grassmannian_hwvs_upsets(mdata,sz):
    M = mdata['M']
    C = mdata['C']
    X,Y,H = mdata['a']
    tot = mdata['tot']
    wtg = mdata['wtg']

    lp = MixedIntegerLinearProgram()
    for p in tot:
        wt, mult = p
        lp.set_min(lp[p],0)
        lp.set_max(lp[p],mult)

    lp.add_constraint(lp.sum(lp[p] for p in tot) == sz)

    for p in tot:

        for k in range(1,len(wtg.outgoing_edges(p))+1):
            for es in combinations(wtg.outgoing_edges(p),k):
                kdim = block_matrix([[xr] for _,q,xr in es]).right_kernel().dimension()
                lp.add_constraint(lp[p] <= lp.sum(lp[q] for _,q,xr in es) + kdim)

        for k in range(2,len(wtg.incoming_edges(p))+1):
            for es in combinations(wtg.incoming_edges(p),k):
                cokdim = block_matrix([[xr.transpose()] 
                    for q,_,xr in es]).right_kernel().dimension()
                lp.add_constraint(p[-1] - lp[p] <= lp.sum(
                    q[-1] - lp[q] for _,q,xr in es) + kdim)

    for p in tot:
        for _,q,x1 in wtg.outgoing_edges(p):
            for _,w,x2 in wtg.outgoing_edges(q):
                kdim = (x2*x1).right_kernel().dimension()
                lp.add_constraint(lp[p] <= lp[w] + kdim)

    from sage.numerical.mip import MIPSolverException
    def dfs(i):
        if i == len(tot):
            yield [(p,int(lp.get_min(lp[p]))) for p in tot if int(lp.get_min(lp[p])) > 0]
            return
        p = tot[i]
        wt, mult = p
        for val in range(0, mult+1):
            lp.set_min(lp[p],val)
            lp.set_max(lp[p],val)
            try:
                lp.solve()
                for up in dfs(i+1):
                    yield up
            except MIPSolverException:
                pass
        lp.set_min(lp[p],0)
        lp.set_max(lp[p],mult)

    # print len(list(dfs(0)))
    return dfs(0)

def AB_grass_restrict_ok(W,adim,rcandidate,fast=False):
    bdim = W.dimensions()[0] // adim
    M = matrix_11_to_21(W,adim)
    eqs = matrix_rank_le_eqs(M, M.dimensions()[0] - rcandidate,0 if fast else 1)
    if 1 in eqs: return None
    M = matrix_11_to_21(transpose_tensor(W,adim),bdim)
    eqs.extend(matrix_rank_le_eqs(M, M.dimensions()[0] - rcandidate,0 if fast else 1))
    # eqs = list(ideal(eqs).interreduced_basis())
    if 1 in eqs: return None
    if fast: return W
    S = W.base_ring().quo(eqs)
    return W.change_ring(S)

# takes a set of vectors I in A \ot B (in the form of an ab x S matrix), and
# returns a matrix whose columns span A*I subset S^2(A) \ot B
def matrix_11_to_21(B,a):
    b = B.dimensions()[0] // a
    S = B.dimensions()[1]
    W = {}
    # i,k < a
    # j < b
    # s < S
    for I,s in B.nonzero_positions():
        i,j = I // b, I % b
        v = B[I,s]
        for k in range(a):
            mi,ma = min(i,k),max(i,k)
            ix = (( binomial(ma+1,2)+mi )*b+j,s*a+k)
            W[ix] = W.get(ix,0) + v
    W = matrix(B.base_ring(),binomial(a+1,2)*b,S*a,W)
    return W

# A is bc x S, B is ca x U, C is ab x V
def matrix_to_111(A,B,C):
    a = int(sqrt(B.dimensions()[0]*C.dimensions()[0]/A.dimensions()[0]))
    b = C.dimensions()[0] // a
    c = B.dimensions()[0] // a
    W = {}
    for i in range(a):
        for I,l in A.nonzero_positions():
            j,k = I // c, I % c
            W[((i*b+j)*c+k, i*A.dimensions()[1]+l)] = A[j*c+k,l]
    for j in range(b):
        for I,l in B.nonzero_positions():
            k,i = I // a, I % a
            W[((i*b+j)*c+k, a*A.dimensions()[1] + j*B.dimensions()[1]+l)] = B[k*a+i,l]
    for k in range(c):
        for I,l in C.nonzero_positions():
            i,j = I // b, I % b
            W[((i*b+j)*c+k, a*A.dimensions()[1] + b*B.dimensions()[1] +\
                    k*C.dimensions()[1]+l)] = C[i*b+j,l]
    W = matrix(A.base_ring(),a*b*c,a*A.dimensions()[1]+b*B.dimensions()[1]+c*C.dimensions()[1],W)
    return W

# given a set of vectors in A \ot B (in the form of an ab x S matrix), returns
# the same set in B \ot A (in the orm of a ba x S matrix)
def transpose_tensor(B,a):
    b = B.dimensions()[0] // a
    S = B.dimensions()[1]
    Bp = {}
    for I,k in B.nonzero_positions():
        i,j = I // b, I % b
        Bp[(j*a+i,k)] = B[I,k]
    return matrix(B.base_ring(),B.dimensions()[0],B.dimensions()[1],Bp)

# returns an iterator over all order ideals by nondecreasing size
# TODO can speed by topological sorting P and only generating new ones which are
# new, removing the need for checking uniqueness (can use a tree generator)
def order_ideal_iterator(P):
    from sage.combinat.backtrack import TransitiveIdealGraded
    g = P.hasse_diagram()
    mi = P.minimal_elements()
    return TransitiveIdealGraded( lambda down:
            [down.union(frozenset([up]))
                for _,up,_ in g.edge_boundary(down) 
                if all(x in down for x,_,_ in g.incoming_edge_iterator([up]))] +\
            [down.union(frozenset([up])) for up in mi if up not in down]
            ,[frozenset()])

# an iterator over all order filters by nonincreasing size, slightly slower than
# above
def order_filter_iterator(P):
    ps = frozenset(P)
    return (ps-I for I in order_ideal_iterator(P))
    # from sage.combinat.backtrack import TransitiveIdealGraded
    # return TransitiveIdealGraded( lambda up:
    #         [up-frozenset((e,)) for e in P.order_filter_generators(up)]
    #         ,[frozenset(P)])

# converts a weight space dictionary created from weight_decomposition into a common field
def M_to_number_field(M):
    M = [(p,(i,j),v.dimensions(),v[i,j]) 
            for p,v in M.items() for i,j in v.nonzero_positions()]
    S,els,_ = number_field_elements_from_algebraics([v for _,_,_,v in M])
    MM = {}
    for a,v in zip(M,els):
        p,ij,dims,_ = a
        MM.setdefault((p,dims),{})[ij] = v
    M = {k[0]:matrix(S,k[1][0],k[1][1],v) for k,v in MM.items()}
    return M,S

def print_sparse_by_columns(m,st=None,header=None,sort=True):
    if sort: 
        m = m.columns()
        m.sort(key = lambda c: len(c.nonzero_positions()),reverse=True)
        m=matrix(m).transpose()
    if st is None: st = lambda i,j,e: '%s[%d]' % (str(e),i)
    if header is None: header = lambda j: ''
    return '\n'.join(header(j) + ' + '.join(st(i,j,m[i,j]) 
        for i in m.nonzero_positions_in_column(j))
        for j in range(m.dimensions()[1]))

def psbc2(n):
    return lambda i,j,e: '%s a%d%d' % (str(e), (i//n)%n, i%n)

def psbc4(n):
    return lambda i,j,e: '%s a%d%db%d%d' %\
            (str(e), i//n^3, (i//n^2)%n, (i//n)%n, i%n)

def psbc6(n):
    return lambda i,j,e: '%s a%d%db%d%dc%d%d' %\
            (str(e), i//n^5, (i//n^4)%n, (i//n^3)%n, (i//n^2)%n, (i//n)%n, i%n)

def refine_candidates(candidates):
    out = []
    for c in candidates:
        R = c.base_ring()
        try:
            var = R.defining_ideal().variety(QQbar)
            for v in var:
                cc = c.apply_map(lambda e: e.lift()).\
                    subs({ R.cover_ring()(x):e for x,e in v.items() })
                try:
                    cc = cc.change_ring(QQ)
                except:
                    pass
                out.append(cc)
        except AttributeError: # positive dimensional variety
            try:
                for v,m in R.modulus().roots(QQbar):
                    cc = c.apply_map(lambda e:
                        e.lift()).subs({R.cover_ring().gen():v})
                    try:
                        cc = cc.change_ring(QQ)
                    except:
                        pass
                    out.append(cc)
            except AttributeError:
                if R.ngens() == 1 and not R.is_field(): 
                    # how to distinguish a polynomial ring from others?
                    R = PolynomialRing(R.base_ring(),R.gens(),implementation='singular')
                    c = c.change_ring(R)
                out.append(c)
    return out

def print_upset(P,up,verbose=0):
    up = dict(up)
    tot={p:i for i,p in enumerate(reversed(P.linear_extension()))}
    def lab(p):
        if verbose==2:
            return '%s,%s:%d/%d' % (str(p[0]),str(p[1]),up.get(p,0),p[-1])
        elif verbose==1:
            return '%s:%d/%d' % (str(p[0]),up.get(p,0),p[-1])
        else:
            return '%d:%d/%d' % (tot[p],up.get(p,0),p[-1])
    # return DiGraph([(lab(a),lab(b)) for a,b in P.cover_relations()])
    return Poset(([lab(a) for a in P],[(lab(a),lab(b)) for a,b in P.cover_relations()]))

def closed_under_raising(reps, ghwv):
    r = matrix_rank(block_matrix([[ghwv]+[x*ghwv for x in reps[0]+reps[2]]]))
    return r == ghwv.dimensions()[1]

# vim: ft=python
