
# f should be injective
def map_poset(P,f):
    t = list(P)
    tix = {e:i for i,e in enumerate(t)}
    t = [f(e) for e in t]
    return Poset((t,[(t[tix[a]],t[tix[b]])
        for a,b in P.cover_relations()]),cover_relations=True)

def cyclicr(t):
    return (t[2], t[0], t[1])

def alphagen():
    k=1
    while True:
        for cs in product(*['abcdefghijklmnopqrstuvwxyz']*k):
            yield ''.join(cs)
        k+=1

# M a full rank matrix over a laurent series ring F. This function finds and
# applies the transformation G (over F) so that A = G*M satisfies
# 1. A has no nonzero coefficients for a negative power of t
# 2. substituting t=0 into A (well defined operation by 1.) yields a full rank matrix
def make_low_independent(M,trans=False,verbose=False):
    t = M.base_ring().gen()
    M = copy(M)
    if trans:
        M = M.augment(identity_matrix(M.base_ring(),M.dimensions()[0],sparse=True),subdivide=True)
    while True:
        ds = [min([e.valuation() for e in r if e!=0]) for r in M.subdivision(0,0)]
        if verbose: print ds
        mlo = matrix([r.apply_map(lambda e: e[d]) for d,r in zip(ds,M.subdivision(0,0))])
        if mlo.rank() == M.dimensions()[0]: break
        rel = mlo.left_kernel().basis_matrix()[0]
        i = next(i for i in range(M.dimensions()[0]) if rel[i]!=0)
        assert rel[i] != 0
        M[i,:] *= rel[i]
        for j in range(i+1,M.dimensions()[0]):
            M.add_multiple_of_row(i,j, rel[j]*t^(ds[i]-ds[j]) )
        assert min([e.valuation() for e in M.subdivision(0,0)[i].list() if e!=0]) > ds[i]
    for i,r in enumerate(M.subdivision(0,0)):
        M[i] *= t^(-min([e.valuation() for e in r if e!=0]))
    if trans:
        return M.subdivision(0,0),M.subdivision(0,1)
    else:
        return M

# same as above but laurent series ring is over inexact ring
def make_low_independent_approx(M,trans=False,verbose=False,geps=None):
    t = M.base_ring().gen()
    M = copy(M)
    if trans:
        M = M.augment(identity_matrix(M.base_ring(),M.dimensions()[0],sparse=True),subdivide=True)
    while True:
        for i,r in enumerate(M.subdivision(0,0)):
            M[i] *= t^(-min([e.valuation() for e in r if e!=0]))

        mlo = np.array(M.subdivision(0,0).apply_map(lambda e: e[0]))
        u,s,vh = scipy.linalg.svd(mlo,overwrite_a=True)
        eps = geps or s.max() * max(mlo.shape) * np.finfo(mlo.dtype).eps
        print 'eps', eps
        print 'singular values'
        print s
        if s[-1] > eps: break

        M = matrix(CDF,u.conj().T) * M

        for i,j in product(*map(range,M.dimensions())):
            M[i,j] = sum(v*t^k for k,v in M[i,j].laurent_polynomial().dict().items() if abs(v) > eps)

        if verbose: 
            print [min([e.valuation() for e in r if e!=0]) for r in M.subdivision(0,0)]
    if trans:
        return M.subdivision(0,0),M.subdivision(0,1)
    else:
        return M


def getvec(fname='gen/out.txt',cx=False):
    cs = map(float,open(fname).read().split())
    if cx: cs = [CDF(*cs[2*i:2*(i+1)]) for i in range(len(cs)//2)]
    return cs

def putvec(cs,fname='gen/start.txt',cx=False):
    if cx: cs = [e for c in cs for e in [CDF(c).real(),CDF(c).imag()]]
    open(fname,'w').write('\n'.join(map(str,cs))+'\n')

def combinations_multi(l,ks):
    def combinations_multi_int(l,ks):
        if len(ks) == 0:
            yield ()
            return
        for a in combinations(l,ks[0]):
            ca = list(set(l) - set(a))
            for b in combinations_multi_int(ca,ks[1:]):
                yield (a,)+b
    return ([[l[j] for j in js] for js in jss] for jss in
            combinations_multi_int(range(len(l)),ks))

def restrict_to_vars(eqs):
    xs = list(reversed(sorted(ideal(eqs).basis.variables())))
    R = PolynomialRing(eqs[0].parent().base_ring(),xs)
    xsma = {str(x):x for x in R.gens()}
    sub = [xsma.get(str(x),R.zero()) for x in eqs[0].parent().gens()]
    return [e(sub) for e in eqs]


# basis is monomials in order (0,0), (0,1), (1,1), (0,2), (1,2)...
# ie, lexicographic order where the left entries are least significant.
# This is the convenient order to compute rank with respect to
def combinations_with_replacement_alt(n,k):
    assert n > 0
    ix = [0]*k
    yield tuple(ix)
    while True:
        for s in range(k):
            if ix[s] < (ix[s+1] if s+1 < k else n-1):
                ix[s] += 1
                for l in range(s):
                    ix[l] = 0
                break
        else:
            return
        yield tuple(ix)

def combinations_with_replacement_alt_rank(n,k,ix):
    return sum(binomial(e+s,s+1) for s,e in enumerate(ix))

# returns the matrix defining the map S^n V x S^m V -> S^(n+m) V, 
# bases ordered as above
def poly_mult_map(n,m,v):
    M = {}
    bdim = binomial(v+m-1,m)
    rank = lambda kx: combinations_with_replacement_alt_rank(v,n+m,kx)
    for i,ix in enumerate(combinations_with_replacement_alt(v,n)):
        for j,jx in enumerate(combinations_with_replacement_alt(v,m)):
            # print i,bdim,j,i*bdim+j
            M[(rank(sorted(ix+jx)),i*bdim+j)] = 1
    # print (binomial(v+n+m-1,n+m),binomial(v+n-1,n)*bdim)
    # print M
    return matrix(QQ,binomial(v+n+m-1,n+m),binomial(v+n-1,n)*bdim,M)

def matrix_rank(M):
    if M.base_ring() not in [QQ,ZZ]:
        return M.rank()
    M = (M*M.denominator()).change_ring(ZZ)
    if M.is_sparse():
        try:
            return M.rank('linbox')
        except:
            pass
        try:
            return sparse_matrix_rank(M)
        except OSError:
            print 'WARNING: cannot call sparse_matrix_rank. Is cpp/sparse_matrix_rank missing?'
    return M.dense_matrix().rank()

# inexact seems to be faster for smaller matrices (p <= 2, or 100x100),
# but much slower for larger (300x300) and has the slight risk of being
# inaccurate
def inexactrank(M):
    # is there a less verbose (functional) way to do this?
    rank=0
    for sigma in M.change_ring(CDF).dense_matrix().singular_values(eps=1e-13):
        if sigma != 0:
            rank += 1
        else:
            return rank
    return rank

# children before root dfs of sage graph
def graph_dfs(g,rs):
    seen = set()
    def dfs(v):
        if v in seen: return []
        seen.add(v)
        return chain.from_iterable([dfs(w) for _,w,_ in g.outgoing_edges(v)]+[[v]])
    return chain.from_iterable(dfs(v) for v in rs)

# tries to find a topological search order which is nice in the sense that 
# the suffixes should define subgraphs with as many edges as possible.
# Using the reverse of such an order is meant to be useful to backtrack in
# problems where edges represent constraints. (Eg, border apolarity)
def nice_topological_order(g):
    g=copy(g)
    def size_subtree_missing(v,missing):
        seen = copy(missing)
        def dfs(v):
            if v in seen: return 0
            seen.add(v)
            return 1 + sum(dfs(w) for _,w,_ in g.outgoing_edges(v))
        return dfs(v)
    def graph_dfs_sort(g,s):
        seen = set()
        def dfs(v):
            if v in seen: return []
            seen.add(v)
            edges = g.outgoing_edges(v)
            edges.sort(key=lambda e: -size_subtree_missing(e[1],seen))
            return chain.from_iterable([dfs(w) for _,w,_ in edges]+[[v]])
        return dfs(s)
    vs = []
    while len(g.vertices()) > 0:
        s = min([v for v in g if g.in_degree(v) == 0],
            key=lambda s: len(list(g.depth_first_search(s))))
        vsa = list(graph_dfs_sort(g,s))
        for w in vsa:
            g.delete_vertex(w)
        vs.extend(vsa)
    vs.reverse()
    return vs

def nice_order_pre(g, vs):
    vso = nice_topological_order(g.subgraph(vs))[::-1]
    vo = []
    # import IPython
    # IPython.embed()
    for i,v in enumerate(vso):
        vo.append(v)
        # bad = set(g.reverse().depth_first_search(vso[i+1:]))
        bad = set(g.depth_first_search(vso[i+1:],neighbors=g.neighbor_in_iterator))
        bad = bad.union(set(vo))
        vo.extend(nice_topological_order(g.subgraph([v for v in g if v not in bad])))
    assert len(vo) == len(g.vertices())
    return vo

# tensor product of rings, with embeddings
# works with multivariate polynomial rings over QQ and their quotients, and QQ
# dont know how to test for univariate rings
def adjoin_rings(Rs):
    Rquos = [(Ri,R) for Ri,R in enumerate(Rs) \
            if sage.rings.quotient_ring.is_QuotientRing(R)]
    # changed base to generix in 234 and 236
    from sage.rings.polynomial.multi_polynomial_ring_generic import MPolynomialRing_generic
    Rpolys = [(Ri,R) for Ri,R in enumerate(Rs) \
            if isinstance(R,MPolynomialRing_generic)]
    if len(Rquos) == 0 and len(Rpolys) == 0:
        return QQ,[lambda x:x for i in range(len(Rs))]


    startis = [0]+[R.ngens() for _,R in Rquos+Rpolys]
    for i in range(1,len(startis)): startis[i] += startis[i-1]

    # currently assume QQ base ring, in future combine number fields?
    Rout = PolynomialRing(QQ,'x',startis[-1])

    eqs = []
    for i,R in enumerate(Rquos):
        _,R = R
        eqs.extend([p(Rout.gens()[startis[i]:startis[i+1]]) 
            for p in R.defining_ideal().gens()])

    if len(eqs) > 0:
        Rout = Rout.quo(Rout.ideal(eqs))

    Rems = [lambda e: Rout(e) for i in range(len(Rs))]
    for i,R in enumerate(Rquos):
        Ri,_ = R
        Rems[Ri] = lambda p,i=i: p.lift()(Rout.gens()[startis[i]:startis[i+1]])
    sh = len(Rquos)
    for i,R in enumerate(Rpolys):
        Ri,_ = R
        Rems[Ri] = lambda p,i=i: p(Rout.gens()[startis[sh+i]:startis[sh+i+1]])

    return Rout,Rems

# trans*m[:,jxs] = M
# res[:r,:r] is upper triangular, and res[r:,r:] contains no legal pivots
def sparse_elimination(M,pivot_legal=None,pivot_badness=None,
        transformation=False,in_place=False):
    if not in_place:
        M = copy(M)
    if pivot_legal is None:
        pivot_legal = lambda e: e.is_unit()
    jxs = range(M.dimensions()[1])
    if transformation:
        trans = identity_matrix(M.base_ring(), M.dimensions()[0], sparse=True)
    r = 0
    while True:
        if r == min(*M.dimensions()): break
        try:
            if pivot_badness is None:
                i,j = next((i+r,j+r) for (i,j),e in M[r:,r:].dict().items() if
                        pivot_legal(e))
            else:
                i,j = min(((i+r,j+r,e) for (i,j),e in M[r:,r:].dict().items() if
                    pivot_legal(e)),key = lambda p: pivot_badness(p[-1]))[:2]
        except ValueError:
            break
        except StopIteration:
            break

        M.swap_rows(r,i)
        M.swap_columns(r,j)
        jxs[r], jxs[j] = jxs[j], jxs[r]
        if transformation:
            trans.swap_rows(r,i)
            trans[r,:] *= M[r,r].inverse_of_unit()
        M[r,:] *= M[r,r].inverse_of_unit()
        for i in M.column(r)[r+1:].nonzero_positions():
            i += r+1
            if transformation:
                trans.add_multiple_of_row(i,r,-M[i,r])
            M.add_multiple_of_row(i,r,-M[i,r])
        r += 1

    if transformation:
        return M,jxs,r,trans
    else:
        return M,jxs,r

# assumes E in echelon form
# reduce rows of m modulo the rows of E, resulting rows will be the unique
# representative with zeros in the pivot rows of E
def echelon_reduce(E, m):
    for i,j in enumerate(E.pivots()):
        for k in m.nonzero_positions_in_column(j):
            m[k] -= m[k,j]*E[i]

# vim: ft=python
