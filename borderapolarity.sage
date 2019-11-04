# In this file, modules of semisimple lie algebras represented in code in the
# following way:
# Distinguish primitive positive root vectors x_1, ..., x_r. These determine 
# corresponding h_i and y_i so that [x_i,y_i] = h_i, [h_i, x_i] = 2x_i and 
# [h_i,y_i] = -2y_i. 

# Representations will be specified by giving the action of these three lists of
# elements. Explicitly, if V is a representation of L of dimension n, V will be
# described by [ [ X_1, ..., X_r], [Y_1, ..., Y_r], [H_1, ..., H_r] ], where
# the X_i, Y_i and H_i are n x n matrices.

# weights will be described by the list of their values on the H_i, which are
# necessarily integers

# In addition to the representations in the above form, routines of this file
# also require the cartan matrix of L, with rows and columns ordered to
# correspond with the x_i.

from itertools import *

# border apolarity specialized to cyclic invariant tensors, so only 110
# candidates are computed. 
def border_apolarity_cycl_inv(T,reps,C,r):
    mdata,em = border_apolarity_110data(T,reps,C)
    adim = len(T)
    cand110 = []
    for ghwv in grassmannian_hwvs(mdata, em.dimensions()[0]-r):
        cand = em*ghwv
        cand = AB_grass_restrict_ok(cand,adim,r)
        if cand is not None:
            cand110.append(cand)
            print 'candidates',len(cand110)

    # For our purposes, assume that there are no parameters in the results
    # (and there is at most one solution for each choice of charts of the
    # grassmannians). If this condition does not hold, the following will throw
    # an exception
    cand110 = [ cand.change_ring(QQ) for cand in cand110 ]

    cand111 = []
    for xs in product(*map(enumerate,[cand110]*3)):
        ixs = tuple(i for i,x in xs)
        if ixs != min([ixs,ixs[1:]+ixs[:1],ixs[2:]+ixs[:2]]):
            # we may skip triples equal to others we check modulo cyclic permutation
            continue
        print ixs,
        sys.stdout.flush()

        W = matrix_to_111(*[x for _,x in xs])
        if W.rank() <= W.dimensions()[0] - r:
            cand110.append(W)
    return cand111,cand110

# computes the weight_decomposition data of the 110, 101, or 011 space
# (depending on missing), and the embedding Tperp -> A \ot B (for 110, eg).
def border_apolarity_110data(T,reps,C,missing=2):
    reps = reps[missing+1:] + reps[:missing]
    print 'Afull..'
    Afull = module_dual(module_product(*reps))

    Tperp = T
    for i in range(missing):
        Tperp = tensor_cycl(Tperp)
    Tperp = matrix(QQ,[m.list() for m in Tperp],sparse=True)
    Tperp = Tperp.right_kernel_matrix().transpose().sparse_matrix()

    M = submodule_from_basis(Afull,Tperp)
    data = weight_decomposition(M,C)
    em = Tperp

    return data,em

# This function returns a generator yielding the set of highest weight vectors 
# of the grassmannian, each possibly parameterized (having entries in a
# polynomial ring modulo an ideal). The grassmannian planes are represented 
# as the column space of matrices.

# mdata is the information for the module as summarized by weight_decomposition
# subdim determines which grassmannian to consider (the dimension of the returned spaces)
def grassmannian_hwvs(mdata,subdim,verbose=True):
    for upi,up in enumerate(grassmannian_hwvs_upsets(mdata,subdim)):
        if verbose: print upi,
        for hwt in grassmannian_hwvs_for_upset(mdata,up,verbose):
            yield hwt

# For a fixed choice of subspace dimensions for each weight space (given by up),
# this function returns a generator yielding the set of highest weight vectors
# for the module given by mdata
def grassmannian_hwvs_for_upset(mdata,up,verbose=False):
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
            print nzsi,
            sys.stdout.flush()

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

            W[wt] = (M[p[0]]*t, m,f)
        cur = block_matrix([[B for B,m,f in W.values()]],subdivide=False)

        # cur will be the result. However, first we must restrict the parameters
        # of cur to those for which the full grassmannian plane is closed under
        # raising operators. In the below, eqs will be populated with the
        # corresponding list of polynomial conditions on the parameters

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

    print len(list(dfs(0)))
    return dfs(0)

# takes a grassmannian candidate, represented as the column space of the matrix
# W, possibly with parameters, and performs the 210 and 120 tests. If no
# parameter value satisfies the test, None is returned, otherwise, W is returned
# in the quotient ring modulo the required equations on the parameters
def AB_grass_restrict_ok(W,adim,r):
    bdim = W.dimensions()[0] // adim
    M = matrix_11_to_21(W,adim)
    eqs = matrix_rank_le_eqs(M, M.dimensions()[0] - r)
    # if 1 in eqs or (M.base_ring() is not QQ and 1 in M.base_ring().ideal(eqs)): 
    if 1 in eqs:
        return None
    M = matrix_11_to_21(transpose_tensor(W,adim),bdim)
    eqs.extend(matrix_rank_le_eqs(M, M.dimensions()[0] - r))
    if 1 in eqs or (M.base_ring() is not QQ and 1 in M.base_ring().ideal(eqs)): 
        return None
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

def module_product(a,b):
    return [[x.tensor_product(identity_matrix(QQ,y.dimensions()[0],sparse=True)) +
            identity_matrix(QQ,x.dimensions()[0],sparse=True).tensor_product(y)
        for x,y in zip(s1,s2)] for s1,s2 in zip(a,b)]

def module_dual(a):
    return [[-x.transpose() for x in s] for s in a]

# B is a basis of column vectors for the desired submodule
# being a submodule is not checked
def submodule_from_basis(a,B):
    return [[restrict_map(m,B,B) for m in s] for s in a]

# B: n x k matrix
# C: m x l matrix
# m: m x n matrix mapping the column space of B to the column space of C
# returns the l x k matrix representing m in the bases B and C
def restrict_map(m,C,B):
    I = C.pivot_rows()
    return C[I,:].solve_right((m*B)[I,:])

# a : a representation of L
# C : the cartan matrix of the semisimple part of L
#
# computes a dictionary summarizing info of this representation. In particular,
# the keys of the returned dictionary are
# 'M': a dictionary mapping weights to a distinguished basis of the corresponding weight space
# 'C': The Cartan matrix of L
# 'a': the representation
# 'tot': a distinguished total ordering of the weights compatible with P; for
#     convenience weights are given along with the multiplicity of the weight space
# 'wtg': directed graph with a vertex for each weight and an edge for each
#     raising operator. The edges are labelled with the raising operator between the
#     corresponding weight spaces in the distinguished bases of those spaces
# 
def weight_decomposition(a,C):
    tot = simultaneous_eigenspace(a[2])

    totkey = lambda wt: tuple(C.solve_right(vector(QQ,wt)).list())
    tot.sort(key=lambda p: totkey(p[1]))

    M = {wt:V for V,wt in tot}
    tot = [(wt,V.dimensions()[1]) for V,wt in tot]

    wtg = DiGraph()
    for p in tot:
        wt, mult = p
        B = M[wt]
        for i,x in enumerate(a[0]):
            wt2 = tuple((C[:,i] + vector(wt).column()).list())
            if wt2 in M:
                wtg.add_edge(p,(wt2,M[wt2].dimensions()[1]),restrict_map(x,M[wt2],B))

    return { 'M' : M, 'tot' : tot, 'C': C, 'a': a, 'wtg': wtg }

# ms : list of matrices over QQ which commute and are diagonalizable
# 
# computes the simultaneous eigenspaces of ms
def simultaneous_eigenspace(ms):
    n = ms[0].dimensions()[0]
    spaces = [(identity_matrix(n,sparse=True) ,())]
    for m in ms:
        nspaces = []
        for B,wt in spaces:
            for wt0,W in restrict_map(m,B,B).eigenspaces_right():
                W = W.basis_matrix().transpose().sparse_matrix()
                nspaces.append((B*W,wt+(wt0,)))
        spaces = nspaces
    return spaces

# T : tensor in A ot B ot C
#
# returns T permuted to lie in B ot C ot A
def tensor_cycl(T):
    m = len(T)
    n,s = T[0].dimensions()
    S = [zero_matrix(QQ,s,m) for i in range(n)]
    for i,j,k in product(range(m),range(n),range(s)):
        S[j][k,i] = T[i][j,k]
    return S

# M : matrix
# r : natural number 
#
# computes generators of the ideal of r+1 by r+1 minors of M
def matrix_rank_le_eqs(M,r):
    if M.base_ring() is QQ:
        return [1] if M.rank() > r else []
    else:
        return list(minors_ideal(M,r+1).basis)

# M : matrix with base ring a polynomial ring over QQ or a quotient of such by an ideal
# nminors: a parameter determining when to revert to the base case in the
#    recursive algorithm. Output should be independent of its value, but
#    performance may be tuned
#
# returns the radical of the ideal of r by r minors of m
def minors_ideal(M,r):
    from collections import Counter
    Rorig = M.base_ring()

    if M.is_zero() or r > min(*M.dimensions()): 
        return Rorig.ideal()
    M,lo = sparse_elimination_by_units(M)
    if r <= lo:
        return Rorig.ideal(1)
    r -= lo
    M = M[lo:,lo:] 
    if M.is_zero():
        return Rorig.ideal()
    M = M[[i for i in range(M.dimensions()[0]) if not M[i].is_zero()], 
          [j for j in range(M.dimensions()[1]) if not M[:,j].is_zero()]]

    # if desired, one can avoid the complexity added by the recursive 
    # procedure and localization ring and immediately return the result 
    # by computing all the minors directly (which in some cases can be 
    # prohibitively too many). To try this, uncomment the following line
    # return Rorig.ideal(minors_sparse(M,r))

    try:
        I = Rorig.defining_ideal()
        R = Rorig.cover_ring()
    except AttributeError:
        I = Rorig.ideal(0)
        R = Rorig

    # M is assumed to have base ring of type MPolyQuoLoc
    def rec(M,r):
        R = M.base_ring()
        # print R,',',r
        if M.is_zero():
            return R._ideal
        if r == 1:
            return R._ideal + [e._p for e in M.coefficients()]

        e = Counter(M.dict().values()).most_common(1)[0][0]

        R1 = MPolyQuoLoc(R.base(), R._ideal+e._p, R._den, nilpotents=False )
        t = M.change_ring(R1)
        t,lo = sparse_elimination_by_units(t)
        I1 = rec(t[lo:,lo:],r-lo) if lo < r else R.base().ideal(1)

        R2 = MPolyQuoLoc(R.base(), R._ideal, R._den*e._p, nilpotents=False )
        t = M.change_ring(R2)
        t,lo = sparse_elimination_by_units(t)
        I2 = rec(t[lo:,lo:],r-lo) if lo < r else R.base().ideal(1)

        return I1.intersection(I2)

    S = MPolyQuoLoc(R,I,1)
    if R is Rorig:
        return rec(M.change_ring(S),r)
    else:
        return rec(M.apply_map(lambda e:
            e.lift()).change_ring(S),r).change_ring(Rorig)

# M : matrix
# r : natural number
#
# computes a list containing all the nonzero r by r minors of matrix M, 
# optimizing for very sparse M
def minors_sparse(M,r):
    M = M[[i for i in range(M.dimensions()[0]) if not M[i].is_zero()], [j for j
        in range(M.dimensions()[1]) if not M[:,j].is_zero()]]
    a,b = M.dimensions()
    out = []
    for ix in combinations(range(a),r):
        for jx in combinations([j for j in range(b) if not M[ix,j].is_zero()],r):
            out.extend([e for e in M[ix,jx].minors(r) if not e.is_zero()])
    return out

# M : matrix
# 
# row reduces M as far as possible pivoting using unit entries. Column swaps are
# also performed. M is modified in place and the number of pivots r is returned. 
# After this operation, M has the form
# [ T A ]
# [ 0 B ]
# where T is r by r and  upper triangular with units along the diagonal, B 
# contains no units
def sparse_elimination_by_units(M):
    M = copy(M)
    def is_unit(e):
        try:
            return e.is_unit()
        except NotImplementedError:
            # This happens if we have an element of a coordinate ring which isnt constant.
            # There is no harm to be conservative here
            return False
    r = 0
    while True:
        if r == min(*M.dimensions()): break
        try:
            i,j = next((i+r,j+r) for i,j in M[r:,r:].nonzero_positions() if
                    is_unit(M[i+r,j+r]))
        except StopIteration:
            break

        M.swap_rows(r,i)
        M.swap_columns(r,j)
        M[r,:] *= M[r,r].inverse_of_unit()
        for i in M.column(r)[r+1:].nonzero_positions():
            i += r+1
            M.add_multiple_of_row(i,r,-M[i,r])
        r += 1
    return M,r

# implementation of the ring R/I localized at a polynomial den, where R is a
# multivariate polynomial ring, I an ideal, and den in R
class MPolyQuoLoc(CommutativeRing):
    def __init__(self, PolyR, I, den, nilpotents=True):
        Ring.__init__(self, base=PolyR, category=CommutativeRings())
        I = PolyR.ideal(I)
        if not nilpotents: I = I.radical()
        J,k = I.saturation(den)
        if 1 in J:
            raise ArithmeticError("invalid multiplicative set, (%s)^%d in %s" %
                    (str(den),k,I))

        self._ideal = PolyR.ideal(J.groebner_basis())
        self._den = self._ideal.reduce(PolyR(den))
        self._populate_coercion_lists_()

    def _repr_(self):
        return "%s mod %s localized at %s" % (self.base(), self._ideal, self._den)

    def _element_constructor_(self, x):
        if isinstance(x, MPolyQuoLocElement):
            return MPolyQuoLocElement(self,x._p,0)*\
                MPolyQuoLocElement(self,x.parent()._den,0).inverse_of_unit()^x._k
        else:
            return MPolyQuoLocElement(self,self.base()(x),0)

    def _coerce_map_from_(self, S):
        if S is self.base():
            return True
        elif isinstance(S, MPolyQuoLoc):
            return self.base() == S.base() \
                    and all(p in self._ideal for p in S._ideal.gens()) \
                    and S._den.divides(self._den)

class MPolyQuoLocElement(CommutativeRingElement):
    def __init__(self, parent, p, k):
        RingElement.__init__(self, parent)
        if k < 0:
            p *= parent._den^(-k)
            k = 0
        J = parent._ideal + parent._den
        while k > 0:
            try:
                p = parent.base()(singular.lift(J,p).sage()[-1,0])
            except TypeError: # p not in J
                break
            k -= 1
        p = parent._ideal.reduce(p)
        self._p = p
        self._k = k

    def _repr_(self):
        if self._k > 0:
            return "(%s)/(%s)^%d" % (self._p, self.parent()._den, self._k)
        else:
            return str(self._p)
 
    def _add_(left, right):
        k = max(left._k, right._k)
        return MPolyQuoLocElement(left.parent(), 
                left._p*left.parent()._den^(k-right._k) + \
                        right._p*left.parent()._den^(k-left._k), k)
 
    def _sub_(left, right):
        k = max(left._k, right._k)
        return MPolyQuoLocElement(left.parent(), 
            left._p*left.parent()._den^(k-right._k) - \
                    right._p*left.parent()._den^(k-left._k), k)

    def _mul_(left, right):
        return MPolyQuoLocElement(left.parent(), left._p * right._p, 
                left._k + right._k)
 
    # c is guaranteed to be in parent().base()
    def _rmul_(self, c):
        return MPolyQuoLocElement(self.parent(), c * self._p, self._k)
 
    def _lmul_(self, c):
        return MPolyQuoLocElement(self.parent(), self._p * c, self._k)
 
    @cached_method
    def is_unit(self):
        try:
            self.inverse_of_unit()
            return True
        except ArithmeticError:
            return False
 
    @cached_method
    def inverse_of_unit(self):
        J = self.parent().base().ideal(list(self.parent()._ideal.groebner_basis())+[self._p])
        S,k = J.saturation(self.parent()._den)
        if 1 not in S:
            raise ArithmeticError("element is not a unit")
        e = self.parent().base()(singular.lift(J,self.parent()._den^k).sage()[-1,0])
        return MPolyQuoLocElement(self.parent(),e, k-self._k)

    def _richcmp_(self, rhs, op):
        from sage.structure.richcmp import richcmp
        return richcmp((self._p,self._k),(rhs._p,rhs._k),op)

    def __hash__(self):
        return hash((self._p,self._k))

# vim: ft=python
