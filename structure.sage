load('borderapolarity.sage')

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
r = 50
upsets = list(grassmannian_hwvs_upsets(data,r))
print len(upsets)
# for v in grassmannian_hwvs_upsets(data,r):
#     print v

# for v in grassmannian_hwvs(data,r):
#     print v

# vim: ft=python
