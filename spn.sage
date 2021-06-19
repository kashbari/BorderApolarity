#load('borderapolarity3.sage')

def e(i,j,n):
    return matrix(QQ,n,n,{(i,j):1})

def spn_basis(i,n):
    assert i < (2*n**2 + n)
    if i < n**2: #A & -A^t
        i,j = i//n, i%n
        return e(i,j,2*n) - e(n+i,n+j,2*n)
    elif i < n**2 + n*(n+1)/2: #B
        i -= n**2;
        i,j = i//n, (n-i)%n
        if i == j:
            return e(i,i+n,2*n)
        else:
            return e(i,j+n,2*n)-e(j,i+n,2*n)
    elif i < n**2 + n*(n+1): #C
        i -= n**2 + n*(n+1)/2
        i,j = i//n, (n-i)%n
        if i == j:
            return e(i+n,i,2*n)
        else:
            return e(i+n,j,2*n)-e(j+n,i,2*n)

