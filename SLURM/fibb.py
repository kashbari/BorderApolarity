def fibb(n):
    a,b = 1,1
    for i in range(n-2):
        a,b = b,a+b
    return b

import sys
print fibb(int(sys.argv[1]))
