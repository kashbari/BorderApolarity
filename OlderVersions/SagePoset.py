N = 3 
 
def Partition_Poset(X):
    P = Poset((SetPartitions(X),lambda q,p: q in p.refinements()))
    return P

def p_label(p):
    out = ""
    for block in p:
        for elm in block:
            out += str(elm)
        out += "|"
    return out[:-1]

Po = Partition_Poset(N)
Po.plot(element_labels = {x:p_label(x) for x in Po},vertex_size=500,vertex_shape=None)
