# Quotient of Multivariate Polynomial Ring in t0, t1, t2, t3, t4, t5, t6, t7, t8 over Rational Field by the ideal (-t0*t3 - t2 + 1, -t1*t7 + t1*t8 - t5 - t8 + 1, -t2*t7 + t2*t8 - t7 - t8, -t3*t7 + t3*t8 - 2*t8 - 2, -t4*t7 + t4*t8 - t6 + t8, -t2*t5 + 2*t1*t7 - t2*t7 + t2 + t5 - t7 - 1, -t3*t5 + 2*t1*t7 - t3*t7 + 2*t1 + t3 + 2*t5 - 4, -t4*t5 + t1*t6 - t1*t7 - t4*t7 + t4 - t5 - t6 + 1, 2*t2*t7 - 2*t3*t7 + 2*t2 + 2*t7 - 2, t2*t6 - t2*t7 - 2*t4*t7 - t6 - t7, t3*t6 - t3*t7 - 2*t4*t7 - 2*t4 - 2*t6 - 2, -t2*t5 + 2*t1*t8 - t2*t8 + t2 - t5 - t8 + 1, -t3*t5 + 2*t1*t8 - t3*t8 + 2*t1 + t3, -t4*t5 + t1*t6 - t1*t8 - t4*t8 + t4, 2*t2*t8 - 2*t3*t8 + 2*t2 + 2*t8 + 2, t2*t6 - t2*t8 - 2*t4*t8 + t6 - t8, t3*t6 - t3*t8 - 2*t4*t8 - 2*t4, -2*t2*t5 + 2*t3*t5 - 4*t1 + 4*t2 - 2*t3 - 2*t5 + 4, t2*t5 + 2*t4*t5 - 2*t1*t6 + t2*t6 - t2 - 2*t4 + t5 + t6 - 1, t3*t5 + 2*t4*t5 - 2*t1*t6 + t3*t6 - 2*t1 - t3 - 4*t4, -2*t2*t6 + 2*t3*t6 - 2*t2 - 4*t4 - 2*t6 - 2, t1*t7 - t1*t8 + t5 + t8 - 1, t2*t7 - t2*t8 + t7 + t8, t3*t7 - t3*t8 + 2*t8 + 2, t2*t5 - 2*t1*t7 + t2*t7 - t2 - t5 + t7 + 1, t3*t5 - 2*t1*t7 + t3*t7 - 2*t1 - t3 - 2*t5 + 4, -2*t2*t7 + 2*t3*t7 - 2*t2 - 2*t7 + 2, t2*t5 - 2*t1*t8 + t2*t8 - t2 + t5 + t8 - 1, t3*t5 - 2*t1*t8 + t3*t8 - 2*t1 - t3, -2*t2*t8 + 2*t3*t8 - 2*t2 - 2*t8 - 2, 2*t2*t5 - 2*t3*t5 + 4*t1 - 4*t2 + 2*t3 + 2*t5 - 4, t2*t7 - t2*t8 + t7 + t8, t3*t7 - t3*t8 + 2*t8 + 2, t4*t7 - t4*t8 + t6 - t8, 2*t2*t7 - 2*t3*t7 + 2*t2 + 2*t7 - 2, t2*t6 - t2*t7 - 2*t4*t7 - t6 - t7, t3*t6 - t3*t7 - 2*t4*t7 - 2*t4 - 2*t6 - 2, 2*t2*t8 - 2*t3*t8 + 2*t2 + 2*t8 + 2, t2*t6 - t2*t8 - 2*t4*t8 + t6 - t8, t3*t6 - t3*t8 - 2*t4*t8 - 2*t4, 2*t2*t6 - 2*t3*t6 + 2*t2 + 4*t4 + 2*t6 + 2, t2*t7 - t2*t8 + t7 + t8, t3*t7 - t3*t8 + 2*t8 + 2, 2*t2*t7 - 2*t3*t7 + 2*t2 + 2*t7 - 2, 2*t2*t8 - 2*t3*t8 + 2*t2 + 2*t8 + 2)



R = PolynomialRing(QQ,9,'t')
R.inject_variables()

E = [ -t0*t3 - t2 + 1, -t1*t7 + t1*t8 - t5 - t8 + 1, -t2*t7 + t2*t8 - t7 - t8, -t3*t7 + t3*t8 - 2*t8 - 2, -t4*t7 + t4*t8 - t6 + t8, -t2*t5 + 2*t1*t7 - t2*t7 + t2 + t5 - t7 - 1, -t3*t5 + 2*t1*t7 - t3*t7 + 2*t1 + t3 + 2*t5 - 4, -t4*t5 + t1*t6 - t1*t7 - t4*t7 + t4 - t5 - t6 + 1, 2*t2*t7 - 2*t3*t7 + 2*t2 + 2*t7 - 2, t2*t6 - t2*t7 - 2*t4*t7 - t6 - t7, t3*t6 - t3*t7 - 2*t4*t7 - 2*t4 - 2*t6 - 2, -t2*t5 + 2*t1*t8 - t2*t8 + t2 - t5 - t8 + 1, -t3*t5 + 2*t1*t8 - t3*t8 + 2*t1 + t3, -t4*t5 + t1*t6 - t1*t8 - t4*t8 + t4, 2*t2*t8 - 2*t3*t8 + 2*t2 + 2*t8 + 2, t2*t6 - t2*t8 - 2*t4*t8 + t6 - t8, t3*t6 - t3*t8 - 2*t4*t8 - 2*t4, -2*t2*t5 + 2*t3*t5 - 4*t1 + 4*t2 - 2*t3 - 2*t5 + 4, t2*t5 + 2*t4*t5 - 2*t1*t6 + t2*t6 - t2 - 2*t4 + t5 + t6 - 1, t3*t5 + 2*t4*t5 - 2*t1*t6 + t3*t6 - 2*t1 - t3 - 4*t4, -2*t2*t6 + 2*t3*t6 - 2*t2 - 4*t4 - 2*t6 - 2, t1*t7 - t1*t8 + t5 + t8 - 1, t2*t7 - t2*t8 + t7 + t8, t3*t7 - t3*t8 + 2*t8 + 2, t2*t5 - 2*t1*t7 + t2*t7 - t2 - t5 + t7 + 1, t3*t5 - 2*t1*t7 + t3*t7 - 2*t1 - t3 - 2*t5 + 4, -2*t2*t7 + 2*t3*t7 - 2*t2 - 2*t7 + 2, t2*t5 - 2*t1*t8 + t2*t8 - t2 + t5 + t8 - 1, t3*t5 - 2*t1*t8 + t3*t8 - 2*t1 - t3, -2*t2*t8 + 2*t3*t8 - 2*t2 - 2*t8 - 2, 2*t2*t5 - 2*t3*t5 + 4*t1 - 4*t2 + 2*t3 + 2*t5 - 4, t2*t7 - t2*t8 + t7 + t8, t3*t7 - t3*t8 + 2*t8 + 2, t4*t7 - t4*t8 + t6 - t8, 2*t2*t7 - 2*t3*t7 + 2*t2 + 2*t7 - 2, t2*t6 - t2*t7 - 2*t4*t7 - t6 - t7, t3*t6 - t3*t7 - 2*t4*t7 - 2*t4 - 2*t6 - 2, 2*t2*t8 - 2*t3*t8 + 2*t2 + 2*t8 + 2, t2*t6 - t2*t8 - 2*t4*t8 + t6 - t8, t3*t6 - t3*t8 - 2*t4*t8 - 2*t4, 2*t2*t6 - 2*t3*t6 + 2*t2 + 4*t4 + 2*t6 + 2, t2*t7 - t2*t8 + t7 + t8, t3*t7 - t3*t8 + 2*t8 + 2, 2*t2*t7 - 2*t3*t7 + 2*t2 + 2*t7 - 2, 2*t2*t8 - 2*t3*t8 + 2*t2 + 2*t8 + 2 ]
 
Q = QuotientRing(R,ideal(E))

L = GF(101)
S = PolynomialRing(L,9,'t')
# Q.defining_ideal()
W = QuotientRing(S,ideal(E))

