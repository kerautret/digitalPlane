from NormalComputer.PointVector import PointVector
from sympy.matrices import *

m = Matrix( [1,2,5] )
#m = m.transpose()
n = Matrix( m )
for x in m:
    print(x)
print(m)
print("point")
p = PointVector( m )#[1,2,5] )
for x in p:
    print(x)
print(p)
print(p.tolist())
print(p.rows)
for paire in zip(p,m):
    print(paire)
print(m.dot(p))
print(p.dot(m))
print(p.norm("L1"))
print(p.norm("L2"))
print(p.norm("Linf"))
print(p.square())

v1 = Matrix([1,2,3])
v2 = Matrix([1,2,5])
print(v1,v2, v1+v2)

o = PointVector([0]*3)
e0 = PointVector([1,0,0])
e1 = PointVector([0,1,0])
e2 = PointVector([0,0,1])
q = PointVector([1]*3)
s = PointVector(q)

print(o,e0)
print(o+e0)

#starting triangle
print(o + e0+e1, o + e1+e2, o + e2+e0)
