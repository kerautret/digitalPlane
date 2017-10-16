from sympy.matrices import *
from pattern import Word123, Endomorphism123, oneEndomorphism123FromMatrix, PointedFace, DualMap

#w = ( ('2',-1), ('3',1) )
# w = { '2':-1, '3':1 }
# w2 = wordFromString("1231231231")
# print (w, abelianizationMap123(w))
# print (w2, abelianizationMap123(w2))

w = Word123("11223344556611")
print (w, w.abelianizationMap())
print(w.split('1'))

w3 = Word123("541")
print(w3)
w3.simplify()
print(w3)
w4 = Word123("25")
print(w4)
w4.simplify()
print(w4)

rauzySubstitution = Endomorphism123( { '1': Word123("12"), '2': Word123("13"), '3': Word123("1") } )
print(rauzySubstitution.map)
w2 = rauzySubstitution(Word123("42"))
print(w2, w2.group())

m = rauzySubstitution.matrix()
print(m)
print(oneEndomorphism123FromMatrix(m))


dm = DualMap(rauzySubstitution)
tile = PointedFace( Matrix( [[0], [0], [0]] ) , "1") 
tileset = dm( tile )
print(tileset)
tile = PointedFace( Matrix( [[0], [0], [0]] ) , "2") 
tileset = dm( tile )
print(tileset)
tile = PointedFace( Matrix( [[0], [0], [0]] ) , "3") 
tileset = dm( tile )
print(tileset)