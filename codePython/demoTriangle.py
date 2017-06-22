#!/usr/bin/python
import functools
from NormalComputer.PointVector import PointVector
from NormalComputer.DigitalPlane import DigitalPlane, isReduced, additiveReduction
from NormalComputer.TriangleComputer import TriangleComputer, squaredRadiusOfSphere

from sys import argv

#-------------------------------------------------------------
def usage(script):
    """Describes the usage of the script."""
    print(script + " a b c (mode)")
    print("where:")
    print("(a,b,c) are the components of the normal")
    print("mode {H, R(default)} is the algorithm mode")


#-------------------------------------------------------------    
if len(argv) < 4:
    usage(argv[0])
    exit(0)

#components
param = [int(x) for x in argv[1:4]]
#non negative
if not functools.reduce( lambda x,y: x and y, [x >= 0 for x in param] ):
    print("Components must be non negative")
    exit(0)
#sum of zeros
zeros = 0
for component in param:
    if component == 0:
        zeros += 1
if zeros >= 2:
    print("Only one component over three is allowed to be equal to zero")
    exit(0)
#normal
n = PointVector(param)
#predicate
plane = DigitalPlane(n)

#mode
mode = "R"
if len(argv) > 4:
    mode = argv[4]
    
#basis
o = PointVector([0]*3)
e0 = PointVector([1,0,0])
e1 = PointVector([0,1,0])
e2 = PointVector([0,0,1])
q = PointVector([1]*3)
s = PointVector(q)

#starting triangle
startingTriangle = [o + e0+e1, o + e1+e2, o + e2+e0]

#normal computer
nc = TriangleComputer(startingTriangle, q, s, plane )
nc.setMode(mode)
print("#", nc.mode)
print("# step [triangle vertices] (normal) [vertices remainder] neighborhood isReduced?")
print("#", 0, nc.v, nc.getNormal(), [ plane.remainder(x) for x in nc.v ], nc.printNeighborhood(), nc.isReduced())
#main loop
c = 1
counterNonReduced = 0
while nc.advance():
    basis = nc.getBasis()
    if not isReduced(basis[0], basis[1]): 
        counterNonReduced += 1
    print("#", c, nc.v, nc.getNormal(), [ plane.remainder(x) for x in nc.v ], nc.printNeighborhood(), nc.isReduced())
    c += 1

if mode == "H":
    print("Reduction of the result of the H-algorithm")
    nc.reduction()
    
b = nc.getBasis()
print(b, "is reduced?", isReduced(b[0], b[1]))
vec1, vec2, nbRed = additiveReduction(b[0], b[1])
print("Reduced basis: ", (vec1, vec2))
print("Number of reductions: ", nbRed)
print("Number of times where the temporary basis is not reduced:", counterNonReduced)
print("Output:", nc.getNormal())
