import argparse

from fractions import gcd

from NormalComputer.PointVector import PointVector
from NormalComputer.DigitalPlane import DigitalPlane, isReduced, additiveReduction
from NormalComputer.TriangleComputer import TriangleComputer

#-------------------------------------------------------------

#parse command line    
parser = argparse.ArgumentParser(description="R-algorithm testing for a normal vector whose greatest component is below a threshold")
parser.add_argument("-t", "--threshold",help="maximal value for the greatest component of the normal vector",type=int, default=10)

#-------------------------------------------------------------    

args = parser.parse_args()
maxTests = args.threshold
print(maxTests)

o = PointVector([0]*3)
e0 = PointVector([1,0,0])
e1 = PointVector([0,1,0])
e2 = PointVector([0,0,1])
s = PointVector([1]*3)
q = s

counter = 0
for a in range(1,maxTests):
    for b in range(a,maxTests):
        for c in range(b,maxTests):
            d = gcd(a,b)
            if gcd(d,c) == 1:
                counter += 1
                n = PointVector([a,b,c])
                print("#", n)
                plane = DigitalPlane(n)
                nc = TriangleComputer([o+e0+e1, o+e1+e2, o+e2+e0], q, s, plane )
                while nc.advance():
                    pass
                nest = nc.getNormal()
                basis = nc.getBasis()
                assert(nest == n)
                assert(isReduced(basis[0], basis[1]))


            

