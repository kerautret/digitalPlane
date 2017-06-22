#!/usr/bin/python

from fractions import gcd

from NormalComputer.PointVector import PointVector
from NormalComputer.DigitalPlane import DigitalPlane, isReduced, additiveReduction
from NormalComputer.TriangleComputer import TriangleComputer

from sys import argv

def usage(script):
    print("usage: ")
    print(script + " <max> ")

if len(argv) < 2:
    usage(argv[0])
    exit(0)

maxTests = int(argv[1])
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


            

