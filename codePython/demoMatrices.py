import functools
import argparse

from NormalComputer.PointVector import PointVector
from NormalComputer.DigitalPlane import DigitalPlane
from sympy.matrices import *
from pattern import ContinuedFraction3dByTriangleComputer

#-------------------------------------------------------------

#parse command line    
parser = argparse.ArgumentParser(description="computation of the normal vector of a digital plane")
parser.add_argument("x",help="x-component of the normal",type=int)
parser.add_argument("y",help="y-component of the normal",type=int)
parser.add_argument("z",help="z-component of the normal",type=int)
parser.add_argument("-m", "--algorithmMode",help="mode of the algorithm (default R)",
                    choices=["H","R"],default="R")

#-------------------------------------------------------------    

args = parser.parse_args()
param = [args.x, args.y, args.z]

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
if args.algorithmMode == "H":
    mode = "H"
    
#-------------------------------------------------------------    

algo = ContinuedFraction3dByTriangleComputer(n, mode)

matrices = []
while algo.advance():
    matrices.append( algo.getLastMatrix() )

#print( [ m.transpose() for m in matrices ] )

#------------------------------------------------------------
basis = [ Matrix([1,0,0]), Matrix([0,1,0]), Matrix([0,0,1]) ]

u = Matrix([c for c in n])
print("heights of ", u)

#---
#Soient M0...Mn la seq de matrices qui transforme les vecteurs de base ek en les vecteurs montant mk^n
#1) les vecteurs colonnes de la matrice M0...Mn donnent les vecteurs montant mk^n
#leur produit scalaire avec le vecteur u (la vraie normale) donne les hauteurs (donc de a,b,c, à 1,1,1)
M = eye(3)
print ([ [c for c in M.col(indexCol)] for indexCol in range(3) ], [ M.col(indexCol).dot(u) for indexCol in range(3)] )
for m in matrices:
    M = M * m
    print ([ [c for c in M.col(indexCol)] for indexCol in range(3) ], [ M.col(indexCol).dot(u) for indexCol in range(3)] )
print(M)
print(M.transpose())    

#---
#2) on utilise l'identité <a,b> = <M a,(M^T)^-1 b>
#ce qui donne <M0..n ek, u> (la hauteur de la k-ieme colonne de M0..n) est égale à
# < ek, (M0..n)^T u> (la k-ieme composante de (M0..n)^T u)
M = eye(3)
print ([ c for c in M.transpose().dot(u) ])
for m in matrices:
    M = M * m
    print (M.transpose())
    print ([ c for c in M.transpose().dot(u) ]) 
    
#---
#3) la suite des normales est ( M0...Mn e0 ^ M0...Mn e1 + M0...Mn e1 ^ M0...Mn e2 + M0...Mn e2 ^ M0...Mn e0)
#TODO a simplifier
def getNormal(M):
    #return (M.col(1) - M.col(0)).cross(M.col(2) - M.col(0))
    #or equivalently
    return (M.col(0).cross(M.col(1)) + M.col(1).cross(M.col(2)) + M.col(2).cross(M.col(0)))

print("normal")
M=eye(3)
print ( [c for c in getNormal(M)] )
for m in matrices:
    M = M * m
    print ( [c for c in getNormal(M)] )
    