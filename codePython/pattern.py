from sympy.matrices import *


def abelianizationMap123(aWord):
    """ function that returns the number of occurences of each letter of a given word"""
    occurences = { '1' : 0, '2' : 0, '3' : 0 }
    for l in aWord:
        if l in occurences:
            occurences[l] += 1
    return Matrix( [ occurences['1'],
                     occurences['2'],
                     occurences['3'] ] )


class Substitution123(object):
    """ class that implements a substitution over the alphabet {1,2,3}"""

    def __init__(self,aMap):
        """ param: aMap, a dictionnary that associates letters 1, 2, 3 to words """
        self.map = aMap

    def __call__(self,aWord):
        """
        param: aWord, any word
        return: another word after substitution of its letters by the associated words
        """
        res = ""
        for l in aWord:
            res += self.map[l]
        return res

    def __repr__(self):
        return "" + self.map.__repr__()

    def matrix(self):
        return abelianizationMap123( self.map['1'] )\
            .row_join( abelianizationMap123( self.map['2'] ) )\
            .row_join( abelianizationMap123( self.map['3'] ) )
        
def permutation(v):
    """ computes a permutation matrix to sort the input vector components.

    param: v, any vector of 3 components
    return: permutation matrix such that the components of Mv are sorted
    """
    components = [ v[index,0] for index in range(3) ]
    ranks = [pair[0] for pair in sorted(enumerate(components), key=lambda pair: pair[1]) ]
    res = Matrix(3,3,[0]*9)
    for (index, rank) in enumerate(ranks):
        res[index,rank] = 1
    return res
    

class ContinuedFraction3d(object):
    """class that implements one iteration of a 3d continued fraction algorithm
    as a functor"""
    
    def __init__(self, aTransformationMatrix, aTerminalSet, aVector):
        """ initialization with a transformation matrix
        acting on a vector of increasing components"""
        self.M = aTransformationMatrix
        self.terminalSet = aTerminalSet
        self.lastMatrix = None
        self.V = aVector

    def advance(self):
        """ updates the current vector by another one with smaller components
        if not in the terminal set """
        if self.V in self.terminalSet:
            return False
        else: 
            self.lastMatrix = self.getNextMatrix()
            self.V = self.lastMatrix*self.V
            return True 

    def getNextMatrix(self):
        P = permutation(self.V)
        return P.inv()*self.M*P

    def getLastMatrix(self):
        return self.lastMatrix


def oneSubstitution123FromMatrix(aMatrix):
    print(aMatrix)
    indexToLetters = { 0: '1', 1: '2', 2: '3' }
    d = { '1' : "", '2' : "", '3' : "" }
    for indexCol in range(3):
        k = indexCol
        d[indexToLetters[indexCol]] += (indexToLetters[k] * aMatrix[k,indexCol])
        for _ in range(2):
            k = (k+1)%3
            d[indexToLetters[indexCol]] += (indexToLetters[k] * aMatrix[k,indexCol])
        # for indexCol in range(3):
        #     d[indexToLetters[indexCol]] = (indexToLetters[indexRow] * aMatrix[indexRow,indexCol]) + d[indexToLetters[indexCol]]
    return Substitution123( d )

def splitWords(aWord, aLetter):
    """ returns the list of pairs p, s such that aWord = p aLetter s """
    res = []
    p = '' #empty word
    s = aWord #whole word
    for l in aWord:
        s = s[1:]
        if l == aLetter:
            res.append( [p,s] )
        p += l
    return res
            
    
class PointedFace(object):

    def __init__( self, aPoint, aType ):
        """ initialization of a pointed face by a point and a face type """
        self.p = aPoint
        self.i = aType

    def __repr__(self):
        return "(" + self.p.__repr__() + "," + self.i + "*)"

class GeneralizedSubstitution(object):
    """ class that implements a generalized substitution """

    def __init__(self, aSigma):
        """ param: aSigma, a substitution over the alphabet {1,2,3} """
        self.sigma = aSigma

    def __call__(self, aPointedFace):
        """ Function that sustitutes a pointed face by several others
        param: aPointedFace, any pointed face
        return: a list of pointed faces
        """
        res = []
        for j in ['1', '2', '3']:
            for p,s in splitWords(self.sigma(j), aPointedFace.i):
                point = self.sigma.matrix().inv() * ( aPointedFace.p - abelianizationMap123(p) )
                res.append( PointedFace( point, j ) )
        return res



#-----------------------------------------------

from NormalComputer.TriangleComputer import TriangleComputer
from NormalComputer.DigitalPlane import DigitalPlane
from NormalComputer.PointVector import PointVector

class ContinuedFraction3dByTriangleComputer(object):
    """class that implements a 3d continued fraction algorithm
    based on our triangle computer """
    
    def __init__(self, aVector):
        """ initialization with a transformation matrix
        acting on a vector of increasing components"""
        o = PointVector([0]*3)
        e0 = PointVector([1,0,0])
        e1 = PointVector([0,1,0])
        e2 = PointVector([0,0,1])
        q = PointVector([1]*3)
        s = PointVector(q)
        self.V = PointVector( aVector )
        self.plane = DigitalPlane( self.V )
        self.nc = TriangleComputer([o + e1+e2, o + e2+e0, o + e0+e1], q, s, self.plane)
        
    def advance(self):
        """ do one iteration, the vector has smaller compondents """
        res =  self.nc.advance()
        self.V = PointVector( [self.plane.remainder(self.nc.q - x) for x in self.nc.v]) 
        return res
        
    def getLastMatrix(self):
        """ returns matrix of the previous iteration """
        k, alpha, beta = self.nc.operations[-1] #last operation
        mat = eye(3)
        mat[(k+2)%3,k] = alpha
        mat[(k+1)%3,k] = beta
        return mat.inv()
        
