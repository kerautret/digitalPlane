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
    
    def __init__(self, aTransformationMatrix, aTerminalSet):
        """ initialization with a transformation matrix
        acting on a vector of increasing components"""
        self.M = aTransformationMatrix
        self.terminalSet = aTerminalSet

    def __call__(self, aVector):
        """ call operator that takes a vector and returns another one,
        with smaller components """
        return self.matrix(aVector)*aVector

    def matrix(self, aVector):
        P = permutation(aVector)
        return P.inv()*self.M*P

    def isInTerminalSet(self, aVector):
        return (aVector in self.terminalSet)
        

def oneSubstitution123FromMatrix(aMatrix):
    indexToLetters = { 0: '1', 1: '2', 2: '3' }
    d = { '1' : "", '2' : "", '3' : "" }
    for indexRow in range(3):
        for indexCol in range(3):
            d[indexToLetters[indexCol]] += indexToLetters[indexRow] * aMatrix[indexRow,indexCol]
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
                point = self.sigma.matrix().inv() * ( aPointedFace.p + abelianizationMap123(s) )
                res.append( PointedFace( point, j ) )
        return res

                
# #----------------------------------------------
    
# sigma = Substitution123( { '1': "12", '2': "2", '3': "3"} )
# w = sigma("123")
# print(w)
# print(abelianizationMap123(w))
# print(abelianizationMap123('1'))
# print(abelianizationMap123('2'))
# print(abelianizationMap123('3'))
# print(abelianizationMap123(''))
# print()
# print(sigma.matrix())
# print()
# print(permutation( Matrix( [[3], [1], [2]] ) ))
# print()

# print("poincare")
# terminalSet = [ Matrix( [[1], [0], [0]] ),
#                 Matrix( [[0], [1], [0]] ),
#                 Matrix( [[0], [0], [1]] ) ]
# poincare = ContinuedFraction3d( Matrix( [[1,0,0],[-1,1,0], [0,-1,1]] ))
# u = Matrix( [[2], [2], [3]] )
# print(u)
# while u not in terminalSet:
#     print(poincare.matrix(u).inv())
#     print(oneSubstitution123FromMatrix( poincare.matrix(u).inv() ))
#     u = poincare(u)
#     print(u)


# print()
# print( splitWords("132521","2") )
# print( splitWords("123","3") )
# print( splitWords("123","1") )
# print( splitWords("123","5") )
# print( splitWords("","3") )
# print()

# origin = Matrix( [[0], [0], [0]] )
# sigma1 = Substitution123( { '1': "123", '2': "23", '3': "3"} )
# e1 = GeneralizedSubstitution(sigma1)
# print (e1( PointedFace(origin, "3") ))
# print (e1( PointedFace(origin, "2") ))
# print (e1( PointedFace(origin, "1") ))