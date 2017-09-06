from sympy.matrices import *

class Word123(object):
    """ class that implements a word over the alphabet {'1','2','3'}.
    NB. '4', '5', '6' are conjugate letters of '1','2','3' respectively """

    def __init__(self,aString):
        """ init method
        :param: aString, any string
        """
        self.s = aString

    def __add__(self,other):
        return Word123(self.s + other.s)

    def __iadd__(self,other):
        self.s += other.s
        return self
        
    def abelianizationMap(self):
        """ :return: a column vector mapping (index of) letters
        to number of occurences """        
        d = { '1': 0, '2': 0, '3':0 }
        for c in self.s:
            if c == '1' or c == '2' or c == '3':
                d[c] += 1
            if c == '4':
                d['1'] -= 1
            if c == '5':
                d['2'] -= 1
            if c == '6':
                d['3'] -= 1
        return Matrix( [ d['1'], d['2'], d['3'] ] )        

    def split(self, aLetter):
        """ returns the list of pairs (p, s) such that the current word
        is equal to the concatenation of p . aLetter . s
        :param: self, current word 
        :param: aLetter, any letter
        :return: list of word pairs """
        res = []
        p = ""          #empty string
        s = self.s     #whole word string
        for l in self.s:
            s = s[1:]  #update suffix
            if l == aLetter:
                res.append( (Word123(p), Word123(s)) )
            p += l     #update prefix
        return res

        
    def group(self):
        """ Returns the description of the current word
        as a list of pairs letter-power (negative power for conjugate letters) """
        res = list()
        
        if self.s:
            l = self.s[0] #letter
            c = 1         #counter
            k = 1         #index
            while k < len(self.s):
                if self.s[k] == l:
                    c += 1
                else:
                    res.append( (l,c) )
                    l = self.s[k]
                    c = 1
                k += 1
            res.append( (l,c) )

        return res

    def simplify(self):
        res = ""
        
        if self.s: #if non empty word
            l = self.s[0] #letter
            k = 1         #index
            toAdd = True  #flag
            while k < len(self.s):
                if self.s[k] == Word123.conjugateLetter(l):
                    #a) if two consecutive letters cancel
                    l = self.s[k]
                    k += 1
                    if k >= len(self.s):
                        break
                else:
                    #b) otherwise
                    res += l
                #in both case we advance
                l = self.s[k]
                k += 1
            else: #store the last letter in case b)
                res += l
        
        self.s = res
        
    def conjugateLetter(aLetter):
        if aLetter == '1':
            return '4'
        elif aLetter == '2':
            return '5'
        elif aLetter == '3':
            return '6'
        elif aLetter == '4':
            return '1'
        elif aLetter == '5':
            return '2'
        elif aLetter == '6':
            return '3'

    def conjugate(self):
        """ function that returns the conjugate word """
        res = ""
        for l in reversed(self.s):
            res += Word123.conjugateLetter(l)
        return Word123(res)
        
    def __repr__(self):
        return self.s

    
class Endomorphism123(object):
    """ class that implements an endomorphism over the alphabet {1,2,3}"""

    def __init__(self, aMap):
        """ init method
        :param: aMap, map between letters and words
        """
        self.map = aMap
        self.map['4'] = aMap['1'].conjugate()
        self.map['5'] = aMap['2'].conjugate()
        self.map['6'] = aMap['3'].conjugate()
        
    def __call__(self,aWord):
        """
        call operator
        :param: aWord, any word
        :return: another word after the replacement of each letter by their associated words
        """
        s = ""
        for l in aWord.s:
            s += str(self.map[l])
        res =  Word123(s)
        res.simplify()
        return res
        
    def __repr__(self):
        return "" + self.map.__repr__()

    def matrix(self):
        return self.map['1'].abelianizationMap()\
                            .row_join( self.map['2'].abelianizationMap() )\
                            .row_join( self.map['3'].abelianizationMap() )
    
def oneEndomorphism123FromMatrix(aMatrix):
    indexToLetters = { 0: '1', 1: '2', 2: '3' }
    d = { '1' : Word123(""), '2' : Word123(""), '3' : Word123("") }
    
    for indexCol in range(3):
        for indexRow in range(3):
            d[indexToLetters[indexCol]] += Word123(indexToLetters[indexRow] * aMatrix[indexRow,indexCol])
    # for indexCol in range(3):
    #     k = indexCol
    #     d[indexToLetters[indexCol]] += Word123(indexToLetters[k] * aMatrix[k,indexCol])
    #     for _ in range(2):
    #         k = (k+1)%3
    #         d[indexToLetters[indexCol]] += Word123(indexToLetters[k] * aMatrix[k,indexCol])
    return Endomorphism123( d )

class PointedFace(object):
    """ class that implements a pointed face, ie. a pair point - letter,
    where the letter identifies the face type """
    
    def __init__(self, aPoint, aFaceLetter):
        """ init method
        :param: aPoint, any point
        :param: aFaceLetter, any letter of the alphabet
        """
        self.p = aPoint
        self.l = aFaceLetter

    def __hash__(self):
        return hash( tuple([ tuple([c for c in self.p]), self.l ]) )

    def __repr__(self):
        return "("+ str([c for c in self.p]) + "," + self.l + "*)"
    
class DualMap(object):
    """ class that implements a dual associated with a unimodular morphism """

    def __init__(self, aSigma):
        """ param: aSigma, an morphism over the alphabet {'1','2','3'} """
        self.sigma = aSigma
        self.m = self.sigma.matrix().inv()
        self.e = { "1": Matrix( [ [1],[0],[0] ] ),
                   "2": Matrix( [ [0],[1],[0] ] ),
                   "3": Matrix( [ [0],[0],[1] ] ) }
        
    def __call__(self, aPointedFace):
        """ Function that sustitutes a pointed face by several others
        param: aPointedFace, ie. a pair (point, letter) 
        return: a dict of pointed faces
        """
        aPoint = aPointedFace.p
        aFaceLetter = aPointedFace.l
        
        d = dict()

        #positive part
        for j in ['1', '2', '3']:
            w = self.sigma.map[j]
            for p,s in w.split(aFaceLetter):
                point = self.m * ( aPoint - p.abelianizationMap() )
                newPointedFace = PointedFace( point, j )
                if newPointedFace in d: 
                    d[ newPointedFace ] += 1
                else:
                    d[ newPointedFace ] = 1

        #negative part
        for j in ['1', '2', '3']:
            w = self.sigma.map[j]
            for p,s in w.split(Word123.conjugateLetter(aFaceLetter)): 
                point = self.m * ( aPoint - p.abelianizationMap() + e[aFaceLetter] )
                newPointedFace = PointedFace( point, j )
                if newPointedFace in d: 
                    d[ newPointedFace ] -= 1
                else:
                    d[ newPointedFace ] = -1

        return d

#-----------------------------------------------

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
        mat[(k+2)%3,k] = -alpha #k-1
        mat[(k+1)%3,k] = -beta  #k-2
        # print(k, alpha, beta)
        # print(mat)
        return mat
        
