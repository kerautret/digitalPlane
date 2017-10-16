from sys import stderr
import numpy
from scipy.spatial import ConvexHull
from sympy.matrices import *
import math

from NormalComputer.PointVector import PointVector
from NormalComputer.DigitalPlane import isReduced

#------------------------------------------------------------------

class TriangleComputer(object):
    """Class for computing the normal vector onto a digital plane.

    It implements both the H- and R- algorithm.
    """
    
    def __init__(self, triangle, bezout, shift, predicate):
        """Construction of the computer.
        
        :param triangle: starting triangle
        :param bezout: Bezout point lying above the starting triangle.
        :param shift: shift vector (+-1,+-1,+-1) that provides the orthant
        of the normal vector.
        :param predicate: point predicate returning 'True' if the given
        point is in the digital plane, 'False' otherwise.
        """
        self.v = triangle
        self.q = bezout
        self.s = shift
        self.predicate = predicate
        self.mode = "R" #can be "H"
        self.operations = []
        assert(self.isValid())

    def setMode(self, mode):
        """Sets the algorithm to use (either H-, R-, CH-)
        :param mode: mode (one of the following strings: 'H', 'R', 'CH')
        :raises: ValueError for any unknown mode
        """
        if not (mode == "R" or mode == "H" or mode == "CH"):
            raise ValueError
        self.mode = mode

    def getHexagon(self):
        """ Returns the hexagon, lying above the triangle and around the Bezout point
        """
        m = [ self.q - x for x in self.v ]
        a = [ self.v[k] + m[k-1] for k in range(len(m)) ]
        b = [ self.v[k] + m[k-2] for k in range(len(m)) ]
        #list interlacing
        return [ item for sublist in zip(a,b[1:]+b[:1]) for item in sublist ]

    def getRays(self):
        """ Returns the rays as a list of pairs:
        - point (starting point of the ray)
        - vector (direction vector of the ray).
        """
        m = [ self.q - x for x in self.v ]
        a = [ (self.v[k] + m[k-1], m[k-2]) for k in range(len(m)) ]
        b = [ (self.v[k] + m[k-2], m[k-1]) for k in range(len(m)) ]
        #list interlacing
        return [ item for sublist in zip(a,b[1:]+b[:1]) for item in sublist ]

    def neighborhood(self):
        """ Returns the neighborhood as a list of triples:
        - index (of the vertex to update: 0, 1 or 2)
        - point (starting point of the ray)
        - vector (direction vector of the ray).
        """
        m = [ self.q - x for x in self.v ]
        a = [ (k, self.v[k] + m[k-1], m[k-2]) for k in range(len(m)) ]
        b = [ (k, self.v[k] + m[k-2], m[k-1]) for k in range(len(m)) ]
        #list interlacing
        return [ item for sublist in zip(a,b[1:]+b[:1]) for item in sublist ]

    def printNeighborhood(self):
        """Returns a string representation of the current neighborhood.

        Neighborhood points lying in the digital plane are depicted with crosses,
        whereas other points are depicted with circles.
        """
        res = ""
        for x in self.getHexagon():
            if self.predicate(x):
                res += " x "
            else:
                res += " o "
        return res
        
    def isValid(self):
        """ returns 'True' if data are consistent and valid, 'False' otherwise."""

        res = len(self.v) == 3
        #1) triangle inside the plane
        for x in self.v:
            res = res and self.predicate(x)
        assert(res)
        #2) remainder of self.q
        nest = self.getNormal()
        res = res and ( (nest.dot(self.q) - nest.dot(self.v[0])) == 1)
        assert(res)
        #3) projection
        res = res and self.isProjectedInside(self.v)
        assert(res)
        #4) mode
        res = res and (self.mode == "R" or self.mode == "H" or self.mode == "CH")
        assert(res)

        return res

    def isProjectedInside(self, triangle):
        """ returns 'True' if the Bezout point projects
        (in the shift vector direction)
        inside the given triangle, 'False' otherwise."""
        assert(len(triangle) >= 3)
        res = True
        m = [ self.q - x for x in triangle ]
        for k in range(len(m)):
            nk = m[k-1].cross(m[k])
            res = res and ( nk.dot(-self.s) <= 0 )
        return res

    def circularEquality(l1, l2):
        """ returns 'True' if the two lists l1 and l2
        are circularly identical, 'False' otherwise."""
        if len(l1) != len(l2):
            return False
        else:
            res = False
            for k in range(len(l1)):
                if l2 == ( l1[-k:]+l1[:-k] ):
                    res = True
                    break
            return res
        
    def equalsTo(self, triangle):
        """ returns 'True' if the given triangle is
        circularly identical to self.v, 'False' otherwise. """
        return TriangleComputer.circularEquality(triangle,self.v)

    def ray(self, aStartingPoint, aVector):
        """Returns a closest point along a ray given by
        a starting point and a direction vector.

        :warning: the starting point must belong to the digital plane """
        assert(self.predicate(aStartingPoint))

        index = 0
        previousX = aStartingPoint
        if self.mode == "R": 
            currentX = aStartingPoint + aVector
            while self.predicate(currentX) and \
                  ( squaredRadiusOfSphere(self.v + [currentX]) < squaredRadiusOfSphere(self.v + [previousX]) ):
                index += 1
                previousX = currentX
                currentX = previousX + aVector
            assert(self.predicate(previousX))

        return (previousX, index)

    def advanceByCHalgorithm(self):
        """Updates the triangle according to the CH-algorithm."""
        res = False
        
        pointSet = [ x for x in self.getHexagon() if self.predicate(x) ]
        if len(pointSet) > 0:
            #0) filter
            closest = min(pointSet, key=lambda x: squaredRadiusOfSphere(self.v + [x]))
            minValue = squaredRadiusOfSphere(self.v + [closest])
            allClosest = [ x for x in pointSet if squaredRadiusOfSphere(self.v + [x]) == minValue ]
            #1) input data
            pointSet = allClosest + self.v
            n = len(pointSet)
            d = 3
            inputData = numpy.ndarray(shape=(n,d), dtype=int)
            for i in range(n):
                for j in range(d):
                    inputData[i,j] = pointSet[i][j]
            print(inputData)
            #2) convex hull computation        
            hull = ConvexHull(inputData)
            #3) find the upper facet where self.q projects into
            for s in hull.simplices:
                facet = [ PointVector(inputData[index].tolist()) for index in s ]
                print(facet)
                rFacet = list(reversed(facet)) #to take into account both orientations
                if not TriangleComputer.circularEquality(facet, self.v) and \
                   not TriangleComputer.circularEquality(rFacet, self.v):
                    if self.isProjectedInside(facet):
                        self.v = facet
                        res = True
                        break
                    elif self.isProjectedInside(rFacet):
                        self.v = rFacet
                        res = True
                        break
            if not res:
                raise RuntimeError("no facet found")
        return res

    def advanceByHorRalgorithm(self):
        """Updates the triangle according to the H- or R-algorithm."""
        res = False

        innerPoints = []

        m = [ self.q - x for x in self.v ]
        neighbors = [ (k, self.v[k] + m[k-1], m[k-2]) for k in range(len(m)) ]
        for triple in neighbors:
            if self.predicate(triple[1]):
                x, beta = self.ray(triple[1], triple[2])
                innerPoints.append( (triple[0], x, 1, beta) )
        neighbors = [ (k, self.v[k] + m[k-2], m[k-1]) for k in range(len(m)) ]
        for triple in neighbors:
            if self.predicate(triple[1]):
                x, alpha = self.ray(triple[1], triple[2])
                innerPoints.append( (triple[0], x, alpha, 1) )

        if innerPoints:
            self.update(min(innerPoints, key=lambda x: squaredRadiusOfSphere(self.v + [x[1]])))
            res = True
            assert(self.isValid())

        return res
            
            
    def advance(self):
        """Updates the triangle by the neighbors belonging to the digital plane."""
        
        if self.mode == 'CH':
            return self.advanceByCHalgorithm()
        else:
            return self.advanceByHorRalgorithm()
 
        return res
        
    def update(self, aTuple):
        """Updates the underlying triangle wrt 'aTuple', which is a pair index/point. """
        k = aTuple[0]
        self.v[k] = aTuple[1]
        self.operations.append( (k,aTuple[2],aTuple[3]) )

    def getBasis(self):
        """Returns the two shortest triangle edges as the basis. """
        u = self.v[1]-self.v[0]; 
        v = self.v[2]-self.v[1]; 
        w = self.v[0]-self.v[2];
        assert(w == (-u-v)); 

        if u.square() < v.square(): 
            if u.square() < w.square(): 
                #u has minimal length
                if (-w).square() < v.square():
                    return (u, -w)
                else:
                    return (u, v)
            else:
                #w has minimal length
                if (-v).square() < u.square():
                    return (w, -v)
                else:
                    return (w, u)
        else:
            if v.square() < w.square(): 
                #v has minimal length
                if (-u).square() < w.square():
                    return (v, -u)
                else:
                    return (v, w)
            else:  
                #w has minimal length
                if (-v).square() < u.square():
                    return (w, -v)
                else:
                    return (w, u)

    def getNormal(self):
        """Returns the normal vector of the underlying triangle. """
        b = self.getBasis()
        return b[0].cross(b[1])

    def isReduced(self):
        """Returns 'True' if the triangle is reduced, 'False' otherwise. """
        d = [ self.v[k] - self.v[k-1] for k in range(len(self.v)) ]
        d = d[1:] + d[:1]
        return (d[0].square() <= (d[1] - d[2]).square()) and \
            (d[1].square() <= (d[2] - d[0]).square()) and \
            (d[2].square() <= (d[0] - d[1]).square())
        
    def isReduced2(self):
        """Returns 'True' if the basis extracted from the triangle is reduced,
        'False' otherwise. """
        b = self.getBasis()
        return isReduced(b[0], b[1])

    def reductionOneStep(self):
        """Reduces the triangle by one step"""
        d = [ self.v[k] - self.v[k-1] for k in range(len(self.v)) ]
        d = d[1:] + d[:1]
        for k in range(len(self.v)):
            diagonal = d[k-1] - d[k-2] 
            if diagonal.square() < d[k].square():
                #if the diagonal is (strictly) shorter, we update one vertex
                newVertex = self.v[k-1] + diagonal
                if self.isProjectedInside([ self.v[k], newVertex, self.v[k-1] ]):
                    self.v[k-2] = newVertex
                else:
                    assert(self.isProjectedInside([ newVertex, self.v[k-2], self.v[k-1] ]))
                    self.v[k] = newVertex
                break
            
    def reduction(self):
        """Reduces the triangle until it is reduced"""
        counter = 0
        while not self.isReduced():
            self.reductionOneStep()
            counter += 1
        return counter
       
#------------------------------------------------------------------
#criteria for minimality
#------------------------------------------------------------------

def squaredRadiusOfSphere(lst):
    """ Function that returns the squared radius of the sphere passing
    by the four points of the given list.

    :raises: ValueError if the 4 points are coplanar.
    """
    assert(len(lst) == 4)

    lstV = [ lst[0] - X for X in lst[1:] ]
    lstOppV = [ lst[3] - lst[2], lst[1] - lst[3], lst[2] - lst[1] ]

    M = Matrix( [ [c for c in v] for v in lstV ])
    volume = M.det()
    if volume == 0:
        raise ValueError("ERR: 4 copanar points")

    a2, b2, c2 = [ x.square() for x in lstV ]
    ap2, bp2, cp2 = [ x.square() for x in lstOppV ]
    M2 = Matrix( [ [0, a2, b2, c2],
                  [a2, 0, cp2, bp2],
                  [b2, cp2, 0, ap2],
                  [c2, bp2, ap2, 0] ] )
    return -M2.det() / (16*volume**2)

def radiusOfSphere(lst):
    """ Function that returns the radius of the sphere passing
    by the four points of the given list.

    :raises: ValueError if the 4 points are coplanar.
    """
    return math.sqrt(float(squaredRadiusOfSphere(lst)))
       
