import math

class PointVector(object):
    """Class that represents a nd point or vector."""

    #-----------basic services--------------------------

    def __init__(self, param):
        """Construction of the vector.

        Initialization by either a list of components,
        a tuple of components or another vector.
        :param param: components given as a list, tuple or PointVector.
        :raises TypeError for any other types. 
        """
        if isinstance(param, list):
            self.coords = tuple(param)
        elif isinstance(param,tuple):
            self.coords = param
        elif isinstance(param, PointVector):
            self.coords = param.coords
        else:
            raise TypeError
        
    def __repr__(self):
        return "("+",".join([str(x) for x in self.coords])+")"

    def __str__(self):
        return self.__repr__()

    def toTuple(self):
        """Returns the vector represented as a tuple of components"""
        return self.coords

    def dim(self):
        """Returns the dimension, ie. the number of components"""
        return len(self.coords)

    def norm(self, type):
        """Returns the norm of the vector.

        The returned norm can be either the L1, L2 or Linf one.
        :param type: norm type (one of the following strings: L1, L2, Linf).
        :raises ValueError for any unknown norm type.
        """
        if type == "L1":
            return sum([ abs(x) for x in self.coords ])
        elif type == "L2":
            return math.sqrt(sum([ x**2 for x in self.coords ]))
        elif type == "Linf":
            return max([ abs(x) for x in self.coords ])
        else:
            raise ValueError

    #------------lexicographic total order-----------------------

    def __lt__(self, other):
        assert isinstance(other, PointVector)
        assert self.dim() == other.dim()
        return self.coords < other.coords
        
    def __eq__(self, other):
        return not self < other and not other < self
    def __ne__(self, other):
        return not self == other
    def __gt__(self, other):
        return other < self
    def __ge__(self, other):
        return not self < other
    def __le__(self, other):
        return not other < self

    def __hash__(self):
        return hash( self.coords )
        
    #-------------basic operations--------------------------------

    def __add__(self, other):
        assert isinstance(other, PointVector) and self.dim() == other.dim()
        return PointVector([ x+y for (x,y) in zip(self.coords, other.coords) ])

    def __neg__(self):
        return PointVector([ -x for x in self.coords ])
        
    def __sub__(self, other):
        assert isinstance(other, PointVector)
        return self + (-other)

    def __mul__(self,other):
        assert not isinstance(other, PointVector)
        return PointVector([ x*other for x in self.coords])

    def __div__(self,other):
        assert not isinstance(other, PointVector)
        return PointVector([ x/other for x in self.coords])

    def dot(self, other):
        assert isinstance(other, PointVector) and self.dim() == other.dim()
        return sum([ x*y for (x,y) in zip(self.coords, other.coords) ])

    def square(self):
        return self.dot(self)

    def crossProduct(self, other):
        assert isinstance(other, PointVector)
        assert self.dim() == 3 and other.dim() == 3
        return PointVector([ self.coords[1]*other.coords[2] - self.coords[2]*other.coords[1],
                             self.coords[2]*other.coords[0] - self.coords[0]*other.coords[2],
                             self.coords[0]*other.coords[1] - self.coords[1]*other.coords[0] ])
        
def diagonal(aValue, d):
    """ Function that returns an instance of PointVector.

    The d components are all equal to a given value. 
    :param aValue: value of the d components of the vector.
    :param d: dimension, ie. number of components. 
    """
    return PointVector([aValue] * d)
