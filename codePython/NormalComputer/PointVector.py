import math
from sympy.matrices import *

class PointVector(Matrix):

    def square(self):
        return self.dot(self)

    def dim(self):
        """Returns the dimension, ie. the number of components"""
        return self.rows

    def norm(self, type):
        """Returns the norm of the vector.

        The returned norm can be either the L1, L2 or Linf one.
        :param type: norm type (one of the following strings: L1, L2, Linf).
        :raises ValueError for any unknown norm type.
        """
        if type == "L1":
            return sum([ abs(x) for x in self ])
        elif type == "L2":
            return math.sqrt(sum([ x**2 for x in self ]))
        elif type == "Linf":
            return max([ abs(x) for x in self ])
        else:
            raise ValueError

    def __repr__(self):
        return "("+ ",".join( [str(x) for x in self] ) +")"

    def __hash__(self):
        return hash(tuple([x for x in self]))