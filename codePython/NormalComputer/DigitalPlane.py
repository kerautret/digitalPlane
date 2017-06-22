from NormalComputer.PointVector import PointVector
from collections import deque

class DigitalPlane(object):
    """Class that represents a digital plane"""

    def __init__(self, parameters, bound = 0):
        """Construction of a digital plane.

        The components of the normal vector must be given as
        a list of integers or as a PointVector.
        :param parameters: components of the normal vector.
        :param bound: intercept (default, 0).
        :raises TypeError if the type of the normal vector is
        unknown. 
        """
        self.bound = bound
        if isinstance(parameters,list):
            self.normal = PointVector(parameters)
        elif isinstance(parameters,PointVector):
            self.normal = parameters
        else:
            raise TypeError

    def __call__(self, x):
        """Predicate operator.
        
        :param x: any instance of PointVector
        :return: 'True' if x is in the digital plane, 'False' otherwise.
        """
        r = self.remainder(x)
        return self.bound <= r and r < self.normal.norm("L1")

    def remainder(self, x):
        """Returns the remainder of a given point."""
        return self.normal.dot(x) - self.bound

    def __repr__(self):
        return str(self.normal) + ", " + str(self.bound)
        
def neighborhood(current):
    """Returns the 6 neighbors of the given point as a list.
    
    :param current: any point
    :return: list of the 6 neighbors of the given point
    """
    return [ current + PointVector([1,0,0]),
             current + PointVector([-1,0,0]),
             current + PointVector([0,1,0]),
             current + PointVector([0,-1,0]),
             current + PointVector([0,0,1]),
             current + PointVector([0,0,-1])
         ]

def generateDFS(startingPoint, predicate, nodes):
    """Generate a digital point set, basically a piece of digital plane,
    by an iterative depth-first search.
    
    :param startingPoint: seed of the search
    (the predicate must return 'True' for it).
    :param predicate: point predicate that returns 'True' is the input point
    lies inside some bounding box, 'False' otherwise. 
    :param nodes: set of visited points (input/output).
    """
    assert(predicate(startingPoint))

    stack = []
    stack.append(startingPoint)
    while len(stack) > 0:
        x = stack.pop()
        if x not in nodes and predicate(x):
            nodes.add(x)
            for n in neighborhood(x):
                stack.append(n)

def generateBFS(startingPoint, predicate, nodes, edges):
    """Generate the graph elements (nodes and edges) of a digital point set, 
    by an iterative breadth-first search.
    
    :param startingPoint: seed of the search.
    :param predicate: point predicate that returns 'True' is the input point
    lies inside some bounding box, 'False' otherwise. 
    :param nodes: set of visited points (input/output).
    :param edges: set of visited edges, where both ends belong to the set of nodes
    (input/output).
    """
    assert(predicate(startingPoint))
    nodes.add(startingPoint)

    queue = deque()
    queue.append(startingPoint)

    while len(queue) > 0:
        node = queue.pop()
        for v in neighborhood(node):
            if predicate(v):
                if v not in nodes and predicate(v):
                    queue.append(v)
                    nodes.add(v)
                if tuple([node,v]) not in edges and tuple([v,node]) not in edges:
                    edges.add(tuple([node,v]))
                

def isReduced(u, v):
    """Returns 'True' if u and v form a reduced basis."""
    w = u+v
    x = u-v
    return u.square() <= w.square() and u.square()  <= x.square() \
        and v.square() <= w.square() and v.square()  <= x.square()

def additiveReduction(u, v):
    """Performs an additive reduction from the basis (u,v).

    :return: a tuple of 3 elements: u', v' (the vectors of the reduced basis)
    and the number of reduction steps."""

    counter = 0
    flagIsReduced = isReduced(u,v) 
    while not flagIsReduced:
        
        if u.square() < v.square():
            u, v = v, u
        assert(v.square() <= u.square())
        #v is smaller than u
        
        w = u+v
        x = u-v
        minVec = x
        if w.square() < x.square():
            minVec = w

        if minVec.square() < u.square():
            u = minVec
            counter += 1
        else:
            flagIsReduced = True

    assert(isReduced(u,v))
    return u, v, counter