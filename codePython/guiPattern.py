import tkinter as tk
import argparse

from sys import argv
from sympy.matrices import *
from copy import deepcopy

from pattern import ContinuedFraction3d, ContinuedFraction3dByTriangleComputer,\
    Word123, Endomorphism123, oneEndomorphism123FromMatrix, PointedFace, DualMap

#------------------------------------------------
#----------------- helpers ----------------------
#------------------------------------------------

#hexagonalProjection = ( (-0.5,-0.86602540378), (1,0), (-0.5, 0.86602540378) )
#standardProjection = ( (-0.353,-0.353), (1,0), (0,1) )
#NB. we take the opposite for the y-coordinate

def hexagonalProjector(vector3):
    return Matrix( [ [ -0.5, 1, -0.5], [ 0.86602540378, 0, - 0.86602540378] ] ) * vector3; 

def standardProjector(vector3):
    return Matrix( [ [ -0.355, 1, 0], [ 0.353, 0, -1] ] ) * vector3; 

def lowerUnitCube():
    return [ [PointedFace( Matrix( [[-1], [0], [0]] ) , "1")],
             [PointedFace( Matrix( [[0], [-1], [0]] ) , "2")],
             [PointedFace( Matrix( [[0], [0], [-1]] ) , "3")] ]

def upperUnitCube():
    return [ [PointedFace( Matrix( [[0], [0], [0]] ) , "1")],
             [PointedFace( Matrix( [[0], [0], [0]] ) , "2")],
             [PointedFace( Matrix( [[0], [0], [0]] ) , "3")] ]

def firstArg(aX, aY):
    return aX

def secondArg(aX, aY):
    return aY

def report_event(event):    
    """Print a description of an event, based on its attributes.
    """
    event_name = {"2": "KeyPress", "4": "ButtonPress"}
    print("EventType=" + str(event.type), \
        event_name[str(event.type)],\
        "EventWidgetId=" + str(event.widget), \
        "EventKeySymbol=" + str(event.keysym))
    
#------------------------------------------------  
#--------------- main class ---------------------
#------------------------------------------------
    
class PatternApp(object):
    """ Simple class for a tkinter application that draws patterns of digital plane """

    def __init__(self, parent, size, step, projector, corner, colorMode, algo):
        """ Initilization of the application

        :param parent: object of type Tk
        :param normal: plane normal 
        :param size: discrete size of the drawing window
        :param step: grid step of the drawing window
        :param projector: functor that takes a 3d vector and
        returns its 2d projection.
        :param corner: function that returns a starting corner
        :param colorMode: function that selects the relevant type
        """
        #dimensions
        self.discreteSize = size
        self.gridStep = step
        self.dotSize = step / 10.0

        #origin
        mid = self.discreteSize/2
        self.origin = Matrix( [[ mid,mid ]] ).transpose()
        
        #canvas
        realSize = self.gridStep * self.discreteSize
        self.canvas = tk.Canvas(parent, width=realSize, height=realSize, borderwidth=0, highlightthickness=0)
        self.canvas.pack()

        #projector
        self.projector = projector
        #starting corner
        self.unitCube = corner()
        #color
        self.mode = colorMode
        
        #map type <-> vectors / colors
        self.mapTypeVectors = { "1": ( Matrix( [ [0],[1],[0] ] ), Matrix( [ [0],[0],[1] ] ) ),
                                "2": ( Matrix( [ [0],[0],[1] ] ), Matrix( [ [1],[0],[0] ] ) ),
                                "3": ( Matrix( [ [1],[0],[0] ] ), Matrix( [ [0],[1],[0] ] ) ) }
        self.mapTypeVector = { "1": Matrix( [ [1],[0],[0] ] ),
                               "2": Matrix( [ [0],[1],[0] ] ),
                               "3": Matrix( [ [0],[0],[1] ] ) }
        self.mapTypeColor = { "1": "white", "2": "gray60", "3": "gray30" }

        #continued fraction expansion
        self.morphisms = []
        self.M = eye(3,3) #matrix storing the basis vectors
        print(algo.V)
        while algo.advance():
            self.morphisms.append( oneEndomorphism123FromMatrix( algo.getLastMatrix().inv() ))
            self.M *= algo.getLastMatrix()
            print(algo.V)
        print(self.M)
        print()
        
        self.morphisms.reverse()
        for s in self.morphisms:
            print(s)
            print(s.matrix())
            print(s.matrix().inv())
            
        self.start()
            
        #key binding
        self.canvas.bind_all( "<Right>", self.forward )
        self.canvas.bind_all( "<Return>", self.reset )
    
    #------------------------------------------------------------------------------------------

    def transform(self,x):
        """ returns the image of point x after a translation and a scaling"""
        return (x+self.origin)*self.gridStep

    def flatten(self, aListOf2dVectors):
        res = []
        for v in aListOf2dVectors:
            res.append(v[0])
            res.append(v[1])
        return res
        
    def drawParallelogram(self, aF, aK):
        x = aF.p + self.mapTypeVector[aF.l]
        v1, v2 = self.mapTypeVectors[aF.l]
        polygon3d = [ x, x + v1, x + v1 + v2, x + v2 ]
        polygon3dbis = [ self.currentM * x for x in polygon3d ]
        polygon2d = [ self.projector(x) for x in polygon3dbis ]
        coords = self.flatten( [ self.transform(x) for x in polygon2d] )
        self.canvas.create_polygon(coords, 
                                   fill=self.mapTypeColor[self.mode(aF.l, aK)],
                                   outline="black",
                                   width=2,
                                   tags="piece")
                                
    def drawPiece(self):
        self.canvas.delete("piece")
        
        for k,tileSet in enumerate(self.tiles):
            for tile in tileSet:
                self.drawParallelogram( tile, str(k+1) )
                
    def start(self):
        self.index = 0
        self.currentM = self.M
        self.tiles = self.unitCube
        self.drawPiece()

    def reset(self,event):
        self.start()
        
    def forward(self,event):
        report_event(event)
        if self.index < len(self.morphisms):
            #update matrix for embedding
            self.currentM *= self.morphisms[self.index].matrix()
            print(self.currentM)
            #do substitution
            s = DualMap( self.morphisms[self.index] )
            res = []
            for tileSet in self.tiles:
                res2 = []
                for tile in tileSet:
                    res2 += s(tile)
                print(res2)
                res.append(res2)
            self.tiles = res
            self.index += 1
            #draw
            self.drawPiece()
        
#------------------------------------------------  
#--------------- main script --------------------
#------------------------------------------------

#parse command line    
parser = argparse.ArgumentParser(description="draws a piece of digital plane from its normal vector")
parser.add_argument("x",help="x-component of the normal",type=int)
parser.add_argument("y",help="y-component of the normal",type=int)
parser.add_argument("z",help="z-component of the normal",type=int)
parser.add_argument("-w", "--windowSize",help="discrete size of the viewing window (default 10)",
                    type=int,default=10)
parser.add_argument("-u", "--unitSize",help="size of the discrete unit (default 40)",\
                    type=int,default=40)
parser.add_argument("-p", "--projection",help="type of 3d to 2d projection",\
                    choices=["standard","hexagonal"],default="standard")
parser.add_argument("-s", "--startingCorner",help="type of starting corner",\
                    choices=["lower","upper"],default="lower")
parser.add_argument("-c", "--color",help="either a color by tile type or a color by tile group",\
                    choices=["type","group"],default="type")
parser.add_argument("-a", "--algo",help="continued fraction expansion algorithm",\
                    choices=["poincare", "brun", "selmer", "H", "R"],default="poincare")


args = parser.parse_args()
param = [args.x, args.y, args.z]
n = Matrix( [param] ).transpose()

projector = standardProjector
if args.projection == "hexagonal":
    projector = hexagonalProjector
corner = lowerUnitCube
if args.startingCorner == "upper":
    corner = upperUnitCube
mode = firstArg
if args.color == "group":
    mode = secondArg
  
#poincare
algo = ContinuedFraction3d( Matrix( [[1,0,0],[-1,1,0],[0,-1,1]] ),\
                            [ Matrix( [[1], [0], [0]] ), \
                              Matrix( [[0], [1], [0]] ), \
                              Matrix( [[0], [0], [1]] ) ], n )
if args.algo == "brun":
    algo = ContinuedFraction3d( Matrix( [[1,0,0],[0,1,0],[0,-1,1]] ),\
                                [ Matrix( [[1], [0], [0]] ),\
                                  Matrix( [[0], [1], [0]] ),\
                                  Matrix( [[0], [0], [1]] ) ], n ) 
elif args.algo == "selmer":
    algo = ContinuedFraction3d( Matrix( [[1,0,0],[0,1,0],[-1,0,1]] ),\
                                [ Matrix( [[1], [1], [1]] ) ], n ) 
elif args.algo == "H" or args.algo == "R":
    algo = ContinuedFraction3dByTriangleComputer(param, args.algo)
    
#application
root = tk.Tk()
c = PatternApp(root, args.windowSize, args.unitSize, projector, corner, mode, algo)     
root.mainloop()

    

