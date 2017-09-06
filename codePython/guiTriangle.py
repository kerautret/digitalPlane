import tkinter as tk
import argparse

from copy import deepcopy
from sympy.matrices import *

from NormalComputer.PointVector import PointVector
from NormalComputer.DigitalPlane import DigitalPlane, generateBFS, isReduced, additiveReduction
from NormalComputer.TriangleComputer import TriangleComputer, squaredRadiusOfSphere

#------------------------------------------------
#----------------- helpers ----------------------
#------------------------------------------------

#hexagonalProjection = ( (-0.5,-0.86602540378), (1,0), (-0.5, 0.86602540378) )
#standardProjection = ( (-0.353,-0.353), (1,0), (0,1) )
#NB. we take the opposite for the y-coordinate

def hexagonalProjector(vector3):
    return vector3.projection( [ [ -0.5, 1, -0.5], [ 0.86602540378, 0, - 0.86602540378] ] ); 

def standardProjector(vector3):
    return vector3.projection( [ [ -0.355, 1, 0], [ 0.353, 0, -1] ] ); 

def firstArg(aX, aY):
    return aX

def secondArg(aX, aY):
    return aY
    
def _create_circle(self, x, y, r, **kwargs):
    return self.create_oval(x-r, y-r, x+r, y+r, **kwargs)
tk.Canvas.create_circle = _create_circle

def report_event(event):    
    """Print a description of an event, based on its attributes.
    """
    event_name = {"2": "KeyPress", "4": "ButtonPress"}
    print("EventType=" + str(event.type), \
        event_name[str(event.type)],\
        "EventWidgetId=" + str(event.widget), \
        "EventKeySymbol=" + str(event.keysym))

#------------------------------------------------  
#--------------- helper classes -----------------
#------------------------------------------------
class Tile(object):

    colorByType = { PointVector([-1,0,0]): "grey60",
                    PointVector([0,-1,0]): "grey30",
                    PointVector([0,0,-1]): "white", }
    colorByGroup = ["white", "gray60", "gray30"]
    
    def __init__(self, origin, v1, v2, normal):
        self.o = origin
        self.v1 = v1
        self.v2 = v2
        self.n = normal
        
    def origin(self):
        return self.o
        
    def color(self):
        return self.c
        
    def shift(self, v):
        self.o += v

class ShiftedTile(Tile):

    def origin(self):
        return self.o + self.n
        
    def shift(self, v):
        self.o -= v

class PatternMode:
    pass
        
#------------------------------------------------  
#--------------- main class ---------------------
#------------------------------------------------
    
class TriangleComputerApp(object):
    """ Simple class for a tkinter application that draws the normal computer output onto a digital plane """

    #-------------------------------------------------------------------------
    def transform(self,x):
        """ returns the image of point x after a translation and a scaling"""
        return (x+self.origin)*self.gridStep

    def flatten(self, aListOfVectors):
        res = []
        for v in aListOfVectors:
            for c in v: 
                res.append(c)
        return res
    #-------------------------------------------------------------------------
        
    def __init__(self, parent, normal, size, step, projector, patternMode):
        """ Initilization of the application

        :param parent: object of type Tk
        :param normal: plane normal (provided as a PointVector)
        :param size: discrete size of the drawing window
        :param step: grid step of the drawing window
        :param projector: projection from a 3d PointVector
        to a 2d PointVector.
        :param patternMode: specify the color of the tiles or the starting
        set of tiles, useful when drawing pattern is enabled
        """
        #dimensions
        self.discreteSize = size
        self.gridStep = step
        self.dotSize = step / 10.0

        #origin
        mid = self.discreteSize/2
        self.origin = PointVector( [ mid ] *2 )
        
        #canvas
        realSize = self.gridStep * self.discreteSize
        self.canvas = tk.Canvas(parent, width=realSize, height=realSize, borderwidth=0, highlightthickness=0)
        self.canvas.pack()

        #projector
        self.projector = projector
        #mode
        self.mode = patternMode

        #digital plane
        o = PointVector([0] * 3)
        self.plane = DigitalPlane(normal)
        pointSet3 = set()
        edgeSet3 = set()
        generateBFS(o, lambda x: self.plane(x) and x.norm("Linf") <= size,  pointSet3, edgeSet3)
           
        #digital plane as dictionnary of 3d points
        for e3 in edgeSet3:
            e2 = [self.projector(x) for x in e3]
            coords = self.flatten( [ self.transform(x) for x in e2] )            
            self.canvas.create_line(coords,
                                    fill="gray", width=2)

        self.grid = {}
        for p3 in pointSet3:
            p2 = self.transform(self.projector(p3))
            pointKey = self.canvas.create_circle( p2[0], p2[1],
                                                  self.dotSize, fill="black")
            self.canvas.tag_bind( pointKey, "<Button-1>", self.selectStartingPoint )
            self.grid[pointKey] = p3

        #default options
        self.enableH = False
        self.enableR = False
        self.enableT = False
        self.enableP = False  

        #key binding
        self.canvas.bind_all( "<Right>", self.updateForwardAndDraw )
        self.canvas.bind_all( "<Left>", self.updateBackwardAndDraw )
        self.canvas.bind_all( "<Return>", self.reSetTriangle )
        self.canvas.bind_all( "<Key>", self.setOption )

        #data
        self.q = PointVector([1]*3)
        self.s = PointVector([1]*3)
        self.e0 = PointVector([1,0,0])
        self.e1 = PointVector([0,1,0])
        self.e2 = PointVector([0,0,1])

        self.drawBasis()
        self.start()

    #-------------------------------------------------------------------------
    #-------------------------------- actions --------------------------------
    #-------------------------------------------------------------------------

    def selectStartingPoint(self,event):
        report_event(event)
        pointKey = self.canvas.find_withtag(tk.CURRENT)[0] #current item id
        self.q = self.grid[pointKey] + PointVector([1]*3)
        self.start()
        
    def start(self):

        #normal computer to get the list of operations
        o = self.q - self.s
        nc = TriangleComputer([o+self.e0+self.e1, o+self.e1+self.e2, o+self.e2+self.e0],
                              self.q, self.s, self.plane )
        while nc.advance():
            pass

        self.operations = nc.operations

        #other data
        self.initPattern()
        
        #draw first triangle
        self.drawTriangle()
        
    def initPattern(self): 
        self.index = 0
        self.m = [self.e2, self.e0, self.e1]
        if self.mode.corner == "upper":
            self.tiles = [ [ Tile(self.q, -self.m[1], -self.m[2], -self.m[0]) ],
                           [ Tile(self.q, -self.m[2], -self.m[0], -self.m[1]) ],
                           [ Tile(self.q, -self.m[0], -self.m[1], -self.m[2]) ] ]
        elif self.mode.corner == "lower":
            self.tiles = [ [ ShiftedTile(self.q, -self.m[1], -self.m[2], -self.m[0]) ],
                           [ ShiftedTile(self.q, -self.m[2], -self.m[0], -self.m[1]) ],
                           [ ShiftedTile(self.q, -self.m[0], -self.m[1], -self.m[2]) ] ]
        else:
            raise ValueError
        
    def drawTriangle(self): 
        self.canvas.delete("triangle")
        
        if self.enableP:
            self.drawPiece()

        #draw reference points
        ref = self.transform(self.projector(self.q))
        pointKey = self.canvas.create_circle( ref[0], ref[1],
                                              self.dotSize, fill="blue",tags="triangle")

        v = [ self.q - mk for mk in self.m ]
        t = [ self.projector(x) for x in v ]
        coords = self.flatten( [ self.transform(x) for x in t] )
        self.triangleItemId = self.canvas.create_polygon(coords, fill="",outline="blue",width=3, tags="triangle")
                
        if self.enableH:
            self.drawHexagon()
        if self.enableR:
            self.drawRays()
        if self.enableT:
            self.drawBigTriangle()

    def drawHexagon(self):
        self.canvas.delete("hexagon")

        v = [ self.q - mk for mk in self.m ]
        a = [ v[k] + self.m[k-1] for k in range(len(self.m)) ]
        b = [ v[k] + self.m[k-2] for k in range(len(self.m)) ]
        hexagon3d = [ item for sublist in zip(a,b[1:]+b[:1]) for item in sublist ]
        for x in hexagon3d:
            p2 = self.transform(self.projector(x))
            if self.plane(x):
                self.canvas.create_circle( p2[0],p2[1],
                                           self.dotSize, fill="red", outline="red", tags="hexagon")
            else:
                self.canvas.create_circle( p2[0],p2[1],
                                           self.dotSize+2, outline="red", width=1, tags="hexagon")

        hexagon2d = [ self.projector(x) for x in hexagon3d ]
        coords = self.flatten( [ self.transform(x) for x in hexagon2d] )
        self.canvas.create_polygon(coords, fill="",outline="red",width=2, tags="hexagon")


    def drawRays(self):
        self.canvas.delete("rays")

        v = [ self.q - mk for mk in self.m ]
        a = [ (v[k] + self.m[k-1], self.m[k-2]) for k in range(len(self.m)) ]
        b = [ (v[k] + self.m[k-2], self.m[k-1]) for k in range(len(self.m)) ]
        neighborhood = [ item for sublist in zip(a,b[1:]+b[:1]) for item in sublist ]
        for ray in neighborhood:
            
            v = ray[1]
            x0 = ray[0]
            p0 = self.projector(x0)
            p0p = self.transform(p0)
            x = x0
            p = p0
            while self.plane(x):
                pp = self.transform(p)
                self.canvas.create_circle( pp[0], pp[1], 
                                           self.dotSize, fill="red", outline="red", tags="rays")
                x += v
                p = self.projector(x)
                
            pp = self.transform(p)
            self.canvas.create_circle( pp[0], pp[1], 
                                       self.dotSize+2, outline="red", width=1, tags="rays")

            if x != x0:
                pp = self.transform(p)
                self.canvas.create_line(p0p[0], p0p[1], pp[0], pp[1], 
                                        arrow=tk.LAST, fill="red",width=2, tags="rays")

    def drawBigTriangle(self):
        self.canvas.delete("bigTriangle")

        v = [ self.q - mk for mk in self.m ]
        bigTriangle = [ vk + v[k-1] - v[k-2] for k, vk in enumerate(v) ]
        t = [ self.projector(x) for x in bigTriangle ]
        coords = self.flatten( [ self.transform(x) for x in t] )
        self.triangleItemId = self.canvas.create_polygon(coords, fill="", outline="blue", width=3, tags="bigTriangle")

        for x in bigTriangle:
            p = self.transform(self.projector(x))
            if self.plane(x):
                self.canvas.create_circle( p[0], p[1], 
                                           self.dotSize, fill="blue", outline="blue", tags="bigTriangle")
            else:
                self.canvas.create_circle( p[0], p[1], 
                                           self.dotSize+2, outline="blue", width=1, tags="bigTriangle")

    def drawPiece(self):
        self.canvas.delete("piece")
        for k,tileSet in enumerate(self.tiles):
            for tile in tileSet:
                self.drawTile(tile, self.mode.color(Tile.colorByType[tile.n], Tile.colorByGroup[k]))
                
    def drawTile(self, t, c):
        polygon3d = [ t.origin(),
                      t.origin() + t.v1,
                      t.origin() + t.v1 + t.v2,
                      t.origin() + t.v2 ]
        polygon2d = [ self.projector(x) for x in polygon3d ]
        coords = self.flatten( [ self.transform(x) for x in polygon2d] )
        self.canvas.create_polygon(coords, 
                                   fill=c, outline="black",
                                   width=2, tags="piece")
                                
    def drawBasis(self):

        o = PointVector([0]*3)
        op = self.projector(o)
        e0p = self.projector(self.e0)
        e1p = self.projector(self.e1)
        e2p = self.projector(self.e2)
        opp = self.transform(op)
        
        self.canvas.create_line(self.flatten( [ opp, self.transform(e0p) ] ),
                                arrow=tk.LAST, fill="purple", width=2, tags="basis")
        self.canvas.create_line(self.flatten( [ opp, self.transform(e1p) ] ),
                                arrow=tk.LAST, fill="purple", width=2, tags="basis")
        self.canvas.create_line(self.flatten( [ opp, self.transform(e2p) ] ),
                                arrow=tk.LAST, fill="purple", width=2, tags="basis")

        pointKey = self.canvas.create_circle(opp[0], opp[1],
                                             self.dotSize, fill="purple",tags="basis")
        self.canvas.tag_bind( pointKey, "<Button-1>", self.selectStartingPoint )
        self.grid[pointKey] = o
        
    def updateForwardAndDraw(self,event): 
        
        if self.index < len(self.operations):
            k, alpha, beta = self.operations[self.index]

            #tile k does not change
            #other tiles are either sheared if their coefficient is zero
            #or obtained by concatenation
            if alpha == 1: 
                for i in range(0,beta+1):
                    temp = deepcopy(self.tiles[k])
                    for t in temp:
                        t.shift( -self.m[k] + self.m[k-1] + self.m[k-2]*i )
                    self.tiles[k-1] += temp
            if beta == 1: 
                for i in range(0,alpha+1):
                    temp = deepcopy(self.tiles[k])
                    for t in temp:
                        t.shift( -self.m[k] + self.m[k-1]*i + self.m[k-2] )
                    self.tiles[k-2] += temp


            self.m[k] += - self.m[k-1]*alpha - self.m[k-2]*beta
            self.index += 1

            self.drawTriangle()

    def updateBackwardAndDraw(self,event): 

        if self.index > 0:
            self.index -= 1
            k, alpha, beta = self.operations[self.index]

            if alpha == 1:
                nb = len(self.tiles[k])*(beta+1)
                del self.tiles[k-1][-nb:]
            if beta == 1: 
                nb = len(self.tiles[k])*(alpha+1)
                del self.tiles[k-2][-nb:]
            
            self.m[k] -= - self.m[k-1]*alpha - self.m[k-2]*beta

            self.drawTriangle()
            
    def reSetTriangle(self,event): 

        self.initPattern()        
        self.drawTriangle()

    def setOption(self,event): 
        if event.char == "H": 
            self.enableH = not self.enableH
            if self.enableH:
                self.drawHexagon()
            else:
                self.canvas.delete("hexagon")
        elif event.char == "R": 
            self.enableR = not self.enableR
            if self.enableR:
                self.drawRays()
            else:
                self.canvas.delete("rays")
        elif event.char == "T":
            self.enableT = not self.enableT
            if self.enableT:
                self.drawBigTriangle()
            else:
                self.canvas.delete("bigTriangle")
        elif event.char == "P":
            self.enableP = not self.enableP
            if self.enableP:
                self.drawPiece()
            else:
                self.canvas.delete("piece")

        self.drawTriangle()
                
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
parser.add_argument("-c", "--color",help="either a color by tile type or a color by tile group",\
                    choices=["type","group"],default="type")
parser.add_argument("-s", "--startingCorner",help="the pattern can either contain the upper of the lower corner",\
                    choices=["upper","lower"],default="upper")
parser.add_argument("-k", "--showKeybindings", help="print to the standard output the keys you can hit to modify the display",
                    action="store_true")
   
args = parser.parse_args()
param = [args.x, args.y, args.z]
n = PointVector( param )
projector = standardProjector
if args.projection == "hexagonal":
    projector = hexagonalProjector
mode = PatternMode()
mode.color = firstArg
if args.color == "group":
    mode.color = secondArg
mode.corner = "upper"
if args.startingCorner == "lower":
    mode.corner = args.startingCorner
if args.showKeybindings:
   print("H:enables/disables hexagon display")
   print("R:enables/disables rays display")
   print("T:enables/disables circumscribing triangle display")
   print("P:enables/disables display of the underlying digital plane pattern")
   print("(All displays can be superimposed)")
   print("<Right> next triangle")
   print("<Left> previous triangle")
   print("<Return> come back to the starting triangle")

#application
root = tk.Tk()
c = TriangleComputerApp(root, n, args.windowSize, args.unitSize, projector, mode)     
root.mainloop()
    

