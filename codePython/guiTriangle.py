import tkinter as tk
import argparse

from copy import deepcopy

from NormalComputer.PointVector import PointVector
from NormalComputer.DigitalPlane import DigitalPlane, generateBFS, isReduced, additiveReduction
from NormalComputer.TriangleComputer import TriangleComputer, squaredRadiusOfSphere

#------------------------------------------------
#----------------- helpers ----------------------
#------------------------------------------------

#hexagonalProjection = ( (-0.5,-0.86602540378), (1,0), (-0.5, 0.86602540378) )

def hexagonalProjector(vector3):
    vector2 = []
    vector2.append( -0.5 * vector3[0] + vector3[1] - 0.5 * vector3[2] )
    vector2.append( 0.86602540378 * vector3[0] - 0.86602540378 * vector3[2]) #NB. we take the opposite for the y-coordinate
    return tuple(vector2)
    
#standardProjection = ( (-0.353,-0.353), (1,0), (0,1) )

def standardProjector(vector3):
    vector2 = []
    vector2.append( -0.353 * vector3[0] + vector3[1] )
    vector2.append( 0.353 * vector3[0] - vector3[2] ) #NB. we take the opposite for the y-coordinate
    return tuple(vector2)

def _create_circle(self, x, y, r, **kwargs):
    return self.create_oval(x-r, y-r, x+r, y+r, **kwargs)
tk.Canvas.create_circle = _create_circle

class Parallelogram(object):

    def __init__(self, origin, v1, v2, color):
        self.o = origin
        self.v1 = v1
        self.v2 = v2
        self.c = color

    def shift(self, v):
        self.o += v

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
    
class TriangleComputerApp(object):
    """ Simple class for a tkinter application that draws the normal computer output onto a digital plane """

    def embbed(self,x,k):
        """ returns the embedding of the k-th component of point x"""
        return (x[k]+self.origin[k])*self.gridStep
    
    def __init__(self, parent, normal, size, step, projector):
        """ Initilization of the application

        :param parent: object of type Tk
        :param normal: plane normal (provided as a PointVector)
        :param size: discrete size of the drawing window
        :param step: grid step of the drawing window
        :param projector: functor that takes a 3d PointVector and
        returns its projection as a 2-tuple. 
        """
        #dimensions
        self.discreteSize = size
        self.gridStep = step
        self.dotSize = step / 10.0

        #origin
        mid = self.discreteSize/2
        self.origin = (mid,mid)
        
        #canvas
        realSize = self.gridStep * self.discreteSize
        self.canvas = tk.Canvas(parent, width=realSize, height=realSize, borderwidth=0, highlightthickness=0)
        self.canvas.pack()

        #projector
        self.projector = projector
        
        #digital plane
        o = PointVector([0] * 3)
        self.plane = DigitalPlane(normal)
        pointSet3 = set()
        edgeSet3 = set()
        generateBFS(o, lambda x: self.plane(x) and x.norm("Linf") <= size,  pointSet3, edgeSet3)
           
        #digital plane as dictionnary of 3d points
        for e3 in edgeSet3:
            base = self.projector(e3[0].toTuple())
            tip = self.projector(e3[1].toTuple())
            self.canvas.create_line(self.embbed(base,0), self.embbed(base,1),
                                    self.embbed(tip,0), self.embbed(tip,1),
                                    fill="gray", width=2)

        self.grid = {}
        for p3 in pointSet3:
            p2 = self.projector(p3.toTuple())
            pointKey = self.canvas.create_circle( self.embbed(p2,0),
                                                  self.embbed(p2,1),
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

        #draw basis
        self.drawBasis()

        #data
        self.q = PointVector([1]*3)
        self.start()

    #-------------------------------------------------------------------------
    #-------------------------------- actions --------------------------------
    #-------------------------------------------------------------------------

    def start(self):

        #normal computer to get the list of operations
        s = PointVector([1]*3)
        o = self.q - s
        e0 = PointVector([1,0,0])
        e1 = PointVector([0,1,0])
        e2 = PointVector([0,0,1])
        nc = TriangleComputer([o+e1+e2, o+e2+e0, o+e0+e1], self.q, s, self.plane )
        while nc.advance():
            pass

        self.operations = nc.operations

        #other data
        self.index = 0
        self.m = [e0, e1, e2]
        self.tiles = [ [ Parallelogram(self.q, -self.m[1], -self.m[2], "gray60") ],
                       [ Parallelogram(self.q, -self.m[2], -self.m[0], "white") ],
                       [ Parallelogram(self.q, -self.m[0], -self.m[1], "gray30") ] ]
        
        #draw first triangle
        self.drawTriangle()
        
        
    def selectStartingPoint(self,event):
        report_event(event)
        pointKey = self.canvas.find_withtag(tk.CURRENT)[0] #current item id
        self.q = self.grid[pointKey] + PointVector([1]*3)
        self.start()
        
    def drawTriangle(self): 
        self.canvas.delete("triangle")
        
        if self.enableP:
            self.drawPiece()

        #draw reference points
        ref = self.projector(self.q.toTuple())
        pointKey = self.canvas.create_circle( self.embbed(ref,0),
                                              self.embbed(ref,1),
                                              self.dotSize, fill="blue",tags="triangle")

        v = [ self.q - mk for mk in self.m ]
        t = [ self.projector(x.toTuple()) for x in v ]
        self.triangleItemId = self.canvas.create_line(self.embbed(t[0],0), self.embbed(t[0],1),
                                                      self.embbed(t[1],0), self.embbed(t[1],1),
                                                      self.embbed(t[2],0), self.embbed(t[2],1), 
                                                      self.embbed(t[0],0), self.embbed(t[0],1), 
                                                      fill="blue",width=3, tags="triangle")
                
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
        hexagon = [ item for sublist in zip(a,b[1:]+b[:1]) for item in sublist ]
        for x in hexagon:
            p = self.projector(x.toTuple())
            if self.plane(x):
                self.canvas.create_circle( self.embbed(p,0),
                                           self.embbed(p,1),
                                           self.dotSize, fill="red", outline="red", tags="hexagon")
            else:
                self.canvas.create_circle( self.embbed(p,0),
                                           self.embbed(p,1),
                                           self.dotSize+2, outline="red", width=1, tags="hexagon")

        h = [ self.projector(x.toTuple()) for x in hexagon ]
        self.canvas.create_line(self.embbed(h[0],0), self.embbed(h[0],1),
                                self.embbed(h[1],0), self.embbed(h[1],1),
                                self.embbed(h[2],0), self.embbed(h[2],1), 
                                self.embbed(h[3],0), self.embbed(h[3],1), 
                                self.embbed(h[4],0), self.embbed(h[4],1),
                                self.embbed(h[5],0), self.embbed(h[5],1), 
                                self.embbed(h[0],0), self.embbed(h[0],1), 
                                fill="red",width=2, tags="hexagon")


    def drawRays(self):
        self.canvas.delete("rays")

        v = [ self.q - mk for mk in self.m ]
        a = [ (v[k] + self.m[k-1], self.m[k-2]) for k in range(len(self.m)) ]
        b = [ (v[k] + self.m[k-2], self.m[k-1]) for k in range(len(self.m)) ]
        neighborhood = [ item for sublist in zip(a,b[1:]+b[:1]) for item in sublist ]
        for ray in neighborhood:
            
            v = ray[1]
            x0 = ray[0]
            p0 = self.projector(x0.toTuple())
            x = x0
            p = p0
            while self.plane(x):
                self.canvas.create_circle( self.embbed(p,0),
                                           self.embbed(p,1),
                                           self.dotSize, fill="red", outline="red", tags="rays")
                x += v
                p = self.projector(x.toTuple())
                
            self.canvas.create_circle( self.embbed(p,0),
                                       self.embbed(p,1),
                                       self.dotSize+2, outline="red", width=1, tags="rays")

            if x != x0:
                self.canvas.create_line(self.embbed(p0,0), self.embbed(p0,1),
                                        self.embbed(p,0), self.embbed(p,1), 
                                        arrow=tk.LAST, fill="red",width=2, tags="rays")

    def drawBigTriangle(self):
        self.canvas.delete("bigTriangle")

        v = [ self.q - mk for mk in self.m ]
        bigTriangle = [ vk + v[k-1] - v[k-2] for k, vk in enumerate(v) ]
        t = [ self.projector(x.toTuple()) for x in bigTriangle ]

        self.triangleItemId = self.canvas.create_line(self.embbed(t[0],0), self.embbed(t[0],1),
                                                      self.embbed(t[1],0), self.embbed(t[1],1),
                                                      self.embbed(t[2],0), self.embbed(t[2],1), 
                                                      self.embbed(t[0],0), self.embbed(t[0],1), 
                                                      fill="blue",width=3, tags="bigTriangle")

        for x in bigTriangle:
            p = self.projector(x.toTuple())
            if self.plane(x):
                self.canvas.create_circle( self.embbed(p,0),
                                           self.embbed(p,1),
                                           self.dotSize, fill="blue", outline="blue", tags="bigTriangle")
            else:
                self.canvas.create_circle( self.embbed(p,0),
                                           self.embbed(p,1),
                                           self.dotSize+2, outline="blue", width=1, tags="bigTriangle")

    def drawPiece(self):
        self.canvas.delete("piece")

        for bigTile in self.tiles:
            for tile in bigTile:
                #if tile.o != self.q:
                self.drawParallelogram( tile )
                
    def drawParallelogram(self, parallelogram):
        p = [ parallelogram.o,
              parallelogram.o + parallelogram.v1,
              parallelogram.o + parallelogram.v1 + parallelogram.v2,
              parallelogram.o + parallelogram.v2 ]
        p2 = [ self.projector(x.toTuple()) for x in p ]
        self.canvas.create_polygon(self.embbed(p2[0],0), self.embbed(p2[0],1),
                                   self.embbed(p2[1],0), self.embbed(p2[1],1),
                                   self.embbed(p2[2],0), self.embbed(p2[2],1), 
                                   self.embbed(p2[3],0), self.embbed(p2[3],1), 
                                   fill=parallelogram.c, outline="black",
                                   width=2, tags="piece")
                                
    def drawBasis(self):

        o = PointVector([0]*3)
        e0 = PointVector([1,0,0])
        e1 = PointVector([0,1,0])
        e2 = PointVector([0,0,1])
        op = self.projector(o.toTuple())
        e0p = self.projector(e0.toTuple())
        e1p = self.projector(e1.toTuple())
        e2p = self.projector(e2.toTuple())

        self.canvas.create_line(self.embbed(op,0),self.embbed(op,1),
                                self.embbed(e0p,0),self.embbed(e0p,1),
                                arrow=tk.LAST, fill="purple", width=2, tags="basis")
        self.canvas.create_line(self.embbed(op,0),self.embbed(op,1),
                                self.embbed(e1p,0),self.embbed(e1p,1),
                                arrow=tk.LAST, fill="purple", width=2, tags="basis")
        self.canvas.create_line(self.embbed(op,0),self.embbed(op,1),
                                self.embbed(e2p,0),self.embbed(e2p,1),
                                arrow=tk.LAST, fill="purple", width=2, tags="basis")

        pointKey = self.canvas.create_circle(self.embbed(op,0), self.embbed(op,1),
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

        e0 = PointVector([1,0,0])
        e1 = PointVector([0,1,0])
        e2 = PointVector([0,0,1])
        self.m = [e0, e1, e2]
        self.tiles = [ [ Parallelogram(self.q, -self.m[1], -self.m[2], "gray60") ],
                       [ Parallelogram(self.q, -self.m[2], -self.m[0], "white") ],
                       [ Parallelogram(self.q, -self.m[0], -self.m[1], "gray30") ] ]
        self.index = 0
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
parser.add_argument("-s", "--showKeybindings", help="print to the standard output the keys you can hit to modify the display",
                    action="store_true")
   
args = parser.parse_args()
param = [args.x, args.y, args.z]
n = PointVector( param )
projector = standardProjector
if args.projection == "hexagonal":
    projector = hexagonalProjector
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
c = TriangleComputerApp(root, n, args.windowSize, args.unitSize, projector)     
root.mainloop()
    

