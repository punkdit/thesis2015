#!/usr/bin/env python

import sys
import os

from math import *
from random import *

from pyx import canvas, path, deco, trafo, style, text, color, unit, epsfile, deformer, bitmap

from PIL import Image

north = [text.halign.boxcenter, text.valign.top]
northeast = [text.halign.boxright, text.valign.top]
northwest = [text.halign.boxleft, text.valign.top]
south = [text.halign.boxcenter, text.valign.bottom]
southeast = [text.halign.boxright, text.valign.bottom]
southwest = [text.halign.boxleft, text.valign.bottom]
east = [text.halign.boxright, text.valign.middle]
west = [text.halign.boxleft, text.valign.middle]
center = [text.halign.boxcenter, text.valign.middle]


st_dashed = [style.linestyle.dashed]
st_dotted = [style.linestyle.dotted]

st_Thick = [style.linewidth.Thick]


text.set(mode="latex")
text.set(docopt="10pt")
text.preamble(r'\usepackage{amsmath,amsfonts,amssymb}')
#text.preamble(r"\def\I{\mathbb{I}}")
text.preamble(r"\def\ket #1{|#1\rangle}")


rgb = color.rgb
rgbfromhexstring = color.rgbfromhexstring

red, green, blue, yellow, orange = (
    rgbfromhexstring("#d00000"),
    rgbfromhexstring("#006000"),
    rgb.blue, 
    rgb(0.75, 0.75, 0),
    rgb(0.75, 0.55, 0),
    )

blue = rgb(0., 0., 0.8)
pale_blue = rgb(0.7, 0.7, 1.0)
pink = rgb(1., 0.4, 0.4)
white = rgb(1., 1., 1.)
black = rgb(0., 0., 0.)
grey = rgb(0.8, 0.8, 0.8)

brown = rgbfromhexstring("#AA6C39"),

light_shade = rgb(0.85, 0.65, 0.1)

#shade = brown
#shade = orange
shade = rgb(0.8, 0.8, 0.8)


N = 20

def dopath(ps, extra=[], fill=[], closepath=False, smooth=0.0, stroke=True):
    if not ps:
        print "dopath: empty"
        return
    ps = [path.moveto(*ps[0])]+[path.lineto(*p) for p in ps[1:]]
    if closepath:
        ps.append(path.closepath())
    p = path.path(*ps)
    extra = list(extra)
    if smooth:
        extra.append(deformer.smoothed(smooth))
    if fill:
        c.fill(p, extra+fill)
    if stroke:
        c.stroke(p, extra)


stack = []
def push():
    global c
    stack.append(c)
    c = canvas.canvas()

def pop(*args):
    global c
    c1 = stack.pop()
    c1.insert(c, *args)
    c = c1



class Turtle(object):
    def __init__(self, x, y, theta):
        self.x = x
        self.y = y
        self.theta = theta
        self.ps = [(x, y)]
        self.pen = True

    def penup(self):
        self.pen = False
        self.ps = []
        return self

    def pendown(self):
        self.pen = True
        self.ps = [(self.x, self.y)]
        return self

    def fwd(self, d):
        self.x += d*sin(self.theta)
        self.y += d*cos(self.theta)
        if self.pen:
            self.ps.append((self.x, self.y))
        return self

    def reverse(self, d):
        self.fwd(-d)
        return self

    def right(self, dtheta, r=0.):
        theta = self.theta
        self.theta += dtheta
        if r==0.:
            return self
        N = 20
        x, y = self.x, self.y
        x0 = x - r*sin(theta-pi/2)
        y0 = y - r*cos(theta-pi/2)
        for i in range(N):
            theta += (1./(N))*dtheta
            x = x0 - r*sin(theta+pi/2)
            y = y0 - r*cos(theta+pi/2)
            if self.pen:
                self.ps.append((x, y))
        self.x = x
        self.y = y
        return self

    def left(self, dtheta, r=0.):
        self.right(-dtheta, -r)
        return self

    def stroke(self, extra=[], fill=[], closepath=False):
        dopath(self.ps, extra, fill, closepath, smooth=0.)
        self.ps = [(self.x, self.y)]
        return self



#############################################################################
#
#

c = canvas.canvas()


def blob(x, y, R):
    #t = Turtle(0., 0., 0.)
    p = path.circle(x, y, R)
    c.fill(p, [color.transparency(0.8)])
    c.stroke(p)


blob(-0.5, 0.5, 1.0)
blob(-0.5, -0.5, 1.0)
blob(0.5, 0., 1.0)

c.text(-0.8, 0.8, "$A$", center)
c.text(-0.8, -0.8, "$B$", center)
c.text(0.8, 0.0, "$C$", center)

c.writePDFfile("pic-ABC.pdf")




#############################################################################
#
#

c = canvas.canvas()

R = 1.0


def arc(theta0, theta1, radius, N=50):
    theta = theta0
    dtheta = (theta1-theta0)/N
    ps = []
    for i in range(N+1):
        x = radius*cos(theta)
        y = radius*sin(theta)
        ps.append((x, y))
        theta += dtheta
    return ps


def display(x, y, z):
    return (x-0.1*y-0.2*z, 0.7*z-0.8*y)


def render(ps):
    for tx in [
        (lambda x, y : (x, y, y)),
        (lambda x, y : (x, y, -y)),
        (lambda x, y : (y, x, y)),
        (lambda x, y : (y, x, -y)),
        (lambda x, y : (y, -y, x)),
        (lambda x, y : (y, y, x)),
        ]:
        ps1 = [tx(*p) for p in ps]
        ps1 = [display(*p) for p in ps1]
    
        dopath(ps1)
    

#ps = arc(pi/4, 3*pi/4., R)
#render(ps)
#ps = arc(5*pi/4, 7*pi/4., R)
#render(ps)
#
#
#c.writePDFfile("pic-sphere.pdf")


#############################################################################
#
#

c = canvas.canvas()

R = 1.7

X0 = -3.6*R
TX = trafo.translate(X0, 0)
rect = [(-1., -1.), (-1., 1.), (1., 1.), (1., -1.)]

c.fill(path.circle(X0, 0, 1.5*R), [rgb(0.7,0.7,0.7)])

#ps = [(X0+x, y) for (x,y) in rect]
#dopath(ps, closepath=True, smooth=0.2)


def square(x, y, R1, extra=[]):
    theta = 0.1*pi
    vs = [(x, y)]
    t = Turtle(x, y, theta)
    t.left(2*theta, R1).left(0.5*pi-2*theta)
    vs.append((t.x, t.y))
    t.left(2*theta, R1).left(0.5*pi-2*theta)
    vs.append((t.x, t.y))
    t.left(2*theta, R1).left(0.5*pi-2*theta)
    vs.append((t.x, t.y))
    t.left(2*theta, R1).left(0.5*pi-2*theta)
    t.stroke(extra+[TX])
    return vs


lshade = rgb(0.4,0.4,0.4)

x, y = 1.3, -0.7
vs = square(x, y, 2.1*R, [lshade])
print vs

st_back = [TX, lshade]
st_front = [TX]+st_Thick

x, y = vs[0]
t = Turtle(x, y, 0.7*pi)
t.right(0.3,2.)
t.right(0.34,1.2)
t.right(0.5,0.4)
t.stroke(st_back)
t.right(1.0)
t.right(0.4, 0.4)
t.right(0.3, 0.9)
t.right(0.3, 1.1)
t.right(0.25, 1.3)
t.stroke(st_front)

x, y = vs[1]
t = Turtle(x, y, 0.27*pi)
t.fwd(0.555)
t.stroke(st_back)
t.right(1.00*pi)
t.fwd(1.0)
t.stroke(st_front)

x, y = vs[2]
t = Turtle(x, y, -0.27*pi)
t.left(0.3, 1.55)
t.left(0.5, 0.8)
t.left(0.5, 0.3)
t.stroke(st_back)
t.left(0.8)
t.left(0.5, 0.32)
t.left(0.65, 0.8)
t.stroke(st_front)

x, y = vs[3]
t = Turtle(x, y, -0.7*pi)
#t.fwd(1.4)
t.left(0.178, 8.0)
t.stroke(st_back)
t.left(0.93*pi)
#t.fwd(0.57)
t.left(0.07, 8.0)
t.stroke(st_front)



x, y = 1.0, -1.3
vs = square(x, y, 2.4*R, st_Thick)
print vs

c.text(-1.7*R, 0., "$=$", center)


def display(x, y, z):
    return (x-0.2*y-0.0*z, 0.8*z-0.2*y)


back = [
    lambda x, y : (x, y, -1), # bot
    lambda x, y : (-1, y, x), # left
    lambda x, y : (x, -1, y), # back
]
front = [
    lambda x, y : (x, y, 1), # top
    lambda x, y : (1, y, x), # right
    lambda x, y : (x, 1, y), # front
]
faces = back+front

for face in faces:
    
    ps1 = [face(x, y) for (x,y) in rect]
    ps1 = [(x*R, y*R, z*R) for (x,y,z) in ps1]
    ps1 = [display(*p) for p in ps1]

    dopath(ps1, closepath=True, fill=[color.transparency(0.8)])

    if face in back:
        st = [rgb(0.4,0.4,0.4)]
    else:
        st = st_Thick+[style.linejoin.round]
    dopath(ps1, closepath=True, extra=st)


c.writePDFfile("pic-cube.pdf")


#############################################################################
#
#


sys.exit(0)
c = canvas.canvas()

W = 5.
H = 3.

c.stroke(path.line(0.0, -0.2, 0., H), [deco.earrow()])
c.stroke(path.line(-0.2, 0., W, 0), [deco.earrow()])
my = 0.25
c.stroke(path.line(-0.2, -my*0.2, 0.9*W, my*W), [deco.earrow()])

x0, y0 = 0.3*W, 0.3*H
x1, y1 = x0+0.3*W, y0+my*0.3*W


#t = Turtle(x0, y0, pi/4)
#t.penup().fwd(0.2).pendown()
#theta, r = 0.25*pi, 2.0
#t.right(theta, r).right(0.8*pi, 0.2).right(theta, r)
#t.fwd(0.1).stroke([deco.earrow()])

tr = trafo.scale(x=0.5*W, y=0.3*H, sx=1.3, sy=1.5)

c.fill(path.circle(x0, y0, 0.06), [tr])
c.fill(path.circle(x1, y1, 0.06), [tr])

def loop(x0, y0, r1, theta0, theta1, tpy):
    t = Turtle(x0, y0, theta0)
    #theta = 0.55*pi
    theta = theta1
    r2 = 2*r1*sin(theta - 0.5*pi)
    t.penup().fwd(0.1*r1).pendown()
    t.fwd(0.9*r1).right(theta).fwd(r2).right(theta).fwd(0.9*r1) 
    t.stroke([deco.earrow(), 
        deformer.smoothed(4.0), style.linewidth.Thick,
        color.transparency(tpy), tr])

r1 = 2.0
tpy = 0.2
theta = 0.37*pi
theta1 = 0.55*pi
for i in range(10):
    tpy = 0.8
    if i==0 or i==9:
        tpy = 0.
    loop(x0, y0, r1, theta, theta1, tpy)
    #tpy += 0.05
    theta -= 0.06*pi
    #r1 = 0.9*r1
    r1 -= 0.17
    theta1 += 0.01*pi

#c.writePDFfile("pic-monodromy3d.pdf")



