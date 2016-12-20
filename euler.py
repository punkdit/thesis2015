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

st_thick = [style.linewidth.thick]
st_Thick = [style.linewidth.Thick]
st_THick = [style.linewidth.THick]
st_THIck = [style.linewidth.THIck]


text.set(mode="latex")
text.set(docopt="12pt")
#text.preamble(r'\usepackage{amsmath,amsfonts,amssymb}')
text.preamble(r'\usepackage{amsmath,amssymb}')
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
def push(*args):
    global c
    stack.append(c)
    c = canvas.canvas(*args)

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


c = canvas.canvas()




R0 = 4. 
R1 = 1.

X = -6.0
Y = 0.


def torus(u, v, R0=R0, R1=R1):
    x = (R0 + R1*cos(v))*cos(u)
    y = (R0 + R1*cos(v))*sin(u)
    z = R1*sin(v)
    return x, y, z
    

def loop1(u, v0=0., v1=2*pi, N=50):
    v = v0
    dv = (v1-v0)/N
    for i in range(N+1):
        yield torus(u, v)
        v += dv

def loop2(v, u0=0., u1=2*pi, N=50):
    u = u0
    du = (u1-u0)/N
    for i in range(N+1):
        yield torus(u, v)
        u += du


def display(x, y, z):
    return (x-0.2*y-0.0*z+X, 0.8*z-0.3*y+Y)


ni = 60
nj = 60
u = 0.
a = 0.01
du = 2*pi/ni
rects = []
for i in range(ni):
    #ps = [display(*p) for p in loop1(u)]
    #dopath(ps)

    v = 0.
    dv = 2*pi/nj
    for j in range(nj):
        p0 = torus(u-a, v-a)
        p1 = torus(u-a, v+dv+a)
        p2 = torus(u+du+a, v+dv+a)
        p3 = torus(u+du+a, v-a)
        #dopath([p0, p1, p2, p3], fill=[grey])
        cl = 0.3*(cos(-0.4*pi+v+0.0*u)+1.) + 0.15
        rects.append([p0, p1, p2, p3, cl])
        v += dv

    u += du


def order(rect):
    y = sum(p[1] for p in rect[:4])
    return y

rects.sort(key = order)

for rect in rects:
    cl = rect[4]
    cl = rgb(cl, cl, cl)
    rect = [display(*p) for p in rect[:4]]
    dopath(rect, fill=[cl], stroke=False)


#st_front = [green]+st_Thick+[style.linecap.round]
#st_back = [green]#+st_dotted+st_Thick

st_front = [black]+[style.linewidth.Thick]+[style.linecap.round]
st_back = [rgb(0.4,0.4,0.4)]+[style.linewidth.THin]

t0, t1 = -0.30*pi, 0.63*pi
dopath([display(*p) for p in loop1(0.2*pi, t1, t0+2*pi)], extra=st_back)

t0, t1 = 0.4*pi, 1.3*pi
dopath([display(*p) for p in loop1(1.2*pi, t1, t0+2*pi)], extra=st_back)

t0, t1 = -0.05*pi, 0.95*pi
dopath([display(*p) for p in loop2(0., t1, t0+2*pi)], extra=st_back)

t0, t1 = -0.125*pi, 0.995*pi
dopath([display(*p) for p in loop2(pi, t0, t1)], extra=st_back)

t0, t1 = -0.30*pi, 0.63*pi
dopath([display(*p) for p in loop1(0.2*pi, t0, t1)], extra=st_front)

t0, t1 = 0.4*pi, 1.3*pi
dopath([display(*p) for p in loop1(1.2*pi, t0, t1)], extra=st_front)

t0, t1 = -0.05*pi, 0.95*pi
dopath([display(*p) for p in loop2(0., t0, t1)], extra=st_front)

t0, t1 = -0.125*pi, 0.995*pi
dopath([display(*p) for p in loop2(pi, t1, t0+2*pi)], extra=st_front)



c.text(0, 0, "$=$", center)

W = 4.
H = 4.
x = 1.0
y = -0.5*H

c.fill(path.rect(x, y, W, H), [grey])

st_front = [black]+[style.linewidth.Thick]

c.stroke(path.line(x+0.4*W, y, x+0.4*W, y+H), st_front)
c.stroke(path.line(x+0.8*W, y, x+0.8*W, y+H), st_front)

c.stroke(path.line(x, y+0.2*H, x+W, y+0.2*H), st_front)
c.stroke(path.line(x, y+0.6*H, x+W, y+0.6*H), st_front)

for st in [[white]+st_thick,st_dashed+st_thick]:
    c.stroke(path.line(x, y, x+W, y), st)
    c.stroke(path.line(x, y+H, x+W, y+H), st)

for st in [[white]+st_thick,st_dotted+st_thick]:
    c.stroke(path.line(x, y, x, y+H), st)
    c.stroke(path.line(x+W, y, x+W, y+H), st)

c.writePDFfile("pic-torus.pdf")


#############################################################################
#
#


c = canvas.canvas()


W = 2.0
H = 2.0
m = 0.3
r = 0.10


for i in range(2):
    X = [-1.2*W, +1.2*W][i]

    pos = [trafo.translate(X, 0.)]
    
    p = path.rect(0, 0, W, H)
    c.fill(p, [grey]+pos)
    c.stroke(p, pos)
    
    st = [style.linewidth.THick] + pos
    c.stroke(path.line(-m, 0, -m, H), st)
    c.stroke(path.line(W+m, 0, W+m, H), st)
    c.stroke(path.line(0, -m, W, -m), st)
    c.stroke(path.line(0, H+m, W, H+m), st)
    
    st = [] + pos
    c.fill(path.circle(-m, -m, r), st)
    c.fill(path.circle(-m, H+m, r), st)
    c.fill(path.circle(W+m, H+m, r), st)
    c.fill(path.circle(W+m, -m, r), st)
    
    st = [deco.earrow(size=0.2)] + pos

    for theta in [0., 90., 180., 270.]:

        rot = [trafo.rotate(theta, 0.5*W, 0.5*W)]

        if i==0:
            p = path.line(0.75*W, 0.5*H, 1.09*W, 0.5*H)
            c.stroke(p, rot+pos+[
                deco.earrow(size=0.35), 
                style.linewidth.THick, black])
    
            p = path.line(0.77*W, 0.5*H, 1.05*W, 0.5*H)
            c.stroke(p, rot+pos+[deco.earrow(size=0.2), white])

            c.text(0.5*W, -4*m, r"$\partial_2$", center+pos)

        else:

            m1 = 1.6*m
            p = path.line(1.0*W+m1, 0.65*H, 1.0*W+m1, 0.95*H)
            c.stroke(p, rot+pos+[
                deco.earrow(size=0.35), 
                style.linewidth.THick, black])
    
            p = path.line(1.0*W+m1, 0.67*H, 1.0*W+m1, 0.91*H)
            c.stroke(p, rot+pos+[deco.earrow(size=0.2), white])
        
            p = path.line(1.0*W+m1, H-0.65*H, 1.0*W+m1, H-0.95*H)
            c.stroke(p, rot+pos+[
                deco.earrow(size=0.35), 
                style.linewidth.THick, black])
    
            p = path.line(1.0*W+m1, H-0.67*H, 1.0*W+m1, H-0.91*H)
            c.stroke(p, rot+pos+[deco.earrow(size=0.2), white])
    
            c.text(0.5*W, -4*m, r"$\partial_1$", center+pos)



c.writePDFfile("pic-bdy.pdf")


#############################################################################
#
#




W = 4.
H = 4.
x = 0.0
y = 0.0

m = 0.05

c = canvas.canvas([canvas.clip(path.rect(-m, -m, W+2*m, H+2*m))])

w = 0.5*W
h = 0.5*H

x0 = 0.7*w
y0 = 0.3*h

m = 0.18

r = 0.07

st_edge = [style.linewidth.THick]
for i in range(-1, 2):
  for j in range(-1, 2):
    p = path.rect(x0 + i*w + m, y0 + j*h + m, w-2*m, h-2*m)
    c.fill(p, [grey])
    c.stroke(p)

    p = path.line(x0 + i*w + m, y0 + j*h, x0 + (i+1)*w - m, y0+j*h)
    c.stroke(p, st_edge)

    p = path.line(x0 + i*w, y0 + j*h + m, x0 + i*w, y0+(j+1)*h-m)
    c.stroke(p, st_edge)

    c.fill(path.circle(x0 + i*w, y0 + j*h, r))


st_front = [black]+[style.linewidth.Thick]

#c.stroke(path.line(x+0.4*W, y, x+0.4*W, y+H), st_front)
#c.stroke(path.line(x+0.8*W, y, x+0.8*W, y+H), st_front)
#
#c.stroke(path.line(x, y+0.2*H, x+W, y+0.2*H), st_front)
#c.stroke(path.line(x, y+0.6*H, x+W, y+0.6*H), st_front)


# border:

m = 0.03
for st in [[white]+st_THIck, st_dashed+st_thick]:
    c.stroke(path.line(x-m, y-m, x+W+m, y-m), st)
    c.stroke(path.line(x-m, y+H+m, x+W-m, y+H+m), st)

for st in [[white]+st_THIck, st_dotted+st_thick]:
    c.stroke(path.line(x-m, y-m, x-m, y+H+m), st)
    c.stroke(path.line(x+W+m, y+m, x+W+m, y+H+m), st)



c.writePDFfile("pic-torus-hom.pdf")



