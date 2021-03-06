#!/usr/bin/env python

import sys
from math import *
from random import *
import numpy

from pyx import canvas, path, deco, trafo, style, text, color, deformer

text.set(mode="latex") 
text.set(docopt="12pt")
text.preamble(r"\usepackage{amsmath,amsfonts,amssymb}")

rgb = color.rgb
rgbfromhexstring = color.rgbfromhexstring

red, green, blue, yellow = (
    rgbfromhexstring("#d00000"),
    rgbfromhexstring("#006000"),
    rgb.blue, rgb(0.75, 0.75, 0)) 

black = rgb(0., 0., 0.) 
blue = rgb(0., 0., 0.8)
lred = rgb(1., 0.4, 0.4)
white = rgb(1., 1., 1.) 

shade = rgb(0.75, 0.55, 0)
grey = rgb(0.75, 0.75, 0.75)

light_shade = rgb(0.85, 0.65, 0.1)
light_shade = rgb(0.9, 0.75, 0.4)


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

def arrow(x0, y0, x1, y1, extra=[]):
    c.stroke(path.line(x0, y0, x1, y1),
        extra+[deco.earrow(size=0.1)])

#c.stroke(
#    path.curve(0, 0, 0, 4, 2, 4, 3, 3),
#    [style.linewidth.THICK, style.linestyle.dashed, color.rgb.blue,
#    deco.earrow([deco.stroked([color.rgb.red, style.linejoin.round]),
#                       deco.filled([color.rgb.green])], size=1)])

def line(x0, y0, x1, y1, extra=[]):
    c.stroke(path.line(x0, y0, x1, y1), extra)

def varrow(x, y0, y1, extra=[], label=None):
    arrow(x, y0, x, y1, extra)
    if label:
        c.text(x-0.1, (y0+y1)/2., label, [text.halign.boxright])

def harrow(x0, x1, y, extra=[], label=None):
    arrow(x0, y, x1, y, extra)
    if label:
        c.text((x0+x1)/2., y+0.1, label,
            [text.valign.bottom, text.halign.boxcenter])


def anyon(x, y, r=0.07):
    c.fill(path.circle(x, y, r), [white])
    c.stroke(path.circle(x, y, r), [style.linewidth.thick])


N = 20

def dopath(ps, extra=[], fill=False, closepath=True):
    ps = [path.moveto(*ps[0])]+[path.lineto(*p) for p in ps[1:]]
    if closepath:
        ps.append(path.closepath())
    p = path.path(*ps)
    if fill:
        c.fill(p, [deformer.smoothed(0.3)]+extra+[color.rgb.white])
    c.stroke(p, [deformer.smoothed(0.3)]+extra)

def ellipse(x0, y0, rx, ry, extra=[], fill=False):
    ps = []
    for i in range(N):
        theta = 2*pi*i / (N-1)
        ps.append((x0+rx*sin(theta), y0+ry*cos(theta)))
    dopath(ps, extra, fill)


#############################################################################
#
#

#sys.exit(0)
#text.set(docopt="10pt")

dashed = []
dotted = [style.linestyle.dashed]

c = canvas.canvas()

old_c = None
def push():
    global c, old_c
    assert old_c is None
    old_c = c
    c = canvas.canvas()
    
def pop(extra=[]):
    global c, old_c
    old_c.insert(c, extra)
    c = old_c
    old_c = None


# --------------------------------------------------------------------

w, h = 1.0, 1.0
m = 0.1

#x0, y0 = 0.6, h + 10*m
x0, y0 = 0., 0.

push()
#c.text(x0-0.8, y0, "(a)")

c.stroke(path.rect(x0, y0, w, h), dashed)
c.stroke(path.rect(x0+w+m, y0, w, h), dotted)

grarrow = [green, style.linewidth.THick, deco.earrow(size=0.2)]

y = y0+0.3*h
c.stroke(path.line(x0+0.5*w, y, x0+1.5*w+m, y), grarrow)
anyon(x0+0.8*w, y)
anyon(x0+m+1.2*w, y)
    
y = y0+0.7*h
c.stroke(path.line(x0+1.5*w+m, y, x0+0.5*w, y), grarrow)
anyon(x0+0.8*w, y)
anyon(x0+m+1.2*w, y)

pop([trafo.rotate(-90), trafo.translate(0.7, 2.5*h)])
c.text(0., 0.5*h, "(a)")

# --------------------------------------------------------------------

#x0, y0 = 0.6, 0.
    
#c.text(x0-0.8, y0, "(b)")
push()

c.stroke(path.rect(x0, y0, w, h), dashed)
c.stroke(path.rect(x0+w+m, y0, w, h), dotted)


p = path.path(
    path.moveto(x0+0.5*w, y0-0.3*h), 
    path.lineto(x0+0.5*w, y0+0.3*h), 
    path.lineto(x0+1.5*w+m, y0+0.3*h),
    path.lineto(x0+1.5*w+m, y0+0.7*h), 
    path.lineto(x0+0.5*w, y0+0.7*h),
    path.lineto(x0+0.5*w, y0+1.3*h),
)
c.stroke(p, [deformer.smoothed(0.3)]+grarrow)

y = y0+0.3*h
anyon(x0+0.8*w, y)
anyon(x0+m+1.2*w, y)
    
y = y0+0.7*h
anyon(x0+0.8*w, y)
anyon(x0+m+1.2*w, y)

pop([trafo.rotate(-90), trafo.translate(2.*w+0.7, 2.5*h)])
c.text(2.*w, 0.5*h, "(b)")


# --------------------------------------------------------------------

push()
w1, h1 = w, h
w, h = 2., 2.
#w, h = 1.5, 1.5

#x0, y0 = 2.8*w + 0*m, -6*m

dx = 0.8*w
dy = 12*m

r = 0.5*h
ys = [0.2*h, 0.7*h, 1.2*h, 1.7*h]


def braid(x0, y0, x1, y1):
    p = path.line(x0, y0, x1, y1)
    c.stroke(p, [color.rgb.white, style.linewidth.THICk])
    c.stroke(p, [red])

def braid(x0, y0, x1, y1):
    x, y = x0, y0
    dx = (x1-x0)/N
    dy = (y1-y0)/N
    ps = [(x, y)]
    for i in range(N):
        theta = i*2*pi/(N)
        x += dx*(cos(theta)+2) / (2.)
        y += dy
        ps.append((x, y))
    dopath(ps, [color.rgb.white, style.linewidth.THICK], closepath=False)
    dopath(ps, [red], closepath=False)



x1, y1 = x0+dx, y0+dy
perm = [2, 0, 1, 3]
for i in [0, 1, 2, 3]:
    braid(x0, y0+ys[i], x1, y1+ys[perm[i]])

#c.text(x0-1.5, 0, "(c)")


y2 = y0+0.95*h
ps = []
for i in range(N):
    theta = 2*pi*i / (N-1)
    ps.append((x0+0.4*r*sin(theta), y2+0.9*r*cos(theta)))

dopath(ps, dotted, fill=True)

ps = []
for i in range(N):
    theta = pi + pi*i / (N-1)
    ps.append((x0+0.5*h*sin(theta), y2+0.9*h*cos(theta)))

y1 = ys[-1]+y0
r1 = ps[-1][1]-y1
for i in range(N):
    theta = pi*i / (N-1)
    ps.append((x0+r1*sin(theta), y1+r1*cos(theta)))

r2 = ps[-1][1]-y2
for i in range(N):
    theta = pi - pi*i / (N-1)
    ps.append((x0-0.5*r2*sin(theta), y2-r2*cos(theta)))

y1 = ys[0]+y0
#r1 = ps[-1][1]-y1
for i in range(N):
    theta = pi*i / (N-1)
    ps.append((x0+r1*sin(theta), y1+r1*cos(theta)))

dopath(ps, dashed, fill=True)

c.stroke(path.line(x0, y0, x0, y0+2*h), grarrow)
for y in ys:
    anyon(x0, y0+y)


# --------------------------------------------------------------------

#c.text(x0+1.2, 0.3*dy, "(d)")

x0 += dx
y0 += dy

c.stroke(path.line(x0, y0, x0, y0+2*h), grarrow)

for y in ys:
    anyon(x0, y0+y)

y2 = y0+0.45*h
ps = []
for i in range(N):
    theta = 2*pi*i / (N-1)
    ps.append((x0+0.4*r*sin(theta), y2+0.8*r*cos(theta)))

dopath(ps, dotted)

y2 = y0+1.45*h
ps = []
for i in range(N):
    theta = 2*pi*i / (N-1)
    ps.append((x0+0.4*r*sin(theta), y2+0.8*r*cos(theta)))

dopath(ps, dashed)

#pop([trafo.rotate(-90)])
pop([trafo.rotate(-90), trafo.translate(2.*w+1.0, 2*h1)])
c.text(4.3*w1, 2*h1, "(c)")
c.text(5.4*w1, 0.5*h1, "(d)")

c.writePDFfile("pic-syndrome.pdf")

