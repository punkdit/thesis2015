#!/usr/bin/env python3

import sys, os

from numpy import array, dot, allclose, identity


H = array([
    [3, 1, 1, 1, 0, 0, 0, 0], 
    [1, 1, 0, 0, 1, 1, 0, 0], 
    [1, 0, 1, 0, 1, 0, 1, 0], 
    [1, 0, 0, 1, 0, 1, 1, 0], 
    [0, 1, 1, 0, -1, 0, 0, 1], 
    [0, 1, 0, 1, 0, -1, 0, 1], 
    [0, 0, 1, 1, 0, 0, -1, 1], 
    [0, 0, 0, 0, 1, 1, 1, -3]])

E = array([
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1]])

E1 = array([
    [1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, -1, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, -1, 1, 0, 0, 0, 1]])

F = array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 1, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 1, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]])

assert allclose(dot(E, F), identity(4))
assert allclose(dot(E1, F), identity(4))

HA = dot(dot(E, H), F)

print(HA)

assert allclose(HA, dot(dot(E1, H), F))

HAt = HA.transpose()
print("normal:", allclose(dot(HAt, HA), dot(HA, HAt)))



