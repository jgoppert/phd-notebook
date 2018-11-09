import numpy as np
from math import pi, sin


def dyn_const(x_array, t, A):
    x = np.matrix(x_array).T
    xdot = A*x
    return np.array(xdot.T)[0]


def dyn_switch(x_array, t, A1, A2, period):
    x = np.matrix(x_array).T
    if sin((2*pi/period)*t) > 0:
        xdot = A1*x
    else:
        xdot = A2*x
    return np.array(xdot.T)[0]

# vi:ts=4:sw=4:expandtab
