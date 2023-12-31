import sys
import numpy as np
import sympy as sp
from sympy import *

x = symbols('x')


def g(f, xi):
    # f is of type symbols
    x = symbols('x')
    d = diff(f, x)
    xip1 = xi - f.evalf(subs={'x': xi}) / float(d.evalf(subs={'x': xi}))
    print(f.evalf(subs={'x': xi}))
    return xip1


def newton_raphson_solver(f, x0, es, iter_max):
    xip1 = x0
    i = 0
    ea = sys.maxsize
    while i < iter_max and ea > es:
        xi = xip1
        xip1 = g(f, xi)
        if not(xip1 == 0):
            ea = abs(((xip1 - xi) / xip1) * 100)
        i += 1
    print("ea", ea)
    return xip1


fx = x ** 2
print(newton_raphson_solver(fx, 2, 0.0005, 500))
