import sys
import numpy as np
import sympy as sp
from sympy import *
import math


def g(f, xi, m):
    d = diff(f, x)
    xip1 = xi - m * (f.evalf(subs={'x': xi}) / float(d.evalf(subs={'x': xi})))
    return xip1


def newton_raphson_solver(f, x0, es, iter_max, n_sig, m):
    xip1 = x0
    i = 0
    ea = sys.maxsize
    while i < iter_max and ea > es:
        xi = xip1
        xip1 = g(f, xi, m)
        xip1 = round(xip1, -int(math.floor(math.log10(abs(xip1)))) + (n_sig - 1))
        print(i, xip1)
        if not (xip1 == 0):
            ea = round(abs(((xip1 - xi) / xip1) * 100), 5)
            print("ea of i =", i, "is equal ", ea)
        i += 1
    print("ea final = ", ea)
    return xip1


x = symbols('x')
fx = x ** 3 - 0.165 * (x ** 2) + 3.993E-4
print(newton_raphson_solver(fx, 0.05, 0.05, 10, 5, 1))
