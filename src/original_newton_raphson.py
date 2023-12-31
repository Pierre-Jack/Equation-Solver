import sys
import numpy as np
import sympy as sp
from sympy import *
import math


def check_original_newton_raphson(f):
    # ans_array = sp.solve(f, x)
    # print(ans_array)
    # np_ans_array = np.array(ans_array)
    # ans = np.max(np_ans_array)
    # print(ans)
    # f1 = diff(f, x)
    # f2 = diff(f1, x)
    # return True if abs(f2.evalf(subs={'x': ans})/float(2 * f1.evalf(subs={'x': ans}))) < 1 else False
    return true


def g(f, xi):
    # f is of type symbols
    d = diff(f, x)
    xip1 = xi - f.evalf(subs={'x': xi}) / float(d.evalf(subs={'x': xi}))
    return xip1


def newton_raphson_solver(f, x0, es, iter_max, n_sig):
    if not check_original_newton_raphson(f):
        return "not solvable"
    xip1 = x0
    i = 0
    ea = sys.maxsize
    while i < iter_max and ea > es:
        xi = xip1
        xip1 = g(f, xi)
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
print(newton_raphson_solver(fx, 0.05, 0.05, 10, 5))
