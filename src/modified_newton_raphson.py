import sys
from sympy import *
import math


def g1(f, xi, m):
    d = diff(f, x)
    if d == 0:
        return "false"
    return xi - m * (f.evalf(subs={'x': xi}) / float(d.evalf(subs={'x': xi})))


def g2(f, xi):
    d1 = diff(f, x)
    d2 = diff(d1, x)
    a = f.evalf(subs={'x': xi})
    b = d1.evalf(subs={'x': xi})
    c = d2.evalf(subs={'x': xi})
    if float((b**2)-a*c) == 0:
        return "false"
    return xi - (a*b)/float((b**2)-a*c)


def modified1_newton_raphson_solver(f, x0, es, iter_max, n_sig, m):
    xip1 = x0
    i = 0
    ea = sys.maxsize
    while i < iter_max and ea > es:
        xi = xip1
        xip1 = g1(f, xi, m)
        if xip1 == "false":
            return "Not Solvable by Modified 1 Newton Raphson", ea, i
        if not xip1 == 0:
            xip1 = round(xip1, -int(math.floor(math.log10(abs(xip1)))) + (n_sig - 1))
        print(i, xip1)
        if not (xip1 == 0):
            ea = round(abs(((xip1 - xi) / xip1) * 100), 5)
            print("ea of i =", i, "is equal ", ea)
        if abs(f.evalf(subs={'x': xip1})) < 1e-9:
            if abs(xip1 - int(xip1)) < 10**(-n_sig):
                xip1 = int(xip1)
            return xip1, ea, i
        i += 1
    if i == iter_max and ea > es:
        return "Solution is not found upon the given tolerance.", ea, i
    print("ea final = ", ea)
    if abs(xip1 - int(xip1)) < 10**(-n_sig):
        xip1 = int(xip1)
    return xip1, ea, i


def modified2_newton_raphson_solver(f, x0, es, iter_max, n_sig):
    xip1 = x0
    i = 0
    ea = sys.maxsize
    while i < iter_max and ea > es:
        xi = xip1
        xip1 = g2(f, xi)
        if xip1 == "false":
            return "Not Solvable by Modified 1 Newton Raphson", ea, i
        if not xip1 == 0:
            xip1 = round(xip1, -int(math.floor(math.log10(abs(xip1)))) + (n_sig - 1))
        print(i, xip1)
        if not (xip1 == 0):
            ea = round(abs(((xip1 - xi) / xip1) * 100), 5)
            print("ea of i =", i, "is equal ", ea)
        if abs(f.evalf(subs={'x': xip1})) < 1e-9:
            if abs(xip1 - int(xip1)) < 10**(-n_sig):
                xip1 = int(xip1)
            return xip1, ea, i
        i += 1
    if i == iter_max and ea > es:
        return "Solution is not found upon the given tolerance.", ea, i
    print("ea final = ", ea)
    if abs(xip1 - int(xip1)) < 10 ** (-n_sig):
        xip1 = int(xip1)
    return xip1, ea, i


x = symbols('x')
# fx = x**5 - 11*x**4 + 46*x**3 - 90*x**2 + 81*x - 27
# print("Modified 1")
# print(modified1_newton_raphson_solver(fx, 1.3, 0.00001, 50, 5, 2))
# print("\nModified 2")
# print(modified2_newton_raphson_solver(fx, 0, 0.00001, 50, 5))
