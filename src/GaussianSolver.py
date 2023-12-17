from lib2to3.pygram import Symbols
from sympy import sympify, N
from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
)
import numpy as np
import math

# import simpy as simpy

# eq1 = [1, 2, 3, 4, 5]
# eq2 = [6, 7, 8, 9, 10]
# eq3 = [11, 11, 12, 14, 15]
# eq4 = [16, 16, 16, 19, 20]
#
# a = [eq1, eq2, eq3, eq4]
# b = [eq1, eq2, eq3, eq4]
# m = np.array([eq1, eq2, eq3])


# m = np.matrix([eq1, eq2, eq3])
def rnd(x, sig):
    if isinstance(x, str):
        return x
    if isinstance(x, list):
        result = []
        for i in x:
            result.append(rnd(i, sig))
        return result
    if abs(x) < 1*10**(-sig+1):
        return 0
    return round(x, -int(math.floor(math.log10(abs(x)))) + (sig - 1))

def get(x, idx):
    if isinstance(x, np.ndarray):
        return x[idx]
    if isinstance(x, float) or isinstance(x, int):
        return x
    if isinstance(x, list):
        return x[idx]


def inv(x):
    if isinstance(x, float) or isinstance(x, int):
        return 1 / x
    if isinstance(x, str):
        return "1/(" + x + ")"


def neg(x):
    if isinstance(x, float) or isinstance(x, int):
        return -x
    if isinstance(x, str):
        return "(-" + x + ")"


def scale(x, c):
    if isinstance(c, float) or isinstance(c, int):
        if isinstance(x, np.ndarray):
            return c * x
        if isinstance(x, float) or isinstance(x, int):
            return c * x
        if isinstance(x, str):
            return "(" + str(c) + "*" + x + ")"
        if isinstance(x, list):
            result = []
            for i in x:
                result.append(scale(i, c))
            return result

    if isinstance(c, str):
        if isinstance(x, list):
            result = []
            for i in x:
                result.append(scale(i, c))
            return result
        return "(" + c + "*" + str(x) + ")"


def add(x, y):
    if isinstance(x, np.ndarray) or isinstance(y, np.ndarray):
        return x + y
    if (isinstance(x, float) or isinstance(x, int)) and (isinstance(y, float) or isinstance(y, int)):
        return x + y
    if isinstance(x, str) or isinstance(y, str):
        return "(" + str(x) + "+" + str(y) + ")"
    if isinstance(x, list) or isinstance(x, list):
        result = []
        for i in range(len(x)):
            result.append(add(get(x, i), get(y, i)))
        return result


def swap(a, idx1, idx2):
    if isinstance(a, list):
        a[idx1], a[idx2] = a[idx2], a[idx1]
    elif isinstance(a, np.ndarray):
        a[[idx1, idx2]] = a[[idx2, idx1]]


def maxIdx(a, col):
    max = col
    for i in range(col, len(a)):
        if not isinstance(a[i][col], str):
            if not (a[i][len(a[i]) - 1] == 0 or a[max][len(a[max]) - 1] == 0):
                if (isinstance(a[max][col], str) or a[i][col] / a[i][len(a[i]) - 1] > a[max][col] / a[max][
                    len(a[max]) - 1]) and a[i][col] > 1e-3:
                    max = i
    return max


def gaussian(a, nsignificant):
    n = len(a)
    for i in range(n):
        # for row in range(i, n):
        #     tmp = a[row][len(a[row])-1]
        #     if tmp != 0 and tmp != 1:
        #         a[row] = scale(a[row], 1/tmp)
        swap(a, i, maxIdx(a, i))
        for k in range(i + 1, n):
            tmp1 = a[i][i]
            tmp2 = a[k][i]
            a[k] = add(a[k], scale(scale(a[i], inv(tmp1)), neg(tmp2)))
            a[k] = rnd(a[k], nsignificant)
            a[k] = [0 if idx == i else a[k][idx] for idx in range(len(a[k]))]


def backSub(a, nsignificant):  # applies backward substitution to a matrix in the upper triangular form
    solution = []
    n = len(a)
    m = len(a[0]) - 1  # cols of the coefficient matrix, col m is the augmented
    for i in range(n - 1, -1, -1):
        tmp = a[i][m]
        for j in range(m - 1, i, -1):
            tmp = add(tmp, neg(scale(a[i][j], solution[j-i-1])))
            tmp = rnd(tmp, nsignificant)

        if not isinstance(tmp, str) and abs(tmp) < 1e-10:
            tmp = 0
        solution = [rnd(scale(tmp, inv(a[i][i])), nsignificant)] + solution
    for i in range(len(solution)):
        if isinstance(solution[i], str):
            e = solution[i]
            e = parse_expr(e)
            solution[i] = N(e, 5)

    return solution


def gaussJordan(a, nsignificant):
    n = len(a)
    solution = []
    for i in range(n):
        swap(a, i, maxIdx(a, i))
        for k in range(n):
            if k == i:
                continue
            tmp1 = a[i][i]
            tmp2 = a[k][i]
            a[k] = add(a[k], scale(scale(a[i], inv(tmp1)), neg(tmp2)))
            a[k] = rnd(a[k], nsignificant)
            a[k] = [0 if idx == i else a[k][idx] for idx in range(len(a[k]))]
    for i in range(n):
        a[i] = rnd(scale(a[i], inv(a[i][i])), nsignificant)
        solution.append(a[i][-1])
    for i in range(len(solution)):
        if isinstance(solution[i], str):
            e = solution[i]
            e = parse_expr(e)
            solution[i] = N(e, 5)
    return solution


# m = m*1.0
# gaussJordan(a)

# print(a)
# print(np.array(backSub(a)))
# print(swap(a, 0, 2))
# p = Symbols('p')
# print(sympify(e))
# print(e)
# print(parse_expr(e))
# print(simpy.simplify(backSub(a)[0]))
# print(np.array(a))
# print(backSub(a))
