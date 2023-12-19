import numpy as np
import math
from numpy import linalg as LA

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

def doolittle_lu(a, sig_figs):
    n = np.shape(a)[0]

    l = np.zeros((n, n))
    u = np.zeros((n, n))

    for j in range(n):
        l[j][j] = 1

        for i in range(j + 1):
            temp = sum(u[k][j] * l[i][k] for k in range(i))
            u[i][j] = a[i][j] - temp

            u[i][j] = rnd(u[i][j], sig_figs)

        for i in range(j, n):
            temp = sum(u[k][j] * l[i][k] for k in range(j))
            l[i][j] = (a[i][j] - temp) / u[j][j]

            l[i][i] = rnd(l[i][j], sig_figs)
    return l, u


def crout_lu(a, sig_figs):
    n = np.shape(a)[0]

    l = np.zeros((n, n))
    u = np.eye(n)



    for i in range(n):
        l[i][0] = a[i][0]
    for j in range(1, n):
        u[0][j] = a[0][j] / l[0][0]

        u[0][j] = rnd(u[0][j], sig_figs)

    for j in range(1, n - 1):
        for i in range(j, n):
            temp = sum(l[i][k] * u[k][j] for k in range(j))
            l[i][j] = a[i][j] - temp

            l[i][j] = rnd(l[i][j], sig_figs)

        for k in range(j + 1, n):
            temp = sum(l[j][i] * u[i][k] for i in range(j))
            u[j][k] = (a[j][k] - temp) / l[j][j]
            u[j][k] = rnd(u[j][k], sig_figs)

    l[n - 1][n - 1] = a[n - 1][n - 1] - sum(l[n - 1][k] * u[k][n - 1] for k in range(n))
    l[n - 1][n - 1] = rnd(l[n - 1][n - 1], sig_figs)

    return l, u


def check_if_valid_for_cholesky(a):
    if not np.array_equal(a, a.transpose()):
        return False
    elif not a.shape[0] == a.shape[1]:
        return False
    else:
        for x in LA.eig(a)[0]:
            if x < 0:
                return False
    return True


def cholesky_lu(a, sig_figs):
    n = np.shape(a)[0]

    l = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1):
            temp = 0
            if (j == i):
                for k in range(j):
                    temp += pow(l[j][k], 2)
                    temp = rnd(temp, sig_figs)
                l[j][j] = int(math.sqrt(a[j][j] - temp))
                l[j][j] = rnd(l[j][j], sig_figs)
            else:
                for k in range(j):
                    temp += (l[i][k] * l[j][k])
                    temp = rnd(temp, sig_figs)

                if (l[j][j] > 0):
                    l[i][j] = int((a[i][j] - temp) / l[j][j])
                    l[i][j] = rnd(l[i][j], sig_figs)
    u = l.transpose()

    return l, u


def solve_lu(l, u, b, sig_figs):
    y = np.matmul(LA.inv(l), b)
    for i in range(y.shape[0]):
        y[i][0] = rnd(y[i][0], sig_figs)
    x = np.matmul(LA.inv(u), y)
    for i in range(x.shape[0]):
        x[i][0] = rnd(x[i][0], sig_figs)
    return x

