import sys
import math
max_iter = 100
sig = 5
es = 0.00001

def rnd(x):
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


def f(x):
    y= x**2 - 2
    return rnd(y)

def secant(Xi0, Xi1):
    ea = sys.maxsize
    i = 0
    print("i\tXi-1\t Xi \tXi+1\tEa")
    while i < max_iter-1 and ea > es:
        Xi2 = Xi1 - f(Xi1)*(Xi1 - Xi0) / (f(Xi1) - f(Xi0))
        Xi1 = rnd(Xi1)
        Xi2 = rnd(Xi2)
        ea = abs((Xi2-Xi1)/Xi2)
        # ea = rnd(ea)
        print(i,"\t",Xi0,"\t",Xi1,"\t",Xi2,"\t",ea)
        Xi0 = Xi1
        Xi1 = Xi2
        i += 1
    if ea <= es:
        print("Root found at X=", Xi1, "with tolerance", ea)
        print("Number of iterations needed: ", i)
        # print("Runtime: ")
    elif ea > es and i >= max_iter:
        print("Maximum Iterations reached(", max_iter, "). No root found under given tolerance")
        
secant(0.5, 1.0)

