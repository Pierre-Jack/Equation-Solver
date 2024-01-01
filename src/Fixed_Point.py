import math


def Fixed_Point(x0, eps, numOfIterations, g_function, significant_figures):
    x_new = x0
    k = 1
    # g_function = lambda x: eval(g_function)
    while (k<=numOfIterations) :
        x_old = x_new
        x_new = g_function(x_old)
        x_new  = round(x_new, -int(math.floor(math.log10(abs(x_new)))) + (significant_figures - 1))
        print(x_new)
        if x_new == 0:
            x_new = 1e-10
        abs_error  = abs((x_new-x_old)/x_new)*100

        if abs_error<=eps :
            print(k)
            print(x_new)
            return x_new
        k=k+1
    if k == (numOfIterations + 1):
        print(k)
        print('Cannot Converge.')
        return 'Cannot Converge.'


x0 = 4
eps = 0.0001
numOfIterations = 50
significant_figures = 5
g_function_str = input("Enter the g(x) function as a Python expression, using 'x' as the variable: ")
g_function = g_function = lambda x: eval(g_function_str, {'x': x, 'e': math.e})
x = Fixed_Point(x0, eps, numOfIterations, g_function, significant_figures)



# eq1: (2*x +3)**0.5
# equ2: e**-x