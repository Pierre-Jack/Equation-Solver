from tkinter import ttk
from tkinter import *
import tkinter as tk
import numpy as np
from numpy import linalg as LA
import itertools
from timeit import default_timer as timer
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
NavigationToolbar2Tk)
from lib2to3.pygram import Symbols
from sympy import sympify, N
from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
)
from time import sleep
from numpy import *
import math
import sys
from sympy import *

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

def secant(Xi0, Xi1, es, max_iter):
    ea = sys.maxsize
    i = 0
    print("i\tXi-1\t Xi \tXi+1\tEa")
    while i <= max_iter-1 and ea > es:
        Xi2 = Xi1 - rnd(f(Xi1))*(Xi1 - Xi0) / (rnd(f(Xi1)) - rnd(f(Xi0)))
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
        return Xi1, i
        # print("Runtime: ")
    elif ea > es and i >= max_iter:
        return "Maximum Iterations reached. No root found under given tolerance", max_iter


def bsct(u, l, tol, maxI, sig):
    if u < l:
        tmp=u
        u=l
        l=tmp
    prev = 0
    r = (l+u)/2
    count = 0
    while(abs(r-prev) > tol or count==0):
        if f(r)*f(l) < 0:
            l, u = l, r
        else:
            l, u = r, u
        prev = r
        r = (l+u)/2
        r = round(r, -int(math.floor(math.log10(abs(r)))) + (sig - 1))
        count = count+1
        if count > maxI:
            return "max iterations exceeded", count
    print("num of iterations = " + str(count))
    return r, count

def rf(u, l, tol, maxI, sig):
    if u < l:
        tmp=u
        u=l
        l=tmp
    prev = 0
    r = (l*f(u) - u*f(l))/(f(u)-f(l))
    count = 0
    while(abs(r-prev) > tol or count==0):
        if f(r)*f(l) < 0:
            l, u = l, r
        else:
            l, u = r, u
        prev = r
        r = (l*f(u) - u*f(l))/(f(u)-f(l))
        r = round(r, -int(math.floor(math.log10(abs(r)))) + (sig - 1))
        count = count+1
        if count > maxI:
            return "max iterations exceeded", count
    print("num of iterations = " + str(count))
    return r, count

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
            return x_new, k
        k=k+1
    if k == (numOfIterations + 1):
        print(k)
        print('Cannot Converge.')
        return 'Cannot Converge.', k


def g1(f, xi, m):
    x = symbols('x')
    d = diff(f, x)
    if d == 0:
        return "false"
    return xi - m * (f.evalf(subs={'x': xi}) / float(d.evalf(subs={'x': xi})))


def g2(f, xi):
    x = symbols('x')
    d1 = diff(f, x)
    d2 = diff(d1, x)
    a = f.evalf(subs={'x': xi})
    b = d1.evalf(subs={'x': xi})
    c = d2.evalf(subs={'x': xi})
    if float((b**2)-a*c) == 0:
        return "false"
    return xi - (a*b)/float((b**2)-a*c)


def modified1_newton_raphson_solver(f, x0, es, iter_max, n_sig, m):
    x = symbols('x')
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
        i += 1
    if i == iter_max and ea > es:
        return "Solution is not found upon the given tolerance.", ea, i
    print("ea final = ", ea)
    return xip1, ea, i


def modified2_newton_raphson_solver(f, x0, es, iter_max, n_sig):
    x = symbols('x')
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
        i += 1
    if i == iter_max and ea > es:
        return "Solution is not found upon the given tolerance.", ea, i
    print("ea final = ", ea)
    if xip1 - int(xip1) < 10E-9:
        xip1 = int(xip1)
    return xip1, ea, i


# x = symbols('x')
# fx = x ** 3 - 5 * (x**2) + 7 * x - 3
# print("Modified 1")
# print(modified1_newton_raphson_solver(fx, 0, 0.00001, 50, 5, 2))
# print("\nModified 2")
# print(modified2_newton_raphson_solver(fx, 0, 0.00001, 50, 5))

# x0 = 4
# eps = 0.0001
# numOfIterations = 50
# significant_figures = 5
# g_function_str = input("Enter the g(x) function as a Python expression, using 'x' as the variable: ")
# g_function = g_function = lambda x: eval(g_function_str, {'x': x, 'e': math.e})
# x = Fixed_Point(x0, eps, numOfIterations, g_function, significant_figures)


def get_solution():
    l = float(guesses[0])
    try:
        u = float(guesses[1])
    except:
        pass

    if methodType == "Newton-Raphson(Original)" or methodType == "Newton-Raphson(Modified-1)" or methodType == "Newton-Raphson(Modified-2)":
        code = """def f(x):
                        return """ + func_f.replace("e", "(math.e)")
        exec(code, globals())
    else:
        code = """def f(x):
                        return """ + func_f.replace("e", "(math.e)").replace("sin", "math.sin").replace("cos","math.cos")
        exec(code, globals())

    if methodType == "Newton-Raphson(Original)":
        x = symbols('x')
        fx = f(x)
        s,ea,i = modified1_newton_raphson_solver(fx, l, eps, maxI, sig, 1)
        return s, i

    if methodType == "Newton-Raphson(Modified-1)":
        x = symbols('x')
        fx = f(x)
        s,ea,i = modified1_newton_raphson_solver(fx, l, eps, maxI, sig, u)
        return s, i


    if methodType == "Newton-Raphson(Modified-2)":
        x = symbols('x')
        fx = f(x)
        s,ea,i = modified2_newton_raphson_solver(fx, l, eps, maxI, sig)
        return s, i

    if methodType == "Bisection":
        return bsct(u, l, eps, maxI, sig)
    if methodType == "False-Position":
        return rf(u, l, eps, maxI, sig)
    if methodType == "Fixed point":
        g_function_str = func_g
        g_function = lambda x: eval(g_function_str, {'x': x, 'e': math.e, 'sin': math.sin, 'cos': math.cos})
        return Fixed_Point(l, eps, maxI, g_function, sig)
    if methodType == "Secant Method":
        return secant(l, u, eps, maxI)

    return 1.5


def solve():
    global func_f
    func_f = input_function.get()
    print("f(x)=", func_f)

    global func_g
    if methodType == "Fixed point":
        func_g = g.get()
        print("g(x)=", func_g)


    try:
        global guesses
        guesses = InitialTextBox.get().split(",")
        print("guesses=", guesses)

    except ValueError:
        label = Label(root, text="error")
        label.grid(row=3, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return

    label_input_function.grid_forget()
    input_function.grid_forget()
    labelInitial.grid_forget()
    InitialTextBox.grid_forget()
    plot_button.grid_forget()
    parm_button.grid_forget()

    if methodType == "Fixed point":
        g.grid_forget()
        label_g.grid_forget()
    global solution

    start_time = timer()
    solution, iterations = get_solution()
    end_time = timer()
    time_taken = round(((end_time - start_time) * 1000000), 5)
    print(time_taken)

    new_label = Label(root, text="root" + ': ' + str(solution))
    new_label.grid(row=1, columnspan=2, sticky='W')

    new_label = Label(root, text="time taken: " + str(time_taken) + " µs.")
    new_label.grid(row=2, columnspan=2, sticky='W')

    new_label = Label(root, text="iterations reached: " + str(iterations))
    new_label.grid(row=3, columnspan=2, sticky='W')



def plot():
    global func
    func = input_function.get()
    code = """def f(x):
                    return """ + func.replace("e", "(math.e)")
    exec(code, globals())
    global func_p
    func_p = func
    code = """def p(x):
                        return """ + func_p.replace("e", "(math.e)").replace("sin", "np.sin").replace("cos", "np.cos")
    exec(code, globals())

    fig = Figure(figsize=(4, 4),dpi=100)

    x = np.arange(-10, 10, 0.1)
    y = p(x)

    plot1 = fig.add_subplot(111)

    plot1.plot(x, y)

    plot1.axhline(y=0, color='k')

    canvas = FigureCanvasTkAgg(fig, master=root)

    canvas.draw()

    canvas.get_tk_widget().grid(row=9, column=0, columnspan=3)

    toolbar = NavigationToolbar2Tk(canvas,root)

    toolbar.update()

    canvas.get_tk_widget().grid(row=10, column=0, columnspan=3)



def get_parm():
    try:
        global eps
        eps = e.get()
        if eps == '':
            eps = 0.00001
        eps = float(eps)
    except ValueError:
        label = Label(root, text="You're supposed to enter a number, Try again")
        label.grid(row=1, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return

    try:
        global sig
        sig = sig_fig.get()
        if sig == '':
            sig = 10
        sig = int(sig)
    except ValueError:
        label = Label(root, text="You're supposed to enter an integer, Try again")
        label.grid(row=3, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return

    try:
        global maxI
        maxI = N_Iterations.get()
        if maxI == '':
            maxI = 50
        maxI = int(maxI)
    except ValueError:
        label = Label(root, text="You're supposed to enter an integer, Try again")
        label.grid(row=5, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return

    try:
        global methodType
        methodType = (combo.get())
        print(type(methodType))
    except ValueError:
        label = Label(root, text="Choose a method")
        label.grid(row=7, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return
    if methodType == '':
        label = Label(root, text="Choose a method")
        label.grid(row=7, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return

    my_label.grid_forget()
    my_button.grid_forget()
    combo.grid_forget()
    label_combo.grid_forget()
    e.grid_forget()
    sig_fig.grid_forget()
    label_sig_fig.grid_forget()
    N_Iterations.grid_forget()
    N_IterationsLabel.grid_forget()
    print(eps, sig, maxI, methodType)

    global label_input_function
    label_input_function = tk.Label(root, text='Enter f(x)')
    label_input_function.grid(row=1, column=0, pady=10, sticky='nsew')
    global input_function
    input_function = tk.Entry(root, width=10, borderwidth=5)
    input_function.grid(row=2, column=0, pady=10, sticky='nsew')
    global plot_button
    plot_button = Button(master=root,
                         command=plot,
                         height=2,
                         width=10,
                         text="Plot")
    global parm_button
    plot_button.grid(row=8, column=0, pady=5, sticky='nsew')
########################################################################
    global InitialTextBox, labelInitial, Initials
    labelInitial = tk.Label(root, text='Choose initial guess(es) separated by a comma:')
    labelInitial.grid(row=3, column=0, pady=5, sticky='nsew')
    InitialTextBox = tk.Entry(root, width=10, borderwidth=5)
    InitialTextBox.grid(row=4, column=0, pady=5, sticky='nsew')


    if methodType == "Fixed point":
        global label_g
        label_g = tk.Label(root, text='Enter g(x)')
        label_g.grid(row=5, column=0, pady=10, sticky='nsew')
        global g
        g = tk.Entry(root, width=10, borderwidth=5)
        g.grid(row=6, column=0, pady=10, sticky='nsew')



    parm_button = Button(root, text='Next', command=solve)
    parm_button.grid(row=12, column=0)




########################################################################




root = tk.Tk()
root.resizable(width=False, height=False)
root.title('Numerical Project')
root.minsize(250, 200)
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()
x_position = (screen_width - 400) // 2
y_position = (screen_height - 300) // 2
root.geometry(f'+{x_position}+{y_position}')
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)
root.rowconfigure(1, weight=1)
root.rowconfigure(2, weight=1)
root.rowconfigure(3, weight=1)
root.rowconfigure(4, weight=1)
my_label = tk.Label(root, text='EPS')
my_label.grid(row=0, column=0, pady=10, sticky='nsew')
e = tk.Entry(root, width=10, borderwidth=5)
e.grid(row=1, column=0, pady=10, sticky='nsew')

global label_sig_fig
label_sig_fig = tk.Label(root, text='Significant Figures')
label_sig_fig.grid(row=2, column=0, pady=10, sticky='nsew')
global sig_fig
sig_fig = tk.Entry(root, width=10, borderwidth=5)
sig_fig.grid(row=3, column=0, pady=10, sticky='nsew')


global N_Iterations, N_IterationsLabel
N_IterationsLabel = tk.Label(root, text='Maximum number of Iterations:')
N_IterationsLabel.grid(row=4, column=0, pady=10, sticky='nsew')

N_Iterations = tk.Entry(root, width=10, borderwidth=5)
N_Iterations.grid(row=5, column=0, pady=5, sticky='nsew')

label_combo = tk.Label(root, text='Choose a method:')
label_combo.grid(row=6, column=0, pady=10, sticky='nsew')
combo = ttk.Combobox(
    state="readonly",
    values=["Bisection", "False-Position", "Fixed point", "Newton-Raphson(Original)", "Newton-Raphson(Modified-1)", "Newton-Raphson(Modified-2)", "Secant Method"]
)
combo.grid(row=7, column=0, pady=10, sticky='nsew')

my_button = tk.Button(root, text='Next', command=get_parm)
my_button.grid(row=9, column=0, pady=10, sticky='nsew')
root.mainloop()
