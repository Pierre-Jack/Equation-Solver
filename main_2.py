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



def bsct(u, l, tol, maxI):
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
        count = count+1
        if count > maxI:
            return "max iterations exceeded"
    print("num of iterations = " + str(count))
    return r

def rf(u, l, tol, maxI):
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
        count = count+1
        if count > maxI:
            return "max iterations exceeded"
    print("num of iterations = " + str(count))
    return r



def get_solution():
    l = float(guesses[0])
    u = float(guesses[1])

    code = """def f(x):
                    return """ + func_f
    exec(code, globals())

    if methodType == "Bisection":
        return bsct(u, l, eps, maxI)
    if methodType == "False-Position":
        return rf(u, l, eps, maxI)

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
    solution = get_solution()
    end_time = timer()
    time_taken = round(((end_time - start_time) * 1000000), 5)
    print(time_taken)

    new_label = Label(root, text="root" + ' = ' + str(solution))
    new_label.grid(row=1, columnspan=2, sticky='W')

    new_label = Label(root, text="time taken: " + str(time_taken) + " µs.")
    new_label.grid(row=2, columnspan=2, sticky='W')



def plot():
    global func
    func = input_function.get()
    code = """def f(x):
                    return """ + func
    exec(code, globals())


    fig = Figure(figsize=(4, 4),dpi=100)

    x = np.arange(-100, 100, 0.1)
    y = f(x)

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
            sig = 5
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
