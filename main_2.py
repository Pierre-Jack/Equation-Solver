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


def gauss_elimination(A, S):
    for i in range(len(A)):
        A[i].append(S[i][0])
    print(A)
    import src.GaussianSolver as Gs
    global nSignificant, row
    print("Hello gaussElimination")
    row = 0
    for widget in root.winfo_children():
        widget.destroy()
    start_time = timer()
    Gs.gaussian(A, int(nSignificant))
    x = Gs.backSub(A, int(nSignificant))
    end_time = timer()
    time_taken = round(((end_time - start_time) * 1000000), 5)
    print(time_taken)
    for i in range(n):
        new_label = Label(root, text=chr(97 + i) + ' = ' + str(x[i]))
        new_label.grid(row=row, columnspan=2, sticky='W')
        row += 1
    new_label = Label(root, text="time taken: " + str(time_taken) + " µs.")
    new_label.grid(row=row, columnspan=2, sticky='W')
    print(x)

def gauss_jordan(A, S):
    for i in range(len(A)):
        A[i].append(S[i][0])
    print(A)
    import src.GaussianSolver as Gs
    global nSignificant, row
    print("Hello gaussElimination")
    row = 0
    for widget in root.winfo_children():
        widget.destroy()
    start_time = timer()
    x = Gs.gaussJordan(A, int(nSignificant))
    end_time = timer()
    time_taken = round(((end_time - start_time) * 1000000), 5)
    for i in range(n):
        new_label = Label(root, text=chr(97 + i) + ' = ' + str(x[i]))
        new_label.grid(row=row, columnspan=2, sticky='W')
        row += 1
    new_label = Label(root, text="time taken: " + str(time_taken) + " µs.")
    new_label.grid(row=row, columnspan=2, sticky='W')
    print(x)


def cholesky(A, S):
    global nSignificant, row
    import src.LU_Decomposition as lu
    row = 0
    for widget in root.winfo_children():
        widget.destroy()
    if lu.check_if_valid_for_cholesky(A):
        start_time = timer()
        l, u = lu.cholesky_lu(A, int(nSignificant))
        x = lu.solve_lu(l, u, S, int(nSignificant))
        end_time = timer()
        time_taken = round(((end_time - start_time) * 1000000), 5)
    else:
        x = "A is not symmetric positive definite, therefore Cholesky's method cannot be applied to it"
        print(x)
        new_label = Label(root, text=x)
        new_label.grid(row=row, columnspan=5, sticky='W')
        return
    for i in range(n):
        new_label = Label(root, text=chr(97 + i) + ' = ' + str(x[i][0]))
        new_label.grid(row=row, columnspan=2, sticky='W')
        row += 1
    new_label = Label(root, text="time taken: " + str(time_taken) + " µs.")
    new_label.grid(row=row, columnspan=2, sticky='W')


def doolittle(A, S):
    global row
    for widget in root.winfo_children():
        widget.destroy()
    if abs(A[0][0]) < 1e-8:
        x = "The system cannot be solved using Doolittle LU decomposition since A[0][0] = 0 or its value is too small"
        print(x)
        row = 0
        new_label = Label(root, text=x)
        new_label.grid(row=row, columnspan=5, sticky='W')
        return
    else:
        global nSignificant
        import src.LU_Decomposition as lu
        start_time = timer()
        l, u = lu.doolittle_lu(A, int(nSignificant))
        x = lu.solve_lu(l, u, S, int(nSignificant))
        end_time = timer()
        time_taken = round(((end_time - start_time) * 1000000), 5)
        row = 0
        for i in range(n):
            new_label = Label(root, text=chr(97 + i) + ' = ' + str(x[i][0]))
            new_label.grid(row=row, columnspan=2, sticky='W')
            row += 1
        new_label = Label(root, text="time taken: " + str(time_taken) + " µs.")
        new_label.grid(row=row, columnspan=2, sticky='W')
        print(x)


def crout(A, S):
    global nSignificant
    import src.LU_Decomposition as lu
    start_time = timer()
    l, u = lu.doolittle_lu(A, int(nSignificant))
    x = lu.solve_lu(l, u, S, int(nSignificant))
    end_time = timer()
    time_taken = round(((end_time - start_time) * 1000000), 5)
    for widget in root.winfo_children():
        widget.destroy()
    global row
    row = 0
    for i in range(n):
        new_label = Label(root, text=chr(97 + i) + ' = ' + str(x[i][0]))
        new_label.grid(row=row, columnspan=2, sticky='W')
        row += 1
    new_label = Label(root, text="time taken: " + str(time_taken) + " µs.")
    new_label.grid(row=row, columnspan=2, sticky='W')
    print(x)


def jacobi_iteration(A, S):
    import src.JacobiSeidel as JM
    global nSignificant, iterations, relative_error
    B = list(itertools.chain(*S))
    initialGuess = list(map(int, Initials.split(",")))
    start_time = timer()
    x = JM.jacobi(A, B, initialGuess, int(iterations), int(relative_error), int(nSignificant))
    end_time = timer()
    time_taken = round(((end_time - start_time) * 1000000), 5)
    for widget in root.winfo_children():
        widget.destroy()
    global row
    row = 0
    if isinstance(x, str):
        x = "Cannot Converge"
        print(x)
        new_label = Label(root, text=x)
        new_label.grid(row=row, columnspan=5, sticky='W')
        return
    for i in range(n):
        new_label = Label(root, text=chr(97 + i) + ' = ' + str(x[i]))
        new_label.grid(row=row, columnspan=2, sticky='W')
        row += 1
    new_label = Label(root, text="time taken: " + str(time_taken) + " µs.")
    new_label.grid(row=row, columnspan=2, sticky='W')
    print(x)


def gauss_seidel(A, S):
    import src.JacobiSeidel as GS
    global nSignificant, iterations, relative_error
    B = list(itertools.chain(*S))
    initialGuess = list(map(int, Initials.split(",")))
    start_time = timer()
    x = GS.jacobi(A, B, initialGuess, int(iterations), int(relative_error), int(nSignificant))
    end_time = timer()
    time_taken = round(((end_time - start_time) * 1000000), 5)
    for widget in root.winfo_children():
        widget.destroy()
    global row
    row = 0
    if isinstance(x, str):
        x = "Cannot Converge"
        print(x)
        new_label = Label(root, text=x)
        new_label.grid(row=row, columnspan=5, sticky='W')
        return
    for i in range(n):
        new_label = Label(root, text=chr(97 + i) + ' = ' + str(x[i]))
        new_label.grid(row=row, columnspan=2, sticky='W')
        row += 1
    new_label = Label(root, text="time taken: " + str(time_taken) + " µs.")
    new_label.grid(row=row, columnspan=2, sticky='W')
    print(x)


def done():
    global row, n, is_str
    is_str = False
    try:
        A = []
        A_str = []
        S = []
        for record in entries:
            A.append([])
            A_str.append([])
            for i in range(0, len(record) - 1):
                entry = record[i].get()
                if entry:
                    if methodType != "Gauss Elimination" and methodType != "Gauss-Jordan":
                        A[-1].append(float(entry))
                    else:
                        if entry.isalpha():
                            A_str[-1].append(entry)
                            is_str = True
                        else:
                            A_str[-1].append(float(entry))
                            A[-1].append(float(entry))
                else:
                    A[-1].append(0)
                    A_str[-1].append(0)
            entry = record[-1].get()
            S.append([float(entry)])
    except ValueError:
        new_label = Label(root, text='Invalid. Try again!')
        new_label.grid(row=row, columnspan=n * 3, sticky='W')
        new_label.after(1000, lambda: new_label.destroy())
        return

    if not is_str and abs(LA.det(np.array(A))) < 1e-8:
        screen_width = root.winfo_screenwidth()
        screen_height = root.winfo_screenheight()
        x_position = (screen_width - 400) // 2
        y_position = (screen_height - 300) // 2
        root.geometry(f'+{x_position}+{y_position}')
        root.columnconfigure(0, weight=1)
        x = "Matrix is singular"
        for widget in root.winfo_children():
            widget.destroy()
        print(x)
        row = 0
        new_label = Label(root, text=x)
        new_label.grid(row=row, column=0, sticky='W')
        return

    if methodType == "Gauss Elimination":
        gauss_elimination(A_str, S)
    elif methodType == "Gauss-Jordan":
        gauss_jordan(A_str, S)
    elif methodType == "LU-Decomposition":
        if LUType == "Cholesky":
            cholesky(np.array(A), np.array(S))
        elif LUType == "Crout":
            crout(np.array(A), np.array(S))
        elif LUType == "Doolittle":
            doolittle(np.array(A), np.array(S))
    elif methodType == "Jacobi Iteration":
        jacobi_iteration(A, S)
    elif methodType == "Gauss-Seidel":
        gauss_seidel(A, S)


def get_coeff():
    global nSignificant
    nSignificant = '5'
    try:
        nSignificant = (sig_fig.get())
        if nSignificant == '':
            nSignificant = '5'
    except ValueError:
        nSignificant = '5'
        return
    try:
        sig_fig.grid_forget()
        label_sig_fig.grid_forget()
        parm_button.grid_forget()
        if methodType == "LU-Decomposition":
            global LUType
            LUType = (combo2.get())
            combo2.grid_forget()
            label_combo2.grid_forget()
            parm_button.grid_forget()
        elif methodType == "Jacobi Iteration" or methodType == "Gauss-Seidel":
            global InitialTextBox, labelInitial, Initials, N_Iterations, N_IterationsLabel, Relative_error, Relative_errorLabel, iterations, relative_error
            Initials = (InitialTextBox.get())
            iterations = int(N_Iterations.get())
            relative_error = float(Relative_error.get())
            InitialTextBox.grid_forget()
            labelInitial.grid_forget()
            N_Iterations.grid_forget()
            N_IterationsLabel.grid_forget()
            Relative_error.grid_forget()
            Relative_errorLabel.grid_forget()
            parm_button.grid_forget()
    except ValueError:
        label = Label(root, text="Invalid Input")
        label.grid(row=1, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return
    try:
        global n
        n = e.get()
        if n == '':
            n = 0.00001
        n = float(n)
    except ValueError:
        label = Label(root, text="You're supposed to enter a number, Try again")
        label.grid(row=1, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return
    e['state'] = DISABLED
    my_label.grid_forget()
    my_button.grid_forget()
    combo.grid_forget()
    label_combo.grid_forget()
    e.grid_forget()
    global row, entries
    row = 0
    entries = []
    for i in range(n):
        Label(root, text='').grid(row=row, columnspan=3 * n + 1)
        row += 1
        entries.append([])
        col = 0
        for j in map(chr, range(97, 97 + n)):
            entry = Entry(root, width=5, borderwidth=2, justify='right')
            entries[i].append(entry)
            entry.grid(row=row, column=col)
            col += 1
            Label(root, text=j).grid(row=row, column=col, padx=5, sticky='W')
            col += 1
            if col == 3 * n - 1:
                Label(root, text='=').grid(row=row, column=col)
                col += 1
                entry = Entry(root, width=5, borderwidth=2, justify='right')
                entries[i].append(entry)
                entry.grid(row=row, column=col)
            else:
                Label(root, text='+').grid(row=row, column=col)
            col += 1
        row += 1
    Label(root, text='').grid(row=row)
    row += 1
    global new_button
    new_button = Button(root, text='Submit', command=done)
    new_button.grid(row=row, column=3 * n)
    row += 1

# def f(x):
#     return np.sin(x)

# func = "x**2-9"
#

def plot():

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

    canvas.get_tk_widget().grid(row=7, column=0, columnspan=3)

    toolbar = NavigationToolbar2Tk(canvas,root)

    toolbar.update()

    canvas.get_tk_widget().grid(row=8, column=0, columnspan=3)



def get_parm():
    try:
        global n
        n = e.get()
        if n == '':
            n = 0.00001
        n = float(n)
    except ValueError:
        label = Label(root, text="You're supposed to enter a number, Try again")
        label.grid(row=1, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return

    try:
        global methodType
        methodType = (combo.get())
        print(type(methodType))
    except ValueError:
        label = Label(root, text="Choose a method")
        label.grid(row=3, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return
    if methodType == '':
        label = Label(root, text="Choose a method")
        label.grid(row=3, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return
    ###########################################################
    def f(x):
        return 0.001* x**3 + 1
    ###########################################################

    my_label.grid_forget()
    my_button.grid_forget()
    combo.grid_forget()
    label_combo.grid_forget()
    e.grid_forget()
    global label_sig_fig
    label_sig_fig = tk.Label(root, text='Enter number of Significant Figures')
    label_sig_fig.grid(row=2, column=0, pady=10, sticky='nsew')
    global sig_fig
    sig_fig = tk.Entry(root, width=10, borderwidth=5)
    sig_fig.grid(row=3, column=0, pady=10, sticky='nsew')

    global label_input_function
    label_input_function = tk.Label(root, text='Enter f(x)')
    label_input_function.grid(row=4, column=0, pady=10, sticky='nsew')
    global input_function
    input_function = tk.Entry(root, width=10, borderwidth=5)
    input_function.grid(row=5, column=0, pady=10, sticky='nsew')

    plot_button = Button(master=root,
                         command=plot,
                         height=2,
                         width=10,
                         text="Plot")
    plot_button.grid(row=6, column=0, pady=10, sticky='nsew')



    # try:
    #     if methodType == "LU-Decomposition":
    #         global label_combo2, combo2, LUType
    #         if combo.get() == "LU-Decomposition":
    #             label_combo2 = tk.Label(root, text='Choose an LU:')
    #             label_combo2.grid(row=4, column=0, pady=10, sticky='nsew')
    #             combo2 = ttk.Combobox(
    #                 state="readonly",
    #                 values=["Doolittle", "Crout", "Cholesky"]
    #             )
    #             combo2.grid(row=5, column=0, pady=10, sticky='nsew')
    #             global parm_button
    #             parm_button = Button(root, text='Next', command=get_coeff)
    #             parm_button.grid(row=6, column=0)
    #     elif methodType == "Jacobi Iteration" or methodType == "Gauss-Seidel":
    #         global InitialTextBox, labelInitial, Initials, N_Iterations, N_IterationsLabel, Relative_error, Relative_errorLabel
    #         labelInitial = tk.Label(root, text='Choose initials separated by commas:')
    #         labelInitial.grid(row=4, column=0, pady=10, sticky='nsew')
    #         InitialTextBox = tk.Entry(root, width=10, borderwidth=5)
    #         InitialTextBox.grid(row=5, column=0, pady=10, sticky='nsew')
    #
    #         N_IterationsLabel = tk.Label(root, text='Maximum number of Iterations:')
    #         N_IterationsLabel.grid(row=6, column=0, pady=10, sticky='nsew')
    #         N_Iterations = tk.Entry(root, width=10, borderwidth=5)
    #         N_Iterations.grid(row=7, column=0, pady=10, sticky='nsew')
    #
    #         Relative_errorLabel = tk.Label(root, text='Relative Error:')
    #         Relative_errorLabel.grid(row=8, column=0, pady=10, sticky='nsew')
    #         Relative_error = tk.Entry(root, width=10, borderwidth=5)
    #         Relative_error.grid(row=9, column=0, pady=10, sticky='nsew')
    #
    #         parm_button = Button(root, text='Next', command=get_coeff)
    #         parm_button.grid(row=10, column=0)
    #     else:
    #         parm_button = Button(root, text='Next', command=get_coeff)
    #         parm_button.grid(row=4, column=0)
    # except ValueError:
    #     label = Label(root, text="Choose an LU method")
    #     label.grid(row=5, column=0, columnspan=3)
    #     label.after(1000, lambda: label.destroy())
    #     return


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
label_combo = tk.Label(root, text='Choose a method:')
label_combo.grid(row=2, column=0, pady=10, sticky='nsew')
combo = ttk.Combobox(
    state="readonly",
    values=["Bisection", "False-Position", "Fixed point", "Newton-Raphson(Original)", "Newton-Raphson(Modified-1)", "Newton-Raphson(Modified-2)", "Secant Method"]
)
combo.grid(row=3, column=0, pady=10, sticky='nsew')

my_button = tk.Button(root, text='Next', command=get_parm)
my_button.grid(row=6, column=0, pady=10, sticky='nsew')
root.mainloop()
