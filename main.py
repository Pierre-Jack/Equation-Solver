from tkinter import ttk
from tkinter import *
import tkinter as tk
import numpy as np


def gauss_elimination():
    print("Hello gaussElimination")


def gauss_jordan():
    print("Hello gaussJordan")


def cholesky():
    import src.LU_Decomposition as lu
    global row, n
    try:
        A = []
        S = []
        for record in entries:
            A.append([])
            for i in range(0, len(record)-1):
                entry = record[i].get()
                if entry:
                    A[-1].append(int(entry))
                else:
                    A[-1].append(0)
            entry = record[-1].get()
            S.append([int(entry)])
    except ValueError:
        new_label = Label(root, text='Invalid. Try again!')
        new_label.grid(row=row, columnspan=n * 3, sticky='W')
        new_label.after(1000, lambda: new_label.destroy())
        return
    print("Hello Cholesky")
    significant_figs = 5
    l,u = lu.cholesky_lu(A,significant_figs)


def doolittle():
    print("Hello doolittle")


def crout():
    print("Hello crout")


def jacobi_iteration():
    print("Hello JacobiIteration")
    print(Initials)


def gauss_seidel():
    print("Hello GaussSeidel")
    print(Initials)


def done():
    if methodType == "Gauss Elimination":
        gauss_elimination()
    elif methodType == "Gauss-Jordan":
        gauss_jordan()
    elif methodType == "LU-Decomposition":
        if LUType == "Cholesky":
            cholesky()
        elif LUType == "Crout":
            crout()
        elif LUType == "Doolittle":
            doolittle()
    elif methodType == "Jacobi Iteration":
        jacobi_iteration()
    elif methodType == "Gauss-Seidel":
        gauss_seidel()


def get_coeff():
    try:
        if methodType == "LU-Decomposition":
            global LUType
            LUType = (combo2.get())
            combo2.grid_forget()
            label_combo2.grid_forget()
            parm_button.grid_forget()
        elif methodType == "Jacobi Iteration" or methodType == "Gauss-Seidel":
            global InitialTextBox, labelInitial, Initials
            Initials = (InitialTextBox.get())
            InitialTextBox.grid_forget()
            labelInitial.grid_forget()
            parm_button.grid_forget()
    except ValueError:
        label = Label(root, text="Enter Parameters")
        label.grid(row=1, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return
    try:
        global n
        n = int(e.get())
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


def get_parm():
    try:
        global n
        n = int(e.get())
    except ValueError:
        label = Label(root, text="You're supposed to enter a number, Try again")
        label.grid(row=1, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return
    if n < 2:
        label = Label(root, text="At least two variables are required!")
        label.grid(row=1, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return
    try:
        global methodType
        methodType = (combo.get())
    except ValueError:
        label = Label(root, text="Choose a method")
        label.grid(row=3, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return
    my_label.grid_forget()
    my_button.grid_forget()
    combo.grid_forget()
    label_combo.grid_forget()
    e.grid_forget()
    try:
        if methodType == "LU-Decomposition":
            global label_combo2, combo2, LUType
            if combo.get() == "LU-Decomposition":
                label_combo2 = tk.Label(root, text='Choose an LU:')
                label_combo2.grid(row=4, column=0, pady=10, sticky='nsew')
                combo2 = ttk.Combobox(
                    state="readonly",
                    values=["Doolittle", "Crout", "Cholesky"]
                )
                combo2.grid(row=5, column=0, pady=10, sticky='nsew')
                global parm_button
                parm_button = Button(root, text='SubmitParm', command=get_coeff)
                parm_button.grid(row=6, column=0)
        elif methodType == "Jacobi Iteration" or methodType == "Gauss-Seidel":
            global InitialTextBox, labelInitial, Initials
            labelInitial = tk.Label(root, text='Choose initials separated by commas:')
            labelInitial.grid(row=4, column=0, pady=10, sticky='nsew')
            InitialTextBox = tk.Entry(root, width=10, borderwidth=5)
            InitialTextBox.grid(row=5, column=0, pady=10, sticky='nsew')
            parm_button = Button(root, text='SubmitParm', command=get_coeff)
            parm_button.grid(row=6, column=0)
        else:
            get_coeff()
    except ValueError:
        label = Label(root, text="Choose an LU method")
        label.grid(row=5, column=0, columnspan=3)
        label.after(1000, lambda: label.destroy())
        return


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
my_label = tk.Label(root, text='How many variables?')
my_label.grid(row=0, column=0, pady=10, sticky='nsew')
e = tk.Entry(root, width=10, borderwidth=5)
e.grid(row=1, column=0, pady=10, sticky='nsew')
label_combo = tk.Label(root, text='Choose a method:')
label_combo.grid(row=2, column=0, pady=10, sticky='nsew')
combo = ttk.Combobox(
    state="readonly",
    values=["Gauss Elimination", "Gauss-Jordan", "LU-Decomposition", "Jacobi Iteration", "Gauss-Seidel"]
)
combo.grid(row=3, column=0, pady=10, sticky='nsew')

my_button = tk.Button(root, text='SubmitFirst', command=get_parm)
my_button.grid(row=6, column=0, pady=10, sticky='nsew')
root.mainloop()
