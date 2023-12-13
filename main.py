from tkinter import ttk
from tkinter import *
import tkinter as tk


def done():
    print("hello")


def submit():
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
label_combo = tk.Label(root, text='Choose a language:')
label_combo.grid(row=2, column=0, pady=10, sticky='nsew')
combo = ttk.Combobox(
    state="readonly",
    values=["Gauss Elimination", "Gauss-Jordan", "LU-Decomposition", "Jacobi Iteration", "Gauss-Seidel"]
)
combo.grid(row=3, column=0, pady=10, sticky='nsew')
my_button = tk.Button(root, text='Submit', command=submit)
my_button.grid(row=4, column=0, pady=10, sticky='nsew')
root.mainloop()
