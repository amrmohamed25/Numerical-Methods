import tkinter

import matplotlib.pyplot
import numpy
import time
from sympy import sympify, lambdify, diff
from tkinter import *
from tkinter import ttk, filedialog
from matplotlib.pyplot import subplots
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)


def disable(x):
    x["state"] = DISABLED


def enable(x):
    x["state"] = NORMAL


def setTable():
    treev.place(x=250, y=100, height=220, width=1000)
    verscrlbar.place(x=1255, y=100, width=20, height=220)
    horscrlbar.place(x=250, y=330, height=20, width=1000)


def clearTable():
    treev.delete(*treev.get_children())
    treev.place(x=0, y=0, width=0, height=0)
    verscrlbar.place(x=0, y=0, width=0, height=0)
    horscrlbar.place(x=0, y=0, width=0, height=0)


def resetCanvasAndToolbar():
    if 'toolbarFrame' in globals():
        toolbarFrame.grid_remove()
        canvas.get_tk_widget().place_forget()
        matplotlib.pyplot.cla()


def plot(equ_str, prec):
    if not equ_str or not prec or not lower.get() or not higher.get():
        return
    window.geometry("1390x800")
    resetCanvasAndToolbar()
    clearTable()
    if equ_str != equ.getvar("equ"):
        iters.setvar(value=1)
    equ.setvar(name="equ", value=equ_str)
    expr = sympify(equ_str)
    f = lambdify('x', expr)
    # list of squares
    y = [f(i) for i in numpy.linspace(float(lower.get()), float(higher.get()), 100).tolist()]

    # adding the subplot
    l, u, mid, ea = bisection2(float(lower.get()), float(higher.get()), f, float(prec), int(iters.getvar()))
    plot1.axvline(x=l, color='green')
    plot1.axvline(x=u, color='purple')
    plot1.axvline(x=mid, color='red')
    # plotting the graph
    plot1.plot([i for i in numpy.linspace(float(lower.get()), float(higher.get()), 100).tolist()], y)

    # creating the Tkinter canvas
    # containing the Matplotlib figure
    global canvas
    canvas = FigureCanvasTkAgg(fig, master=window)
    canvas.draw()

    # placing the canvas on the Tkinter window
    canvas.get_tk_widget().place(x=250, y=100)

    global toolbarFrame
    toolbarFrame = tkinter.Frame(master=window)
    toolbarFrame.grid(row=9, column=0, rowspan=3)
    # creating the Matplotlib toolbar
    toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)
    toolbar.update()
    if ea is not None and (ea <= float(prec) or int(iters.get()) == int(iters.getvar())):
        disable(plot_button)
    iters.setvar(value=int(iters.getvar()) + 1)


def bisection(l, u, f, prec, iter, old_mid=None):
    mid = (l + u) / 2
    fmid = f(mid)
    ea = None if old_mid is None else abs((mid - old_mid) / mid)
    if fmid == 0:
        ea = "Exact"
    treev.insert("", 'end', values=((iters.getvar()), l, f(l), u, f(u), mid, f(mid), ea))
    iters.setvar(value=int(iters.getvar()) + 1)
    if iter == 1 or fmid == 0:
        return mid, ea
    if old_mid is not None and ea <= prec:
        return mid, ea
    else:
        if f(u) * fmid > 0:
            return bisection(l, mid, f, prec, iter - 1, mid)
        else:
            return bisection(mid, u, f, prec, iter - 1, mid)


def bisection2(l, u, f, prec, iter, old_mid=None):
    mid = (l + u) / 2
    fmid = f(mid)
    ea = None if old_mid is None else abs((mid - old_mid) / mid)
    if fmid == 0:
        return l, u, mid, ea
    if iter == 1:
        return l, u, mid, ea
    if old_mid is not None and ea <= prec:
        return l, u, mid, ea
    else:
        if f(u) * fmid > 0:
            return bisection2(l, mid, f, prec, iter - 1, mid)
        else:
            return bisection2(mid, u, f, prec, iter - 1, mid)


def falseposition(l, u, f, prec, iter, old_mid=None):
    mid = (l * f(u) - u * f(l)) / (f(u) - f(l))
    fmid = f(mid)
    ea = None if old_mid is None else abs((mid - old_mid) / mid)
    if fmid == 0:
        ea = "Exact"
    treev.insert("", 'end', values=(iters.getvar(), l, f(l), u, f(u), mid, f(mid), ea))
    iters.setvar(value=int(iters.getvar()) + 1)
    if iter == 1 or fmid == 0:
        return mid, ea
    if old_mid is not None and ea <= prec:
        return mid, ea
    else:
        if fmid > 0:
            return falseposition(l, mid, f, prec, iter - 1, mid)
        else:
            return falseposition(mid, u, f, prec, iter - 1, mid)


def fixedPoint(xr, g, prec, iter):
    xr_old = xr
    xr = g(xr_old)
    ea = abs((xr - xr_old) / xr)
    if g(xr) == 0:
        ea = "Exact"
    treev.insert("", 'end', values=(iters.getvar(), xr_old, g(xr_old), "", "", xr, g(xr), ea))
    iters.setvar(value=int(iters.getvar()) + 1)
    if ea > prec:
        if iter == 1 or g(xr) == 0:
            return xr, ea
        iter = iter - 1
        return fixedPoint(xr, g, prec, iter)
    else:
        return xr, ea


def newton_raphson(x_i_1, f, f_prime, prec, iter):
    x_i = x_i_1 - f(x_i_1) / f_prime(x_i_1)
    ea = abs((x_i - x_i_1) / x_i)
    if f(x_i) == 0:
        ea = "Exact"
    treev.insert("", 'end', values=(iters.getvar(), x_i_1, f(x_i_1), "", "", x_i, f(x_i), ea))
    iters.setvar(value=int(iters.getvar()) + 1)
    if iter == 1 or (str(ea) != "Exact" and ea <= prec) or f(x_i) == 0:
        return x_i, ea
    else:
        return newton_raphson(x_i, f, f_prime, prec, iter - 1)


def secant_method(x_i_1, x_i, f, prec, iter):
    x_new = x_i - f(x_i) * (x_i - x_i_1) / (f(x_i) - f(x_i_1))
    ea = abs((x_new - x_i) / x_new)
    if f(x_new) == 0:
        ea = "Exact"
    treev.insert("", 'end',
                 values=(iters.getvar(), x_i_1, f(x_i_1), x_i, f(x_i), x_new, f(x_new), ea))
    iters.setvar(value=int(iters.getvar()) + 1)
    if iter == 1 or (str(ea) != "Exact" and ea <= prec) or f(x_new) == 0:
        return x_new, ea
    else:
        return secant_method(x_i, x_new, f, prec, iter - 1)


def submit(scanner, equ_str, prec, iter):
    window.geometry("1390x800")
    iters.setvar(value=1)
    try:
        try:
            clearTable()
            if scanner == 0 or scanner == 1 or scanner == 4:
                if float(lower.get()) >= float(higher.get()):
                    root_label.config(text="Lower should be smaller than higher")
                    root_label.grid(row=2, column=2, columnspan=10)
                    return root_str
        except:
            pass
        setTable()
        if equ_str != equ.getvar("equ") and scanner == 0:
            enable(plot_button)
        treev["columns"] = ("0", "1", "2", "3", "4", "5", "6", "7")

        # Defining heading
        treev['show'] = 'headings'

        # Assigning the width and anchor to  the
        # respective columns
        treev.column("0", minwidth=50, width=60, anchor='c', stretch=False)
        treev.column("1", minwidth=50, width=180, anchor='c', stretch=False)
        treev.column("2", minwidth=50, width=180, anchor='c', stretch=False)
        treev.column("3", minwidth=50, width=180, anchor='c', stretch=False)
        treev.column("4", minwidth=50, width=180, anchor='c', stretch=False)
        treev.column("5", minwidth=50, width=180, anchor='c', stretch=False)
        treev.column("6", minwidth=50, width=180, anchor='c', stretch=False)
        treev.column("7", minwidth=50, width=180, anchor='c', stretch=False)

        # Assigning the heading names to the
        # respective columns
        treev.heading("0", text="i")
        treev.heading("1", text="lower")
        treev.heading("2", text="f(lower)")
        treev.heading("3", text="upper")
        treev.heading("4", text="f(upper)")
        treev.heading("5", text="r")
        treev.heading("6", text="f(r)")
        treev.heading("7", text="ea")
        if len(prec) == 0:
            prec = 0.00001
        if len(iter) == 0:
            iter = 50
        expr = sympify(equ_str)
        f = lambdify('x', expr)
        isFixedDiverging = 0;
        prec = float(prec)
        iter = int(iter)
        fl = f(float(lower.get()))
        start = time.time()

        if scanner == 0 or scanner == 1:
            fu = f(float(higher.get()))
            if fl * fu < 0:
                if scanner == 0:
                    enable(plot_button)
                    root, ea = bisection(float(lower.get()), float(higher.get()), f, prec, iter)
                elif scanner == 1:
                    root, ea = falseposition(float(lower.get()), float(higher.get()), f, prec, iter)
            elif fl == 0:
                root = float(lower.get())
                ea = None
            elif fu == 0:
                root = float(higher.get())
                ea = None
            else:
                root = None
                ea = None
        elif scanner == 2:
            f_prime = lambdify('x', diff(expr))
            if abs(f_prime(float(lower.get()))) < 1:
                root, ea = fixedPoint(float(lower.get()), lambdify('x', expr), prec, iter)
            else:
                root = None
                isFixedDiverging = 1
        elif scanner == 3:
            root, ea = newton_raphson(float(lower.get()), lambdify('x', sympify(equ_str)),
                                      lambdify('x', diff(sympify(equ_str))),
                                      prec, iter)
        elif scanner == 4:
            root, ea = secant_method(float(lower.get()), float(higher.get()), lambdify('x', expr), prec, iter)
        end = time.time()
        if root is None:
            if isFixedDiverging == 1:
                root = "Your function is diverging"
            else:
                root = "Can't find the root!!"
        else:
            root = "Root is " + str(root) + ", Relative error is " + str(ea) + ", Number of iterations are " + str(
                iters.getvar() - 1) + ", Execution Time is " + str(end - start)
        resetCanvasAndToolbar()
        iters.setvar(value=1)
        root_label.config(text=root)
    except:
        root_label.config(text="Error occured.")
    root_label.grid(row=2, column=2, columnspan=10)
    return root_str


def validate_equation(inp):
    try:
        expr = sympify(inp)
        lambdify('x', expr)
    except:
        return False
    return True


def validate_float(inp):
    try:
        float(inp)
    except:
        return False
    return float(inp) >= 0


def validate_int(inp):
    try:
        int(inp)
    except:
        return False
    return int(inp) > 0


def get_combo_box(event):
    scanner = event.widget.current()
    clearTable()
    disable(plot_button)
    if scanner == 0:
        enable(plot_button)
        higher.grid()
        higher_label.grid()
    elif scanner == 1:
        root_label.config(text="")
        higher.grid()
        higher_label.grid()
    elif scanner == 2:
        root_label.config(text="")
        higher.grid_remove()
        higher_label.grid_remove()
    elif scanner == 3:
        root_label.config(text="")
        higher.grid_remove()
        higher_label.grid_remove()

    elif scanner == 4:
        root_label.config(text="")
        higher.grid()
        higher_label.grid()
    else:
        print("Incorrect input")


def read():
    try:
        fileName = filedialog.askopenfile(parent=window, mode='rb', title='Choose a file').read().strip()
    except:
        fileName = ""
    return fileName


if __name__ == '__main__':
    window = Tk()
    fig, plot1 = subplots(figsize=(4, 4))

    color = "black"
    fg_color = "white"
    img = PhotoImage(file="background.png")
    label = Label(
        window,
        image=img
    )
    label.place(x=0, y=0)
    # setting the title
    window.title('Plotting in Tkinter')
    window.geometry("1390x800")
    window.resizable(0, 0)

    tkinter.Label(window, text="Equation:", bg=color, foreground=fg_color).grid(row=0)
    tkinter.Label(window, text="Precision:", bg=color, foreground=fg_color).grid(row=1)
    tkinter.Label(window, text="Maximum Number of iterations:", bg=color, foreground=fg_color).grid(row=2)
    tkinter.Label(window, text="Method:", bg=color, foreground=fg_color).grid(row=3)
    equ = tkinter.Entry(window, validate='focusout', vcmd=(window.register(validate_equation), '%P'))
    prec = tkinter.Entry(window, validate='key', vcmd=(window.register(validate_float), '%P'))
    prec.insert(0, "0.00001")
    iters = tkinter.Entry(window, validate='key', vcmd=(window.register(validate_int), '%P'))
    iters.insert(0, "50")
    equ.grid(row=0, column=1)
    prec.grid(row=1, column=1)
    iters.grid(row=2, column=1)
    iters.setvar(value=1)

    lower_label = tkinter.Label(window, text="Lower:", bg=color, foreground=fg_color)
    lower_label.grid(row=0, column=2)
    lower = tkinter.Entry(window, validate='key', vcmd=(window.register(validate_equation), '%P'))
    lower.grid(row=0, column=3)

    higher_label = tkinter.Label(window, text="Higher:", bg=color, foreground=fg_color)
    higher_label.grid(row=0, column=4)
    higher = tkinter.Entry(window, validate='key', vcmd=(window.register(validate_equation), '%P'))
    higher.grid(row=0, column=5)

    n = tkinter.StringVar()
    methodnames = tkinter.ttk.Combobox(window, width=19, textvariable=n, state="readonly")

    # Adding combobox drop down list
    methodnames['values'] = (' Bisection',
                             ' False Position',
                             ' Fixed Point',
                             ' Newton Raphson',
                             ' Secant'
                             )

    methodnames.grid(column=1, row=3)
    methodnames.current(0)
    methodnames.bind("<<ComboboxSelected>>", get_combo_box)

    plot_button = tkinter.Button(master=window,
                                 command=lambda: plot(equ.get(), prec.get()),
                                 height=2,
                                 width=10,
                                 text="Plot")
    plot_button.grid(row=7, pady=3)

    root_str = ""
    file = None

    equ.setvar(name="equ", value=None)
    Button(master=window,
           command=lambda: (equ.delete(0, len(equ.get())), equ.insert(0, read())),
           height=2,
           width=10,
           text="Browse").grid(row=6, pady=3)

    Button(master=window,
           command=lambda: submit(methodnames.current(), equ.get(), prec.get(), iters.get()),
           height=2,
           width=10,
           text="Submit").grid(row=5)
    root_label = tkinter.Label(window, text=root_str, bg=color, foreground=fg_color)
    treev = ttk.Treeview(window, selectmode='browse')

    # Calling pack method w.r.to treeview

    # Constructing vertical scrollbar
    # with treeview
    verscrlbar = ttk.Scrollbar(window,
                               orient="vertical",
                               command=treev.yview)
    horscrlbar = ttk.Scrollbar(window,
                               orient="horizontal",
                               command=treev.xview)
    # run the gui
    window.mainloop()
