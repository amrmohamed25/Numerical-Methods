# Numerical-Methods

Objective:

The aim of this assignment is to compare and analyze the behavior of numerical
methods studied in class {Bisection, False-position, Fixed point, Newton-Raphson,
Secant}.

Description:

You are required to implement a root finder program which takes as an input the
equation, the technique to use and its required parameters (e.g. interval for the bisection
method).

Specification:

The program must contain the following features:

● An interactive GUI that enables the user to enter equations containing different
functions such as: {poly, exp, cos, sin}. Reading from files must be available as
well.

● Differentiation and Parsing is your task.

● A way to choose a method to solve the given equation.

● A way to enter the precision and the max number of iterations otherwise default
values are used,
Default Max Iterations = 50, Default Epsilon = 0.00001;

● The answer for the chosen method indicates the number of iterations, execution
time, all iterations, approximate root, and precision.

● Compute the theoretical bound of the error for the methods.
