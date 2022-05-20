# This program was written by Brendan Halliday
# For Lab 4 Question 3 part b.
# This program compares relaxation to over relaxation 
# for the nonlinear equation from exercise 6.11 in
# 'Computational Physics' by Mark Newman

import numpy as np

def f(x,c):
    """
    RHS of non linear equation under question
    """
    return 1 - np.exp(-c*x)

def relax(f, c, threshold, x, dx):
    """
    This function uses relaxation method to
    calculate the solution of x = f(x)
    The parameters are a threshold value
    an initial guess x and a dx which 
    is just an initialization for 
    the difference between succesive approximations.
    """
    x_list = [x]
    while dx > threshold:
        x_list.append(f(x_list[-1], c))
        dx = np.abs(x_list[-1]-x_list[-2])
    return x_list[-1]

def relax_count(f, c, threshold, x, dx):
    """
    This function calculates the number
    of iterations required for the relaxation
    mehtod to converge to a the true answer within
    a given threshold.
    """
    x_list = [x]
    i = 0
    while dx > threshold:
        x_list.append(f(x_list[-1], c))
        dx = np.abs(x_list[-1]-x_list[-2])
        i += 1
    return i

def over_relax_count(f, c, threshold, x, dx, omega):
    """
    This function uses the over relaxation 
    method.
    """
    x_list = [x]
    i = 0
    while dx > threshold:
        x_list.append((1 + omega) * f(x_list[-1], c) - omega * x_list[-1])
        dx = np.abs(x_list[-1]-x_list[-2])
        i += 1
    return i
if __name__ == "__main__":
    x = 1.0  # initial x guess
    dx = 1.0  # initial distance (just needs to be big)
    threshold = 1e-6  # convergence threshold
    c = 2.0 # value for c
    omega = 0.5 # guess for omega


    print("Number of iterations required for relaxation method for c = 2.0: ")
    print(relax_count(f, c, threshold, x, dx))

    print("Number of iterations required for over relaxation method:")
    print(over_relax_count(f, c, threshold, x, dx, omega))