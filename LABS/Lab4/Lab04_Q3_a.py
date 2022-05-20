# This program was written by Brendan Halliday
# For Lab 4 Question 3 part a.
# This program computes the solution to the nonlinear
# equation and outputs the value of x.

import numpy as np
import matplotlib.pyplot as plt

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

if __name__ == "__main__":
    x = 1.0  # initial x guess
    dx = 1.0  # initial distance (just needs to be big)
    threshold = 1e-6  # convergence threshold


    c_0 = 2.0 # value for c
    print("This is the approximation for c = 2.0")
    print(relax(f, c_0, threshold, x, dx))

    c = np.arange(0.0, 3.01, 0.01) # values for c
    X = [] # initialize empty list

    # calculate the value of x for each value of c
    # and append siad value to X
    for i in c:  
        X.append(relax(f, i, threshold, x, dx)) 

    # This block just plots X as a function of c
    plt.figure(1)
    plt.plot(c, X,'b-')
    plt.xlabel('Values of c')
    plt.ylabel('Values of x')
    plt.title(r'Relaxation method for $1 - e^{-cx} = x$')
    plt.grid()
    plt.xlim(0,3)
    plt.savefig('Lab4_Q3a.png', dpi=300)
    plt.show()