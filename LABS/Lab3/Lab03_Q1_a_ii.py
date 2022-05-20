# This program was written by Brendan Halliday
# This program calculates the value of the integral in question 1 of
# lab 3. It approximates this integral using Trapezoidal,
# Simpson's and Gaussian quadratures.
# Lastly, it outputs the relative errors for each method
# in a plot.
# It also estimates the error for Gaussian quadratures 
# and graphs those points in the same plot

import numpy as np
import gaussxw as gaussxw
import matplotlib.pyplot as plt

# define helper functions for more readable 
# python code

def f(x):
    """
    Define the function under question
    """
    return 4 / (1 + x**2)

def trapezoidal(f, N, a, b):
    """
    This function returns the approximate
    value of the integral of function f
    using the trapezoidal rule for N steps
    on the interval from a to b.
    """
    h = (b - a) / N
    s = 0.5 * f(a) + 0.5 * f(b)
    for k in range(1, N):
        s += f(a + k * h)
    
    return h * s

def simpson(f, N, a, b):
    """
    This function returns the approximate
    value of the integral of function f
    using Simpson's rule for N steps
    on the interval from a to b.
    """

    h = (b - a) / N
    s = f(a) + f(b)

    #for odd terms
    for k in range(1, N, 2):
        s += 4 * f(a + k * h)
    #for even terms
    for k in range(2, N, 2):
        s += 2 * f(a + k * h)

    return (1 / 3) * h * s
 
def Gauss(f, N, a, b):
    """
    This function returns the approximate
    value of the integral of function f
    using Gaussian quadratures for N steps
    on the interval from a to b. The code
    written in this function is from the file
    gaussint.py. The code was simply copy and pasted 
    into a function for ease of use.
    """
    # Calculate the sample points and weights, then map them 
    # to the required integration domain 
    
    # x and w are arrays. if b is an array, we cannot proceed
    x,w = gaussxw.gaussxw(N)
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w

    # Perform the integration 
    s = 0.0 
    #iterate through each point k
    for k in range(N): 
        s += wp[k]*f(xp[k]) 
    return s 

def relerr(I, I_0):
    """
    This function calculates the relative
    error/difference between the expected I
    and the measured I_0
    """
    return abs(I - I_0)/ I

def G_esterr(Gauss, f, N, a, b):
    """
    This function calculates the estimated error
    of Gaussian quadratures as presented in equation
    1 of the lab manual
    """
    I_2N = Gauss(f, 2*N, a, b)
    I_N = Gauss(f, N, a, b)
    return abs(I_2N - I_N)


if __name__ == "__main__":

    #this is simply to store values of N for future use
    N = []
    for i in range(3, 12):
        N.append(2**i)
    #define bounds of integration
    a = 0.0
    b = 1.0

    #initialize empty lists to store vlaues for error
    trap_error = []
    simp_error = []
    Gauss_error_est = []
    Gauss_error = []
    #define exact value of pi
    PI = np.pi
    #iterate for each value of N
    for i in N:
        print('Integration for N = ', i)
        print()
        #calculate and print approximation using trapezoidal rule
        # and print the relative error
        I_0 = trapezoidal(f, i, a, b)
        err = relerr(PI, I_0)
        # stor error in corresponding emty list
        trap_error.append(err)
        print('Trapezoidal rule: ')
        print('I = ', I_0)
        print('relative error in I = ', err)
        #calculate and print approximation using Simpson's rule
        # and print the relative error
        I_0 = simpson(f, i, a, b)
        err = relerr(PI, I_0)
        #store values of error in corresponding list
        simp_error.append(err)
        print("Simpson's rule: ")
        print('I = ', I_0)
        print('relative error in I =', err)
        
        #calculate and print approximation using Gaussian quadratures
        # and print the relative and estimated error
        I_0 = Gauss(f, i, a, b)
        est_err = G_esterr(Gauss, f, i, a, b)
        err = relerr(PI, I_0)
        #store values of error in corresponding list
        Gauss_error_est.append(est_err)
        Gauss_error.append(err)
        print('Gaussian Quadratures: ')
        print('I = ', I_0)
        print('estimated error =', est_err)
        print('relative error = ', err)
        
        print()

    # the following block of code plots the relative 
    # error for each method
    plt.figure(1)
    plt.plot(N, trap_error ,'rx', label='Trapezoidal')
    plt.plot(N, simp_error,'bx', label="Simpson's")
    plt.plot(N, Gauss_error, 'gx', label = 'Gaussian Quadratures')
    plt.plot(N, Gauss_error_est, 'yx', label = 'Gaussian Quadratures estimated error')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('N')
    plt.ylabel('Relative error')
    plt.title('Relative error as a function of N slices')
    plt.legend()
    plt.grid()
    plt.savefig('relative_error.png', dpi=300)
    plt.show()