# This program was written by Brendan Halliday
# for Lab 3 Q1 b. The program calculates the
# Diffracted intensity over undiffracted intensity of
# sound waves as described in exercise 5.11 from
# 'Computational Physics' by Mark Newman. 
# Additionally, the program output's the relative difference
# between my function and the one calculated by
# scipy.special.fresnel


import numpy as np
from Lab03_Q1_a_ii import Gauss, relerr
from scipy.special import fresnel
import matplotlib.pyplot as plt

def I_I0(C, S):
    """
    This function calculates the value
    for I/I_0 given a value of u.
    """
    return (1/8) * ((2*C + 1)**2 + (2*S + 1)**2)

def c(t):
    return np.cos(0.5 * np.pi * t**2)

def s(t):
    return np.sin(0.5 * np.pi * t**2)

def C(c, u, N):
    """
    This function returns the C
    integral of the fresnel function
    for particular value u using 
    Gaussian quadratures for N points
    """
    a = 0.0
    return Gauss(c, N, a, u)

def S(s, u, N):
    """
    This function returns the S
    integral of the fresnel function
    for a particular value u using
    Gaussian quadratures for N points
    """
    a = 0.0

    return Gauss(s, N, a, u)

if __name__ == "__main__":

    
    # create array for x
    x = np.linspace(-5, 5, 1000)
    # transform x into u
    wavelength = 1
    z = 3
    u = x * np.sqrt(2 / (wavelength * z))
    # call Scipy fresnel functions
    Ssci, Csci = fresnel(u)
    # set N = 50 for Gaussian Quadratures
    N = 50
    # this section initializes an empty list 
    # I and then calculates I/I_0 for u 
    # and appends I/I_0(u) for each u
    IG = []
    for i in u:
        IG.append(I_I0(C(c, i , N), S(s, i, N)))

    Isci = I_I0(Csci, Ssci)
    
    # this block plots I/I_0 using my function and Gaussian
    # quadratures and compares it to the one
    # obtained with scipy.special.fresnel
    plt.figure(1)
    plt.plot(x, IG, color = 'r', label='My function $I_{G}$')
    plt.plot(x, Isci, color = 'b', linestyle= '--', label= 'Scipy $I_{SP}$')
    plt.xlabel('metres')
    plt.ylabel(r'$\frac{I}{I_{0}}$')
    plt.title('Diffracted intensity over undiffracted intensity of sound waves')
    plt.legend()
    plt.grid()
    plt.xlim(-5, 5)
    plt.savefig('Diffracted_sound_waves.png', dpi=300)
    plt.show()

    # the following block of code plots the relative 
    # error for each
    plt.figure(2)
    plt.plot(x, relerr(Isci, IG), color = 'r')
    plt.xlabel('metres')
    plt.ylabel(r'$\delta(x)$')
    plt.title('Relative difference between $I_{SP}$ and $I_{G}$ for N = 50')
    plt.grid()
    plt.xlim(-5, 5)
    plt.savefig('Relative_difference_waves.png', dpi=300)
    plt.show()


    # Now we calculate the maximum relative
    # difference for for each N from 3 to 50

    # create an array N values
    N = np.arange(3, 51, 1)
    # initialize empty list MAX
    MAX = []

    # create new array for x from 0 to 5
    x = np.linspace(0, 5, 500)
    #tranform into u
    u = x * np.sqrt(2 / (wavelength * z))
    #evaluate new fresnel functions using scipy
    Ssci, Csci = fresnel(u)
    for k in N:
        # this section initializes an empty list 
        # I and then calculates I/I_0 for u 
        # and appends I/I_0(u) for each u
        IG = []
        for i in u:
            IG.append(I_I0(C(c, i , k), S(s, i, k)))
        
        # find maximum value in I and append said
        # value to MAX to be graphed for each N
        MAX.append(max(IG))

    #this block plots values for the maximum relative difference
    plt.figure(3)
    plt.plot(N, MAX, 'rx')
    plt.xlabel('N')
    plt.ylabel(r'Max $\delta (x)$')
    plt.title(' Max relative difference between $I_{SP}$ and $I_{G}$ for varying N for 0 < x < 5')
    plt.grid()
    plt.savefig('Relative_difference_max.png', dpi=300)
    plt.show()
