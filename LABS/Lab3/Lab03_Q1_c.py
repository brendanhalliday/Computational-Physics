# This program was written by Brendan Halliday
# This program creates 2D contour plots for I/I_0 as a function
# of x and z position for vaious values of lambda

import numpy as np
import gaussxw
from scipy.special import fresnel
import matplotlib.pyplot as plt

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


def I_I0(C, S):
    """
    This function calculates the value
    for I/I_0 for given C(u) and S(u).
    """
    return (1/8) * ((2*C + 1)**2 + (2*S + 1)**2)

def c(t):
    """
    Integrand of C(u)
    """
    return np.cos(0.5 * np.pi * t**2)

def s(t):
    """
    Integrand of S(u)
    """
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

    # create list with different values of wavelength
    wavelength = [1,2,4]
    # create domain for x and z
    x = np.linspace(-3, 10, 100)
    z = np.linspace(1, 5, 100)
    X, Z = np.meshgrid(x, z)

    # iterate this block for each value of lambda
    for wavelength in wavelength:
        # transform X and Z into U
        U = X * np.sqrt(2 / (wavelength * Z))
        

        N = 50
        # transform each value from 2D array U
        # into I/I_0 and then replace said value with corresponding 
        # entry in U
        for i in range(len(z)):
            for k in range(len(x)):
                I = I_I0(C(c, U[i][k], N), S(s, U[i][k], N))
                U[i][k] = I
    
        #this block plotes the 2d contour plotes for each lambda
        plt.figure(1)
        cp = plt.contourf(X, Z, U)
        cbar = plt.colorbar(cp)
        cbar.ax.set_ylabel('Dimensionless Fraction $I/I_{0}$', rotation=270)
        cbar.ax.get_yaxis().labelpad = 15
        plt.xlabel('x meters')
        plt.ylabel('z meters')
        plt.title('Contour plot of $I/I_{0}$ as a function of the x and z position for $\lambda =$ ' + str(wavelength))
        plt.grid()
        
        plt.savefig('Diffracted_sound_waves_' + str(wavelength) +'.png', dpi=300)
        plt.show()