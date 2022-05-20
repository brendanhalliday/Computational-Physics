# This program was written by Brendan Halliday and Nikolaos Rizos
# for PHY407 lab 4 question 2
# This program answers exercise 6.9 from 'Computational Physics'
# by Mark Newman
# This program calculates the entries of the
# Hamiltonian of an asymetric well and outputs a matrix of sive m x n.
# Additionally, this program calculates the eigen energies by finding 
# the eigenvalues of the hamiltonian matrix

## SEE LINE 148 for code specific to part e

import numpy as np
from scipy import constants
from scipy.constants import e, hbar, m_e #elementary charge, reduced plank constant, mass of electron
import matplotlib.pyplot as plt
import gaussxw as gaussxw

# PART B
# define helper functions

def Hmn(m, n):
    """
    This function evaluates Hmn (the Hamiltonian matrix for
    an asymmetric well) for an arbitrary m and n and returns
    its value in joules. Here, m and n are the row and column 
    indices respectively.
    """
    # store special constants
    M = m_e
    a = 10.0*1.6022*10**(-19) # energy in J
    L = 5.0e-10 # width of the well in meters
    PI = constants.pi # pi


    # series of conditionals as defined in the lab manual
    if m != n and (m % 2) == (n % 2):
        hmn = 0

    elif m != n and (m % 2) != (n % 2):
        hmn = -1.0 * (8.0 * a * m * n) / ((PI**2) * (m**2 - n**2)**2) 

    elif m == n:
        hmn = 0.5 * a + ((PI**2 * m**2)/2.0)*(hbar**2/(M * L**2))

    else:
        pass

    return hmn 

def Hmatrix(Hmn, mmax, nmax):
    """
    This function creates a m x n
    matrix of the Hamiltonian H
    where the entries are in electron volts.
    """
    # initialize an mmax x nmax matrix of zeros
    H = np.zeros((mmax, nmax))
    
    for m in range(1, mmax + 1): # for each row
        for n in range(1, nmax + 1): # for each column
            H[m-1, n-1] = Hmn(m, n) # replace zero mth row, nth column entry with Hmn
    return H / e #division by electron charge for units of electron volts

def wave_function(x, psi_n):
    """
    This function outputs the value of an eigen function
    for each x. This function takes in an eigenvector of
    the Hamiltonian and, from the components(Fourier Coefficients),
    outputs the corresponding eigen function as a superposition of 
    sin waves.
    """
    L = 5.0 # Angstroms
    psi = psi_n[0] * np.sin(np.pi * x / L)
    for i in range(1, len(psi_n)):
        psi = psi + psi_n[i] * np.sin(np.pi * (i+1) * x / L)
    return psi

def Gauss(f, N, a, b, psi_eigenvector):
    """
    This function returns the approximate
    value of the integral of function f
    using Gaussian quadratures for N steps
    on the interval from a to b. The code
    written in this function is from the file
    gaussint.py. The code was simply copy and pasted 
    into a function for ease of use. 

    Note: This function has been altered to calculate 
    modulus squared of our wavefunction.
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
        # Here we square f() for units of prob density
        s += wp[k]*(f(xp[k],psi_eigenvector)**2)
    return s 


if __name__ == "__main__":

# PART C
    m = 10
    n = m
    H = Hmatrix(Hmn, m, n) # calculate a 10 x 10 Hamiltonian
    print("Calculate eigen energies using a 10x10 Hamiltonian: ")
    EIGENVALUES = np.linalg.eigvalsh(H) # calculate the eigenvalues of this Hamiltonian
    # print eigen energies for n = m = 10
    print("Ground state energy: ")
    print(EIGENVALUES[0], " eV") 
    for i in range(1, len(EIGENVALUES)):
        print("Excited state energy number " + str(i+1) + ":")
        print(EIGENVALUES[i], " eV")

# PART D 
    m = 100
    n = m
    H_100 = Hmatrix(Hmn, m, n) # calculate a 100 x 100 Hamiltonian
    print("Now calculate eigen energies using a 100x100 Hamiltonian: ")
    EIGENVALUES100 = np.linalg.eigvalsh(H_100)[0:10] # calculate eigenvalues
    # print eigen energies for m = n = 100
    print("Ground state energy: ")
    print(EIGENVALUES100[0], " eV") 
    for i in range(1, len(EIGENVALUES)):
        print("Excited state energy number " + str(i+1) + ":")
        print(EIGENVALUES100[i], " eV")

    #Now let's calculate the relative error 
    DIFF = abs(EIGENVALUES - EIGENVALUES100)
    ERROR = []
    for i in range(len(EIGENVALUES100)):
        ERROR.append(DIFF[i]/EIGENVALUES100[i])
    print("Relative error of the 10 x 10 Hamiltonian eigenvales:")
    print("Ground state energy: ")
    print(ERROR[0]) 
    for i in range(1, len(ERROR)):
        print("Excited state energy number " + str(i+1) + ":")
        print(ERROR[i])



# PART E ##########################################################################
    a = 0.0 # Angstroms
    L = 5.0 # Angstroms
    X = np.linspace(0, L, 100) # x axis array
    N = 11 # number numerical integration steps for wavefunction normalization

    psi_eigenvectors = np.linalg.eigh(H_100)[1] #calculate the eigenvector
    
    # begin plotting probability densities
    plt.figure(1)
    for i in range(0,3):
        
        psi_eigenvector = psi_eigenvectors[i] #calculate the eigenvector
        PSI = wave_function(X, psi_eigenvector) # calculate the corresponfind wave function for said eigenvector
        PSI_2 = PSI**2 # square each value of the wave function for units of probability density
        A = Gauss(wave_function, N, a, L,psi_eigenvector) # calculate the integral of the probability density |psi|^2 
        PSI_2 = PSI_2/A # normalize |psi|^2

        plt.plot(X, PSI_2, label = r'$|\psi_{i}(x)|^{2}$' + ' i = ' + str(i + 1))

    plt.xlabel(r'x ($\AA$)')
    plt.ylabel(r'$|\psi(x)|^{2}$')
    plt.title('Probability Densities of Ground state and first two excited states')
    plt.grid()
    plt.xlim(0,L)
    plt.legend()
    plt.savefig('Lab4_Q2.png', dpi=300)
    plt.show()