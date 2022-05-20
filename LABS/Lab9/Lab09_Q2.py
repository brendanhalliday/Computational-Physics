# This programwas written by Brendan Halliday
# for PHY407 Lab 9 Question 2.
# This program 

import numpy as np
import scipy.constants as constant


attosecond = 1.0e-18 # in seconds
M = constant.electron_mass # in kg
hbar = constant.hbar # in eV(attosecond)
eV = constant.electron_volt
L = 1.0e-8 # in meters

SIGMA = L / 25.0 # constant in units of L
KAPPA = 500.0/L # constant in units of L^-1
P = 1024 # number of segments in x
X = np.linspace(-L/2, L/2, P + 1) # interval of x in units of L


TAU = attosecond # time increment in attoseconds
OMEGA = 3.0e15 # in units of rad/s

V_0 = 6.0e-17/eV # V_0 in electronVolts
X_1 = L/4 # in units of L

DELTA = L/P # in units of L

def initial_psi(x, x_0):
    """
    Unormalized initial condition of psi at t = 0
    """
    a = -((x - x_0)**2)/(4*SIGMA**2)
    b = 1j*KAPPA*x
    upstairs = a + b
    return np.exp(upstairs)

def squarewell(x):
    """
    Return potential for infinite square well
    defined on the domain of x. x should be a
    numpy array.
    """
    return 0.0*x 

def QHO(x):
    """
    Return the potential for a quantum harmonic
    oscillator on the domain x. x should be a numpy
    array.
    """
    return (1/2)*M*(OMEGA**2)*(x**2)

def doublewell(x):
    """
    Return the potential for a doublewell
    on the domain x. x should be a numpy array.
    """
    return V_0*(x**2/X_1**2 - 1)**2

def A(delta):
    """
    Value of matirx element A of the discretized 
    Hamiltonian for some delta 
    """
    return (-hbar**2)/(2*M*(delta**2))

def B_p(p, delta, V):
    """
    Value of matrix element B_p of discretized
    Hamiltonian for some delta, p, and potential V
    """
    return V(p*delta - L/2) - 2*A(delta)

def HD(P, V, A_value):
    """
    Returns a (P-1)x(P-1) array of the discretized 
    Hamiltonian. A should be computed before using
    A(delta) function. expressed in units of eV
    """
    n = P - 1
    HD = []
    # This will be the first row
    hd = np.zeros(n)
    hd[1] = A_value
    hd[0] = B_p(1, DELTA, V)
    HD.append(hd)
    
    for i in range(n - 2):
        # generate the next row so long as that row
        # is not the bottom row
        hd = np.zeros(n)
        hd[i] = A_value
        hd[i + 1] = B_p(i + 1, DELTA, V)
        hd[i + 2] = A_value
        HD.append(hd)
    # generate the last row (similar to the first row)
    hd = np.zeros(n)
    hd[-1] = B_p(P - 1, DELTA, V)
    hd[-2] = A_value
    HD.append(hd)
    return np.array(HD) # return the array 

def LHS(H):
    """
    Computes the LHS operator for the 
    Crank-Nicolson method.
    """

    I = np.identity(P - 1) # generate a P-1 x P-1 identity matrix
    return I + ((TAU*1j)/(2*hbar))*H

def RHS(H):
    """
    Computes the RHS operator for the 
    Crank-Nicolson method
    """

    I = np.identity(P - 1) # generate a P-1 x P-1 identity matrix
    return I - ((TAU*1j)/(2*hbar))*H

def norm_constant(psi):
    """
    Return the normalization wavefunction
    """
    norm = np.trapz(abs(psi)**2, dx=DELTA) # calculate the normalization constant
    norm_const = 1/np.sqrt(norm)
    return norm_const

def average_x(psi):
    """
    Return the value of expectation value of
    position from the wavefunction.
    """
    psi_dagger = np.conjugate(psi)
    return np.trapz(np.multiply(psi_dagger, np.multiply(X[1:-1],psi)), DELTA)

def energy(H, psi):
    """
    Return the value for energy of the wavefunction 
    """
    psi_dagger = np.conjugate(psi)
    return np.trapz(np.multiply(psi_dagger, H.dot(psi)), dx=DELTA)
    
if __name__ == "__main__":

    # For part (a)
    X_0 = L/5 # in units of L
    N = 3000 # number of time steps
    ENERGY = []
    NORM = []
    AVERAGE = []
    
    # The following lines of code calculates the square integral of the initial wavefunction. 
    # Then, the inverse of the square root of the value is found to give the normalization constant
    INITIAL_PSI = initial_psi(X, X_0) # calculate initial wavefunction 
    INITIAL_PSI[0] = 0.0 # ends must be zero
    INITIAL_PSI[-1] = 0.0 # ends must be zero

    norm = norm_constant(INITIAL_PSI)
    NORM.append(norm)
    INITIAL_PSI = norm * INITIAL_PSI
    # initialize array to store spatial and time data for each point of wavefunction
    psi = np.zeros([N+1, P+1], dtype = "complex_")
    psi[0,:] = INITIAL_PSI # set first row as intitial psi
    A_value = A(DELTA)
    # compute the value of the discretized Hamiltonian once
    # for the square potential well
    HD_MAT = HD(P, squarewell, A_value)
    

    #make sure to calculate the energy after normalization
    E = energy(HD_MAT, INITIAL_PSI[1:-1]) # calculate energy value for initial wavefunction
    ENERGY.append(E) # append this energy to ENERGY list


    #A = average_x(INITIAL_PSI[1:-1])
    #AVERAGE.append(A)

    # Now compute the RHS and LHS matrices for the Crank-Nicolson method
    RHS_MAT = RHS(HD_MAT) 
    LHS_MAT = LHS(HD_MAT)

    # now iterate for each value in time to generate the evolving wave equation
    for n in range(N):

        # calcualte the vector v on the right hand side
        v = RHS_MAT.dot(psi[n,1:-1]) 
        updated_psi = np.linalg.solve(LHS_MAT, v) # solve linear equations

        # additionally, we need to caluclate the energy and ensure that the 
        # wavefunction is normalized throughout
        norm = norm_constant(updated_psi)
        updated_psi = norm * updated_psi

        NORM.append(norm)

        # make sure to calculate the energy of the normalized wavefunction
        E = energy(HD_MAT, updated_psi)
        ENERGY.append(E)
        
        # Lastly, update psi
        psi[n + 1, 1:-1] = updated_psi

    
    ENERGY = np.array(ENERGY)/ENERGY[0]
    NORM = np.array(NORM)
    np.savez('plotting_data', ENERGY=ENERGY, NORM=NORM,psi=psi)

    #part c
    N = 4000

    psi = np.zeros([N+1, P+1], dtype = "complex_")
    psi[0,:] = INITIAL_PSI # set first row as intitial psi
    HD_MAT = HD(P, QHO, A_value)
    RHS_MAT = RHS(HD_MAT) 
    LHS_MAT = LHS(HD_MAT)

    for n in range(N):

        # calcualte the vector v on the right hand side
        v = RHS_MAT.dot(psi[n,1:-1]) 
        updated_psi = np.linalg.solve(LHS_MAT, v) # solve linear equations

        # wavefunction is normalized throughout
        norm = norm_constant(updated_psi)
        updated_psi = norm * updated_psi
        
        # Lastly, update psi
        psi[n + 1, 1:-1] = updated_psi

    np.savez('QHO_data', psi=psi)
