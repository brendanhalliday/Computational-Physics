# This code was modification of Mark Newman's
# squarewell.py. The modification was written by Brendan
# Halliday for PHY407 lab 7 question 2.
# It calculates the first few energy levels of a 
# quantum harmonic and anharmonic oscillator.

import numpy as np
import matplotlib.pyplot as plt

# Constants
m = 9.1094e-31     # Mass of electron
hbar = 1.0546e-34  # Planck's constant over 2*pi
e = 1.6022e-19     # Electron charge
N = 1000
a = 1.0e-11 # in meters
h = a/N # spatial step
V_0 = 50 * e # energy in J

k = 2*V_0/(a**2) 
alpha = np.sqrt(m*k)/hbar

PI = np.pi
X = np.arange(-10*a,10*a,h)


# Potential function
def V(x, n):
    """
    Modified potential for exercise 8.14
    This function is a potential 
    for a quantum oscillator where n is the
    degree on x(should be even).
    """
    return (V_0/a**2) * (x**n)

def f(r,x,E,n):
    psi = r[0]
    phi = r[1]
    fpsi = phi
    fphi = (2*m/hbar**2)*(V(x,n)-E)*psi
    return np.array([fpsi,fphi],float)

def solve(E,n):
    """
    This calculates the wavefunction
    for a particualr energy level
    where n is the degree of the oscillator
    and (-ba, ba) is the interval
    """
    psi = 0.0
    phi = 1.0
    r = np.array([psi,phi],float)

    for x in np.arange(-10*a,10*a,h):
        k1 = h*f(r,x,E,n)
        k2 = h*f(r+0.5*k1,x+0.5*h,E,n)
        k3 = h*f(r+0.5*k2,x+0.5*h,E,n)
        k4 = h*f(r+k3,x+h,E,n)
        r += (k1+2*k2+2*k3+k4)/6

    return r[0]

def solve_mod(E,n):
    """
    Identical to the function above.
    The only difference is that we record 
    each value in psi
    """
    psi = 0.0
    phi = 1.0
    r = np.array([psi,phi],float)
    r_list = []
    #r_list.append(psi)

    for x in np.arange(-10*a,10*a,h):
        k1 = h*f(r,x,E,n)
        k2 = h*f(r+0.5*k1,x+0.5*h,E,n)
        k3 = h*f(r+0.5*k2,x+0.5*h,E,n)
        k4 = h*f(r+k3,x+h,E,n)
        r += (k1+2*k2+2*k3+k4)/6
        r_list.append(r[0])

    return np.array(r_list)

def trapezoidal(func):
    """
    This functions uses the trapezoidal rule
    to calculate the numerical value of a definite 
    intergral.
    """
    s_t = 0.5 * func[0] + 0.5 * func[-1]

    for k in range(1, len(X)):
        s_t += func[k]
    trap = h * s_t
    return trap

def theory_0(X):
    """
    Theoretical ground state of QHO
    """
    theory = []
    for x in X:
        psi = np.exp(-(alpha/2)*(x**2))*(alpha / PI)**(1/4)
        theory.append(psi)
    return theory

def theory_1(x):
    """
    Theoretical first excited state of QHO
    """
    theory = []
    for x in X:
        psi = (np.sqrt(2*alpha) * x) * np.exp(-(alpha/2) * (x**2))*(alpha / PI)**(1/4)
        theory.append(psi)
    return theory

def theory_2(x):
    """
    Theoretical second excited state of QHO
    """
    theory = []
    for x in X:
        psi = ((2*alpha*(x**2) - 1)/np.sqrt(2)) * np.exp(-(alpha/2) * (x**2))*((alpha / PI)**(1/4))
        theory.append(psi)
    return theory

if __name__ == "__main__":

    E1, E2 = [100, 400, 650],[150, 450, 700] # initialize a list of guesses for E1 and E2
    # initialize lists to hold numerical values for energy once calcualted.
    Harmonic = []
    Anharmonic = []
    for n in [2,4]:

        for i in range(len(E1)):
            
            E_1 = E1[i] * e
            E_2 = E2[i] * e
            psi2 = solve(E_1, n)

            target = e/1000
            while abs(E_1-E_2)>target:
                psi1, psi2 = psi2, solve(E_2, n)
                E_1, E_2 = E_2, E_2-psi2*(E_2-E_1)/(psi2-psi1)

            print("E =",E_2/e,"eV")

            # now append the numerical values found for both
            # potentials to their corresponding lists
            if n == 2:
                Harmonic.append(E_2)
            else:
                Anharmonic.append(E_2)

    # now we want to calculate the wavefuction numerically from
    # using the solve() function.
    # additionally, we want to normalize each each have function
    # by calculating its integral numerically using trapezoidal rule
    
    # This section calculates, normalizes and plots the 
    # anharmonic wavefunctions
    for i in range(len(Anharmonic)):
        psi = solve_mod(Anharmonic[i],4) # calculate unnormalized psi
        print(len(psi))
        A = trapezoidal(psi**2) # calculate normalization constant
        psi = psi/np.sqrt(A) # normalize psi
        plt.plot(X, psi, label="Anharmonic level: " + str(i))
    plt.xlabel("$x$")
    plt.ylabel("$\Psi$")
    plt.title('Numerical values for anharmonic QO')
    plt.grid()
    plt.tight_layout()
    plt.legend()
    plt.savefig('Wavefunctions0.png', dpi=300)
    plt.show()
    

    # This section calculates, normalizes and plots the 
    # harmonic wavefunctions 
    for i in range(len(Harmonic)):
        psi = solve_mod(Harmonic[i],2) # calculate unnormalized psi
        A = trapezoidal(abs(psi)**2) # calculate normalization constant
        psi = psi/np.sqrt(A) # normalize psi
        plt.plot(X, psi, label="Harmonic level: " + str(i))

    # this just plots the theoretical curves
    plt.plot(X, theory_0(X), label="Theoretical ground state")
    plt.plot(X, theory_1(X), label="Theoretical first excited state")
    plt.plot(X, theory_2(X), label="Theoretical second excited state")
    
    # below is just some plotting information
    plt.xlabel("$x$")
    plt.ylabel("$\Psi$")
    plt.title('Numerical and theoetical values for a harmonic QO')
    plt.grid()
    plt.tight_layout()
    plt.legend()
    plt.savefig('Wavefunctions1.png', dpi=300)
    plt.show()
    

            

            



