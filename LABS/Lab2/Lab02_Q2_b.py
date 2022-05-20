# This program was written by Brendan Halliday
# For Lab 2 question 2 part b
# The program computes Bessel functions and compares 
# them to scipy.special Bessel functions

import numpy as np
import matplotlib.pyplot as plt
from scipy import special


#steps
N = 1000
#order of Bessel functions we need to calculate
n = [0, 3, 5]
#colours to differentiate the different order Bessel functions
#once calculated
c = ['tab:blue', 'tab:orange', 'tab:green']
# x values for plotting Bessel functions
x = np.linspace(0.0, 20.0, 1000)

def func(phi, n, x):
    """
    Curve under question
    In particular, this function defines the intgrand
    for the Bessel functions
    """
    return np.cos(n * phi - x * np.sin(phi))

def Jnx(N, n, x):
    """
    This function calculates the order-n
    Bessel function of the first kind using
    Simpson's rule for N steps
    """
    #define interval for phi in which 
    #the Bessel function is defined
    phi_a = 0.0
    phi_b = np.pi

    h = (phi_b - phi_a) / N

    #Use simpson's rule
    I = func(phi_a, n, x) + func(phi_b, n, x)

    #for odd terms
    for k in range(1, N, 2):
        I = I + 4 * func(phi_a + k * h, n, x)

    #for even terms
    for k in range(2, N, 2):
        I = I + 2 * func(phi_a + k * h, n, x)
    
    #return the array of the nth order Bessel function
    return (1/3) * h * I * (1/np.pi)


#This block of code plots my Bessel functions 
plt.figure(1)
for i in range(len(n)):
    plt.plot(x, Jnx(N, n[i], x), c[i], label='Order - ' + str(n[i]))
plt.xlabel(r'$x$')
plt.ylabel(r'$J_{n}(x)$')
plt.title('Order-n Bessel Functions of the first kind of Jnx()')
plt.legend()
plt.grid()
plt.xlim((0,20))
plt.savefig('Bessel Graph.jpg')
plt.show()


#This block of code plots my Bessel functions along side the
#Bessel functions calculated using scipy.special
plt.figure(2)
for i in range(len(n)):
    plt.plot(x, Jnx(N, n[i], x), c[i], label='Jnx() Order - ' + str(n[i]))
    plt.plot(x, special.jv(n[i], x), 'b--', label='Scipy Order - ' + str(n[i]))
plt.xlabel(r'$x$')
plt.ylabel(r'$J_{n}(x)$')
plt.title('Bessel Function comparison between special.jv() and Jnx()')
plt.legend()
plt.grid()
plt.xlim((0,20))
plt.savefig('Bessel Graph special.jpg')
plt.show()

#this block of code plotes the fractional error of Jnx()

plt.figure(3)
for i in range(len(n)):
    #this line sets jv() to 1.0 where it is 0.0 to avoid a math error
    jv = np.where(special.jv(n[i], x) == 0.0 ,1.0, special.jv(n[i], x))
    plt.plot(x, np.divide(Jnx(N, n[i], x) - special.jv(n[i], x), jv), c[i], label='Jnx() Order - ' + str(n[i]))
    
plt.xlabel(r'$x$')
plt.ylabel(r'$J_{n}(x)$ error')
plt.title('Bessel Functions - Fractional Error of Jnx()')
plt.legend()
plt.grid()
plt.xlim((0,20))
plt.savefig('Bessel Graph fractional error.jpg')
plt.show()

#this block plots the raw difference between the Jnx() and jv()

plt.figure(4)
for i in range(len(n)):
    plt.plot(x, Jnx(N, n[i], x) - special.jv(n[i], x), c[i], label='Jnx() Order - ' + str(n[i]))
plt.xlabel(r'$x$')
plt.ylabel(r'$J_{n}(x)$')
plt.title('Difference between Jnx() and jv()')
plt.legend()
plt.grid()
plt.xlim((0,20))
plt.savefig('Bessel Graph difference.jpg')
plt.show()
