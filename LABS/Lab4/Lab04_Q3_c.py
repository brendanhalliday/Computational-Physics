# This program was written by Brendan Halliday
# for Lab 4 Question 3 part c which asks to
# answer exercise 6.13 parts a and b from
# 'Computational Physics' by Mark Newman 
# This program defines a binary search function and
# uses it to solve the non linear equation
# derived from the blackbody equation

import numpy as np
from numpy.core.numeric import count_nonzero
from scipy import constants
import Lab04_Q3_a as q3a
import matplotlib.pyplot as plt

# Define helper functions
def f(x):
    """
    Define the nonlinear equation 
    under question.
    """

    return 5 * np.exp(-x) + x - 5

def f_0(x):
    """
    This function is to be used for 
    relaxation method for x = g(x) where
    x = x and g(x) = 5 - 5 * np.exp(-x)
    """
    return 5 - 5 * np.exp(-x)

def f_prime(x):
    """
    This function is the analytical derivative
    of f(x) for each value x
    """
    return 1 - 5 * np.exp(-x) 

def binary_search(f, x_1, x_2, threshold):
    """
    This function uses binary search to
    approximate the zero of a nonlinear equation
    """
    dx = x_2 - x_1
    i = 0
    while dx > threshold:
        i = i + 1
        x_mid = (x_2 - x_1)/2 + x_1
        if np.sign(f(x_2)) != np.sign(f(x_mid)):
            x_1 = x_mid
            x_2 = x_2
        elif np.sign(f(x_1)) != np.sign(f(x_mid)):
            x_1 = x_1
            x_2 = x_mid
        else:
            print("Choose x_1 and x_2 on opposite sides of the zero")
        dx = x_2 - x_1
    return [x_mid, i]



def newtons_method(f, f_prime, x, threshold):
    """
    This function uses Newton's Method
    to approximate the zero of the function
    """
    dx = 1
    i = 0
    while dx > threshold:
        i = i + 1
        xp = x - (f(x)/f_prime(x))
        dx = abs(xp - x)
        x = xp
    return i

def relax_count(f, threshold, x, dx):
    """
    This function calculates the number
    of iterations required for the relaxation
    mehtod to converge to a the true answer within
    a given threshold.
    """
    x_list = [x]
    i = 0
    while dx > threshold:
        x_list.append(f(x_list[-1]))
        dx = np.abs(x_list[-1]-x_list[-2])
        i += 1
    return i


if __name__ == "__main__":

    threshold = 1e-6
    x_2 = 10
    x_1 = 1
    x_b = binary_search(f, x_1, x_2, threshold)[0]
    
    

    print("The zero is locates at: ")
    print("x = ", x_b)

    # displacemetn constant
    h = constants.Planck
    k = constants.Boltzmann
    c = constants.speed_of_light
    b = (c * h) / (k * x_b)
    print("Wien Displacement Constant is: ")
    print("b = ", b)

    # Part c
    # define lambda from the sun
    LAMBDA = 502e-9
    T = b / LAMBDA

    print("The temperature of the sun is estimated to be: ")
    print(T, " K")

    # this next part calculates the number of iterations it takes 
    # each method to converge on the true value
    # it then plots number of iterations as a function of initial guess x
    count_b = []
    
    count_n = []
    
    count_r = []
    
    X = np.linspace(5, 50,100)
    plt.figure(1)



    for i in X:
        count_b.append(binary_search(f, x_1, i, threshold)[1])
    
        count_n.append(newtons_method(f, f_prime,i, threshold))
    
        count_r.append(relax_count(f_0, threshold, i, 1.0))

    
    plt.plot(X, count_b, label="Binary")
    plt.plot(X, count_n, label="Newton")
    plt.plot(X, count_r, label="Relaxation")

    plt.xlabel('Choices for initial x')
    plt.ylabel('Number of iterations')
    plt.title('Iteration comparison between Newtons Method, Binary Search and Relaxation')
    plt.grid()
    plt.xlim(5,50)
    plt.legend()
    plt.savefig('Lab4_Q3c.png', dpi=300)
    plt.show()
    





