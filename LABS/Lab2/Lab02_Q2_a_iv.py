# This program was written by Brendan Halliday
# For Lab 2 question 2 part a iv
# This program simply calculates the 
# practical error estimation using the trapezoidal rule

import numpy as np


def func(x):
    """
    Curve under question
    """
    return 4 / (1 + x ** 2)

PI = np.pi
a = 0.0
b = 1.0
#define N1 and N2 
N_1 = 16
N_2 = 2 * N_1

h_1 = (b - a) / N_1
h_2 = (b - a) / N_2

#trapezoidal rule
I_1 = 0.5 * func(a) + 0.5 * func(b)
for k in range(1, N_1):
    I_1 = I_1 + func(a + k * h_1)
I_1 = h_1 * I_1


I_2 = 0.5 * func(a) + 0.5 * func(b)
for k in range(1, N_2):
    I_2 = I_2 + func(a + k * h_2)
I_2 = h_2 * I_2

print('Trapezoidal Rule:')
print('I_1 is :')
print(I_1)
print('I_2 is : ')
print(I_2)
# calculate estimation of error
error = (1/3)*(I_2 - I_1)
print('Practical error estimation:')
print(error)