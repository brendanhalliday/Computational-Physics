# This program was written by Brendan Halliday
# For Lab 2 question 2 part a ii
# This program compares computes the values for 
# the integral using the Trapezoidal rule and 
# Simpson's rule. It then compares it to the 
# exact value of pi


import numpy as np

# important constants
a = 0.0
b = 1.0
N = 4
h = (b - a)/N

print('Exact Value of pi:')
print(np.pi)

def func(x):
    """
    Curve under question
    """
    return 4 / (1 + x ** 2)

def frac_error(accepted, experimental):
    """
    Calculate percent uncertainty
    """
    return (abs(accepted - experimental)/accepted)

#Trapezoidal Rule

s_t = 0.5 * func(a) + 0.5 * func(b)

for k in range(1, N):
    s_t += func(a + k * h)
print('Trapezoidal rule:')
trap = h * s_t
print(trap)
print('Fractional error:')
print(frac_error(np.pi, trap))




#Simpson's rule 

s_s = func(a) + func(b)

#for odd terms
for k in range(1, N, 2):
    s_s = s_s + 4 * func(a + k * h)
#for even terms
for k in range(2, N, 2):
    s_s = s_s + 2 * func(a + k * h)

simpson = (1 / 3) * h * s_s

print("Simpson's rule:")
print(simpson)

print('Fractional error:')
print(frac_error(np.pi, simpson))