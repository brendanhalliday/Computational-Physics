# This program was written by Brendan Halliday
# For Lab 2 question 2 part a iv
# This program times both methods to see
# how long it takes to achieve and error of less thatn 10e-9
# it outputs the time for each rule along with the number
# of slices required for each rule


import numpy as np
from time import time

def func(x):
    """
    Curve under question
    """
    return 4 / (1 + x ** 2)

#deine import constants
PI = np.pi
a = 0.0
b = 1.0
N = [4,8,16,32,64,128,256,512,1024,2048,4096,8196,16384,32768,65536,131072]
delta = 1e-9

# Trapezoidal rule timing
start_t = time()
for i in N:
    h = (b - a)/i
    s_t = 0.5 * func(a) + 0.5 * func(b)

    for k in range(1, i):
        s_t = s_t + func(a + k * h)
    
    s_t = h * s_t

    # test if value value is close enough
    
    if abs(PI - s_t) < delta:
        end_t = time()
        #print the number of steps required to
        #find approximate value
        print('Number of slices required for Trapezoidal Rule is N = ', i)
        break
    else:
        continue


diff_t = end_t - start_t
print('Trapezoidal rule time: ')
print(diff_t, ' s')
print()

# Simpson's rule timing
start_s = time()
for i in N:
    h = (b - a)/i
    s_s = func(a) + func(b)

    #for odd terms
    for k in range(1, i, 2):
        s_s = s_s + 4 * func(a + k * h)

    #for even terms
    for k in range(2, i, 2):
        s_s = s_s + 2 * func(a + k * h)
    
    s_s = (1/3) * h * s_s


    # test if value value is close enough
    if abs(PI - s_s) < delta:
        end_s = time()
        #print the number of steps required to
        #find approximate value
        print("Number of slices required for Simpson's rule is N = ", i)
        break
    else:
        continue

diff_s = end_s - start_s
print("Simpson's rule time: ")
print(diff_s, ' s')

