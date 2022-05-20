# The following program was written by Brendan Halliday
# for PHY407 lab 7. It uses bits of code from "Computational Physics"
# by Mark Newman and it uses the solution for excersice 8.8
# written by Nico Grisouard

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from time import time

# The following helper functions were written by Brendan Halliday
#------------------------------------------------------------------|
def rho(h, delta, eps):
    """
    This function calucaltes the ratio
    between the desired accuracy and 
    """
    return (h * delta)/eps

def eps(x1, x2, y1, y2):
    """
    This function calcualtes the estimated
    error of x1(approximation using 2 steps of 
    size h) and x2(approximation using 1 step of 2h)
    """
    return ((1/30) * np.abs(x1 - x2))**2 + ((1/30) * np.abs(y1 - y2))**2

def runge_kutta_2d(r, f, h):
    """
    This is a function approximates
    the solution to the ODE using runge-kutta
    given an initial condition r and then returns 
    an updated r for the next time step
    """
    
    k1 = h*f(r)
    k2 = h*f(r + 0.5*k1)  
    k3 = h*f(r + 0.5*k2)
    k4 = h*f(r + k3)
    r += (k1 + 2*k2 + 2*k3 + k4)/6
    return r
#------------------------------------------------------------------|

# Solution to Newman 8.8, Space garbage.
# Author: Nico Grisouard, Univ. of Toronto
# -----------------------------------------------------------------|
def rhs(r):
    """ The right-hand-side of the equations
    INPUT:
    r = [x, vx, y, vy] are floats (not arrays)
    note: no explicit dependence on time
    OUTPUT:
    1x2 numpy array, rhs[0] is for x, rhs[1] is for vx, etc"""
    M = 10.
    L = 2.

    x = r[0]
    vx = r[1]
    y = r[2]
    vy = r[3]

    r2 = x**2 + y**2
    Fx, Fy = - M * np.array([x, y], float) / (r2 * np.sqrt(r2 + .25*L**2))
    return np.array([vx, Fx, vy, Fy], float)

if __name__ == "__main__":
    

    # This next part adapted from Newman's odesim.py --------------------------|
    a = 0.0
    b = 10.0
    N = 10000  # let's leave it at that for now
    h = (b-a)/N

    tpoints = np.arange(a, b, h)

    xpoints = []
    vxpoints = []  # the future dx/dt
    ypoints = []
    vypoints = []  # the future dy/dt
    # below: ordering is x, dx/dt, y, dy/dt
    r = np.array([1., 0., 0., 1.], float)

    time1 = time()

    for t in tpoints:
        xpoints.append(r[0])
        vxpoints.append(r[1])
        ypoints.append(r[2])
        vypoints.append(r[3])
        k1 = h*rhs(r)  # all the k's are vectors
        k2 = h*rhs(r + 0.5*k1)  # note: no explicit dependence on time of the RHSs
        k3 = h*rhs(r + 0.5*k2)
        k4 = h*rhs(r + k3)
        r += (k1 + 2*k2 + 2*k3 + k4)/6
        
    time2 = time()

    print("Time for non-adaptive method: ")
    print(time2 - time1, "s")
    
    # plot Non-adaptive method
    plt.figure()
    plt.plot(xpoints, ypoints, ':', label="non Adaptive")
    # end of Nico Grisouard's solution --------------------------------------------|

    b = 2
    # begining of my main program -------------------------------------------------| 
    delta = 1e-6 # target error in m s^-1
    h0 = 0.01 # reset initial time step
    r0 = np.array([1., 0., 0., 1.], float) # re-initialize r
    
    # re initialize position and velocity lists
    xpoints0 = []
    vxpoints0 = []  
    ypoints0 = []
    vypoints0 = []

    t = a
    time10 = time()
    #repeat this until t equals 10
    while t <= b:
        xpoints0.append(r0[0])
        vxpoints0.append(r0[1])
        ypoints0.append(r0[2])
        vypoints0.append(r0[3])
        
        r_1 = r0 # initialize r_1 as r
        r_2 = r0 # initialize r_2 as r

        # repeat this next section for each point in time
        # while t < tpoints:
        for i in range(2): # calculates r_1 after 2 iterations of R-K using h
            r_1 = runge_kutta_2d(r_1, rhs, h0) # updates r
        x_1 = r_1[0] # extract x component 
        y_1 = r_1[2] # extract y component

        r_2 = runge_kutta_2d(r_2, rhs, 2*h0) # calculate r_2 after 1 iteration of R-K using 2h
        x_2 = r_2[0]
        y_2 = r_2[2]

        # calculate the value of rho
        p = float(rho(h0, delta, eps(x_1, x_2, y_1, y_2)))

        # now adjust step size
        # but we must ensure that h doesn't increase or decrease too much
        upper_lim = 2 * h0
        lower_lim = (1/2) * h0
        if p > 1:
            if h0 * (p**(1/4)) >= upper_lim:
                h0 = upper_lim
            else:
                h0 = h0 * (p**(1/4))
        elif p < 1: # step size is too large
            if h0 * (p**(1/4)) <= lower_lim:
                h0 = lower_lim
            else:
                h0 = h0 * (p**(1/4))

        else:
            h0 = h0

        # Now calculate x and y using RK again but using the updated h
        r0 = runge_kutta_2d(r0, rhs, h0)
        # update time using new h
        t = t + h0
    time20 = time()

    print("Time for adaptive method:")
    print(time20 - time10, " s")

    # plot adaptive method
    plt.plot(xpoints0, ypoints0, ':', label="Adaptive")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.title('Trajectory of a ball bearing around a space rod.')
    plt.axis('equal')
    plt.grid()
    plt.tight_layout()
    plt.legend()
    plt.savefig('Garbage.png', dpi=300)
    plt.show()