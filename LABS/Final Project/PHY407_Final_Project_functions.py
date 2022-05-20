# This program was written by Brendan Halliday
# as a final project for PHY407 in december of 2021

# This program is meant to solve Maxwell's equations 
# for an electromagnetic pulse in one dimenstion.
# The Yee algorithm is used to simultaneously solve the coupled wave
# equtaions of the electric and magnetic fields

import numpy as np

LAMBDA_0 = 0.3 # meters

# Next we define the permeability and permitivity
# for different regions
# region 1
EPS_1 = 1
MU_1 = 1

# region 2
EPS_2 = 2
MU_2 = 1

L_1 = LAMBDA_0/(2*np.sqrt(EPS_2)) # length of the middle dielectric
L_2 = LAMBDA_0/(4*np.sqrt(EPS_2)) # length of the middle dielectric

# region 3
EPS_3 = 4
MU_3 = 1

DX = LAMBDA_0/20 # this will be the spatial step 
# we will need to test the differences in spatial step without
# changing time steps

# DX = LAMBDA_0/10
# DX = LAMDBA_0/5

X = np.arange(0,3, DX)

I_MAX = int(3/DX) # this will be the number of discrete spatial points
I_plus = 40


def v(epsilon, mu):
    """
    The velocity of our wave will change as a function
    of the medium it is propagating through
    """
    return 1/np.sqrt(epsilon * mu)

def dt(velocity):
    """
    Change time step depending on the velocity of the medium
    i.e. the velocity
    """
    return (LAMBDA_0/20)/velocity

# FOR THE INITIAL PULSE ###########################################################
def tau(d_t):
    """
    Define a tau as specified by the exercise from Inan.
    This 
    """
    return (0.8)/d_t 

def E_initial(n, tau):
    """
    Define the source as a Gaussian. This will be yeeted
    from the left boundary. This initial condition is
    based on the one given in the exercise from Inan, with
    a few modifications for better visual output in our
    scenario. 
    """
    return np.exp((-(n - 0.3*tau)**2) / (0.01*tau**2))
# FOR THE INITIAL PULSE ##########################################################

def permittivity_ep(ep1,ep2,ep3,x1,x2):
    """
    This function defines an 1D array that
    defines the permittivity of each point in 
    space. 
    """
    ep = np.ones(I_MAX + I_plus)
    ep[0:x1] = ep1
    ep[x1:x2] = ep2
    ep[x2:] = ep3
    return ep

def permeability_mu(mu1,mu2,mu3,x1,x2):
    """
    This function defines a 1D array that defines 
    the permeability of each point in space.
    """
    mu = np.ones(I_MAX + I_plus)
    mu[0:x1] = mu1
    mu[x1:x2] = mu2
    mu[x2:] = mu3
    return mu

def slab_length(l):
    """
    This function adjusts the desired length of the dielectric slab
    It takes in a desired length from 0 to 3 meters and converts it 
    to indexed bounds of the slab.
    """
    # first calculate the integer number of spatial steps
    N = int(l/DX)

    if N % 2 != 0:
        N = N - 1
        N = N/2
    else:
        N = N/2
    
    lower = I_MAX/2 - N
    upper = I_MAX/2 + N
    # upper and lower are defined in this way so when
    # you choose the bounds of of say E[:,lower:upper]
    # and it returns the array we want
    return int(lower), int(upper)