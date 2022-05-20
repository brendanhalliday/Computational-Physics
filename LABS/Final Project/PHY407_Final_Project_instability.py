# This program was written by Brendan Halliday
# as a final project for PHY407 in december of 2021

# This program is meant to solve Maxwell's equations 
# for an electromagnetic pulse in one dimenstion.
# The Yee algorithm is used to simultaneously solve the coupled wave
# equtaions of the electric and magnetic fields

import matplotlib.pyplot as plt
import numpy as np


LAMBDA_0 = 0.3 # meters

# Next we define the permeability and permitivity
# for different regions
# region 1
EPS_1 = 1
MU_1 = 1

# region 2
EPS_2 = [2, 8, 50, 1000]
MU_2 = 1

# region 3
EPS_3 = [4, 16, 100, 2000]
MU_3 = 1

DX = LAMBDA_0/20 # this will be the spatial step 
# we will need to test the differences in spatial step without
# changing time steps

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


def slab_length_epsilon(eps):
    """
    
    """
    L_1 = LAMBDA_0/(2*np.sqrt(eps)) # length of the middle dielectric
    return L_1


def slab_length(l,dx):
    """
    This function adjusts the desired length of the dielectric slab
    It takes in a desired length from 0 to 3 meters and converts it 
    to indexed bounds of the slab.
    """
    # first calculate the integer number of spatial steps
    N = int(l/dx)

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

if __name__ == "__main__":
    
    DX = LAMBDA_0/20.05 # Change in DX

    I_MAX = int(3/DX)
    T = 200
    XE = np.arange(0, 3 + I_plus*DX, DX)
    E = np.zeros((T, I_MAX + I_plus), dtype=float) 

    # For the B field we make the bounds 1 less than the E field to respect the domain of 3 meters
    B = np.zeros((T, I_MAX + I_plus - 1), dtype=float) 
    
    # find out where we need to stop the emission of the pusle
    TAU = tau(dt(v(EPS_1, MU_1)))
    STOP = 40 #int(6*TAU) + 1 # this is the index where the pulse should stop 321

    narray = np.arange(0, STOP) # create integer array for gaussian pulse
    E_0 = E_initial(narray, TAU) # generate the pulse
    E[:STOP, 0] = E_0 # merge the pulse into the space-time array
    
    L_1 = slab_length_epsilon(EPS_2[0])
    x1, x2 = slab_length(L_1, DX)
    MU = permeability_mu(MU_1, MU_2, MU_3, x1, x2)
    EPS = permittivity_ep(EPS_1, EPS_2[0], EPS_3[0], x1, x2)

    # before we implement the YEE algorithm, we need to ensure dt and v are
    # changing correctly
            
    # firstly, let us figure out v
    V = v(EPS, MU)
    DT = dt(V)

    # define electric conductivity term
    SIG = np.ones(I_MAX + I_plus)
    SIG[0:I_MAX] = 0
    SIG[I_MAX:I_MAX + I_plus] = 1000

    # define the magnetic conductivity term
    SIGM = np.zeros(I_MAX + I_plus)

        
    for n in range(0, T - 1): 
        for i in range(0, I_MAX + I_plus - 2):
                    
            # material coefficients
            C0 = (2*MU[i] - SIGM[i]*DT[i])/(2*MU[i] + SIGM[i]*DT[i])
            C1 = 2*DT[i]/((2*MU[i] + SIGM[i]*DT[i])*DX)
                    
            # firstly, update the magnetic field: 
            B[n + 1, i] = C0 * B[n, i] + C1 * (E[n, i + 1] - E[n, i])
                    
        # we use a second for loop so that the magnetic field can finish updating
        # before we update the electric field:
            
        for i in range(1, I_MAX + I_plus - 1):
                    
            # material coefficients
            C2 = (2*EPS[i] - SIG[i]*DT[i])/(2*EPS[i] + SIG[i]*DT[i])
            C3 = 2*DT[i]/((2*EPS[i] + SIG[i]*DT[i])*DX)
                    
            # now we can update the electric field:
            E[n + 1, i] = C2 * E[n, i] + C3 * (B[n + 1, i] - B[n + 1, i - 1])

        # now we need to update the boundaries
        # first let's do the right boundary condition:
        E[n + 1,-1]  = E[n,-2] + ((V[-1]*DT[-1] - DX)/(V[-1]*DT[-1] + DX))*(E[n + 1,-2] - E[n,-1])

        # then let's do the left boundary condition
        # this condition only to be implement after 
        # a certain time:
        if n >= STOP - 1:
            E[n + 1, 0] = E[n, 1] + ((V[0]*DT[0] - DX)/(V[0]*DT[0] + DX))*(E[n + 1, 1] - E[n, 0])
        else:
            pass

    for t in range(T):
        if t % 1 == 0:
            #continue 
            plt.clf() # clear the plot
            plt.plot(XE[0:-1], E[t,:]) # plot the current frame
            #plt.plot(XB, B[t,:])
            plt.title("EMP Pulse against dielectric interface time step= " + str(t))
            plt.xlabel("X (meters)")
            plt.ylabel("Electric Field")
            plt.ylim(-1, 1)
            plt.xlim(0,3)
            plt.vlines([XE[x1],XE[x2]], ymin=-1, ymax=1, color='grey')
            plt.axvspan(XE[x1], XE[x2], alpha=0.1, color='grey')
            plt.axvspan(XE[x2], XE[I_MAX], alpha=0.5, color='grey')
            plt.grid()
            plt.draw()
            plt.pause(0.0001) # pause to allow a smooth animation

    plt.clf()
    plt.figure(1)
    plt.plot(XE[0:-1], E[115,:]) # plot the current frame
    plt.grid()
    plt.xlabel("X (meters)")
    plt.ylabel("Electric Field")
    plt.title("EMP Pulse against dielectric interface time step= " + str(115))
    plt.ylim(-1, 1)
    plt.xlim(0,3)
    plt.vlines([XE[x1],XE[x2]], ymin=-1, ymax=1, color='grey')
    plt.axvspan(XE[x1], XE[x2], alpha=0.1, color='grey')
    plt.axvspan(XE[x2], XE[I_MAX], alpha=0.5, color='grey')
    plt.savefig("instability.png", dpi=300, bbox_inches = "tight")
    plt.show()
        

   