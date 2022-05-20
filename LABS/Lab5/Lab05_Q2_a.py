# This program was written by Brendan Halliday
# For PHY407 Lab 5 question 2 part a.

# This program simulates the relativistic spring for various initial stretches
# and then takes the Fourier transform of this
# and plotes it in frequency domain.


import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c # speed of light
import os

if __name__ == "__main__":

    directory_name = os.path.dirname(__file__) # find directory of current running python file
    os.chdir(directory_name) # change the working directory to the directory that includes the running file
    
    # defining constants
    m = 1 # in kg mass of spring
    k = 12 # in N/m spring constant
    # initialize list of various initial positions
    x_c = c * np.sqrt(m / k)
    x_0 = [1.0, x_c, 10.0 * x_c] # in metres

    x_0str = ['1', r'$x_{c}$', r'$10x_{c}$'] # This is just for labeling each graph
    # set dt time step and number of steps for each 
    dt = [1.e-4, 1.e-4, 1.e-4] # time increment in seconds
    N = [200000, 220000, 1200000] # number of time increments

    U_1 = [] # list for each position array
    U_2 = [] # list for each velosity array


    
    #for k in range(len(x_0)):
    for j in range(3):
        t = dt[j]*np.arange(0, N[j])

        # initialize emtpy lists for position and velocity
        u_1 = 0 * t
        u_2 = 0 * t
        #set initial stretch as x_c[j]:
        u_1[0] = x_0[j]

        # range is within N[j] - 1 steps since the initial 
        # condition as already been set
        for i in range(N[j]-1):
            # first update the velocity
            u_2[i+1] = u_2[i] - dt[j] * (k / m) * u_1[i] * (1 - u_2[i]**2 / c**2)**(3/2)
            # now update the position using the updated velocity
            u_1[i+1] = u_1[i] + dt[j] * u_2[i+1]
        
        # append these arrays to bigger list
        U_1.append(u_1) 
        U_2.append(u_2)
        # This block plots the position evolution for a relativistic
        # spring for each initial position.
        plt.plot(t, u_1)
        plt.xlabel('time (s)')
        plt.ylabel('$x$ (m)')
        plt.title("Time evolution for relativistic spring with $x_{0}$ = " + x_0str[j] + " meters", y=1.08)
        plt.grid()
        plt.xlim(0, None)
        plt.savefig(x_0str[j] + "spring.png", dpi = 300)
        plt.show()
    
    
    # Calculate the normalized Fourier transform of position
    for u in range(len(U_1)):
        N_points = len(U_1[u])
        # length of interval
        T = N_points*dt[u]
        # convert to (angular) frequency domain
        w_freq = np.arange(N_points/2+1) / T
        # calculate fourier coefficients
        FC = np.abs(np.fft.rfft(U_1[u]))
        FC_max = np.max(FC)
        NORM = FC / FC_max
        plt.plot(w_freq, NORM, label = x_0str[u])
    # plot normalized Fourier transform of position
    plt.xlabel('Angular frequency (Hz)')
    plt.ylabel(r'$|\hat{x}(\omega)|/|\hat{x}(\omega)|_{max}$')
    plt.title("Normalized Fourier Transform of position")
    plt.grid()
    plt.xlim(0, 3)
    plt.legend()
    plt.savefig("FCspring.png", dpi = 300)
    plt.show()
    # calculate normalized Fourier transform of velocity
    for u in range(len(U_2)):
        N_points = len(U_2[u])
        # length of interval
        T = N_points*dt[u]
        # convert to (angular) frequency domain
        w_freq = np.arange(N_points/2+1) / T
        # calculate fourier coefficients
        FC = np.abs(np.fft.rfft(U_2[u]))
        FC_max = np.max(FC)
        NORM = FC / FC_max
        plt.plot(w_freq, NORM, label = x_0str[u])

    # plot normalized Fourier transform for velocity
    plt.xlabel('Angular frequency (Hz)')
    plt.ylabel(r'$|\hat{v}(\omega)|/|\hat{v}(\omega)|_{max}$')
    plt.title("Normalized Fourier Transform of velocity")
    plt.grid()
    plt.xlim(0, 3)
    plt.legend()
    plt.savefig("FCspringvelocity.png", dpi = 300)
    plt.show()

    T1 = [1.8, 2.2, 11.5] # These are estimated periods from the graph outputted in lab 3
                         # for the period of the relativistic spring 

    for u in range(len(U_1)):
        N_points = len(U_1[u])
        # length of interval
        T = N_points*dt[u]
        # convert to (angular) frequency domain
        w_freq = np.arange(N_points/2+1) / T
        # calculate fourier coefficients
        FC = np.abs(np.fft.rfft(U_1[u]))
        FC_max = np.max(FC)
        NORM = FC / FC_max
        plt.plot(w_freq, NORM, label = x_0str[u])
        plt.vlines(float(1.0/(T1[u])), ymin=0, ymax=1, colors="black",label=x_0str[u]  + " frequency")
    # plot normalized Fourier transform of position
    plt.xlabel('Angular frequency (Hz)')
    plt.ylabel(r'$|\hat{x}(\omega)|/|\hat{x}(\omega)|_{max}$')
    plt.title("Normalized Fourier Transform of position")
    plt.grid()
    plt.xlim(0, 3)
    plt.legend()
    plt.savefig("FCspringfrequency.png", dpi = 300)
    plt.show()