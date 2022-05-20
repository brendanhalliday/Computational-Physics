# This program was written by Brendan Halliday
# for PHY407 Lab 5 question 3.

# This program plots the SLP as a function of time and 
# longitude. It also computes the fourier transform of 
# SLP and then extracts only the fourier coefficients 
# associated with longitudinal wave number 3 and 5. 
# The inverse Fourier transform is then used on the
# extracted data to be plotted as filterd data
# in time and longitudinal space.

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
import os

def speed_of_propagation(longitude, time):
    """
    This function calculates the speed of 
    propagation given a longitude in degrees
    and a time in days. This is for latitude of 50
    degrees

    """
    day_in_secs = constants.day # seconds in a day
    time = time * day_in_secs # convert time into seconds
    Re = 6.371e6 # Earth radius in meters
    pi = np.pi
    latitude = 50 * (pi / 180) # latitude into radians
    longitude = longitude * (pi / 180) # longitude into radians

    return Re*np.cos(latitude) * (longitude / time) # angular velocity * radius = linear velocity

if __name__ =="__main__":

    # Like the previous question, this must be here for the code to run
    directory_name = os.path.dirname(__file__) # find directory of current running python file
    os.chdir(directory_name) # change the working directory to the directory that includes the running file

    # These two blocks of code are from the lab manual
    #import data for plotting
    SLP = np.loadtxt('SLP.txt')
    Longitude = np.loadtxt('lon.txt')
    Times = np.loadtxt('times.txt')
    print(np.shape(SLP))

    # plot contour plot of SLP as a function
    # of time and longitude
    plt.contourf(Longitude, Times, SLP)
    plt.xlabel('longitude(degrees)')
    plt.ylabel('days since Jan. 1 2015')
    plt.title('SLP anomaly (hPa)')
    plt.colorbar()
    plt.savefig("SLP.png", dpi = 300)
    plt.show()

    # This part is not from the lab manual.
    # Repeat this next section for longitudinal wavenumber 
    # 3 and 5
    M = [3, 5]
    for i in M:
        
        A = np.fft.fft2(SLP)
            # for each row
        for k in range(len(A)):
                # For each column:
                # These m's are the various wavnumbers starting from m = 0
                # They are assosicated with a Fourier Coefficient in (1 / logitudinal) space
                # 
            for m in range(len(A[0])):
                    # longitudinal wavenumber is not i = [3, 5] then
                    # discard it. This extraction is like the filter 
                    # from the exercise 2
                if m != i:
                    A[k][m] = 0.
                else:
                    A[k][m] = A[k][m]

        # use inverse Fourier transform to convert new coefficients 
        # into longitudinal and time space
        Back = np.fft.ifft2(A)
        
        # plot each extraction
        plt.contourf(Longitude, Times, Back)
        plt.xlabel('longitude(degrees)')
        plt.ylabel('days since Jan. 1 2015')
        plt.title('SLP anomaly (hPa) for m = ' + str(i))
        plt.colorbar()
        plt.grid()
        plt.savefig("SLPFourier_m" + str(i) + ".png", dpi=300)
        plt.show()

    # This is a user input program to estimate 
    # the spped of propagation for m = 5 for
    # different times in the 120 day period
    X = "START"
    while X != 'END':
        print("Input Longitude 2(in degrees):" )
        l2 = input()
        print("Input time 2(in days):" )
        t2 = input()
        print("Input Longitude 1(in degrees):" )
        l1 = input()
        print("Input time 1(in days):" )
        t1 = input()
        print("Speed of propagation is : ")
        print(speed_of_propagation((float(l2)-float(l1)), (float(t2)-float(t1))), " m/s")

        print("To end program, type 'END':")
        print("To continue, type: 'AGAIN':")

        X = input()
    
    








