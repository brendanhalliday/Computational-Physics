# This program was written by Brendan Halliday
# for lab 6 question 2 a

# Here we evaluate a second order system
# of equations describing the motion of
# each floor of an N storey building.

# We also evaluate for the normal modes of
# the system and then initialize the system
# in those modes to plot the motion of the building

import numpy as np
import matplotlib.pyplot as plt
import os

def verlet(v_0, r_0, A, time, dt):
    """
    This function is an algorithm for the 
    verlet method. This function takes in
    two initial vectors for position and velocity.
    It then outputs a list of vectors of r and 
    v for each timestep.
    """
    # Empty lists to store updated vectors for each time step for
    V = [] # Velocity
    R = [] # and Position.

    v = v_0 + 0.5 * dt * A.dot(r_0) # initialize velocity vector. Here we multiply r_0 by 
    # a vector valued function or in our case, the N x N stiffness matrix.

    V.append(v) # append velocity vector to V
    r = r_0 # initialize position vector
    R.append(r) # append position vector to R

    for i in range(len(time)-1):
        r = r + dt * v # new position vector 
        k = dt * A.dot(r) # new k vector
        v = v + k # new velocity vector
        # Remember that each element in r and v
        # correspond to the position and velocity
        # of a specific floor...
        R.append(r)
        V.append(v)
        # ...whereas each array in R and V correspond to
        # the entire motion of the building for a specified
        # time.

    return R, V # Each array in R and V are vecotrs for each time step

def A(k_over_m, N):
    """
    Define stiffness matrix of dimension
    N x N and a  given k/m ratio.
    """
    A = []
    # This will be the first row
    a = np.zeros(N)
    a[0] = -2.
    a[1] = 1.
    A.append(a)
    
    for i in range(N - 2):
        # generate the next row so long as that row
        # is not the bottom row
        a = np.zeros(N)
        a[i] = 1.
        a[i + 1] = -2.
        a[i + 2] = 1.
        A.append(a)
    # generate the last row (similar to the first row)
    a = np.zeros(N)
    a[-1] = -2.
    a[-2] = 1.
    A.append(a)
    return (k_over_m) * np.array(A) # return the array mulitplied by k/m

if __name__ == "__main__":
    directory_name = os.path.dirname(__file__) # find directory of current running python file
    os.chdir(directory_name) # change the working directory to the directory that includes the running file

    # store important constants
    k_over_m = 400 # rad/s k/m
    dt  = 0.001 # seconds
    time = [np.arange(0, 800) * dt, np.arange(0, 1800) * dt] # time array
    x_0 = 0.1 # meters
    N = [3, 10] # number of floors
    plt.rcParams["font.family"] = "Times New Roman" # to make font the same as the report
    
    for n in N:
        # initial conditions
        v_initial = np.zeros(n) # state of rest
        x_initial = np.zeros(n) 
        x_initial[0] = x_0 # initial displacement from the vertical

        # calculate R and V using the verlet function
        for t in range(len(time)):
            R, V = verlet(v_initial, x_initial, A(k_over_m, n), time[t], dt)
                
            # plot the morion of the building for each n
            plt.figure()
            
            for i in range(len(R[0])):

                x_t = [] # position of a particular floor with respect to time

                for k in range(len(R)):
                    x_t.append(R[k][i]) 
                    # Here, we index and loop in this order to extract 
                    # the time evolution of  a single floor(i.e. we extract
                    # the columns of R)
                # motion of the floor over time
                plt.plot(time[t], x_t, label="Floor " + str(i + 1))
            # This is just plotting details
            plt.xlabel("time (s)")
            plt.ylabel("Position from the vertical (meters)")
            plt.title("Motion of each floor for an N = " + str(n) +" storey building", y = 1.05)
            plt.grid()
            plt.xlim(0, time[t][-1])
            plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
            ncol=3, fancybox=True, shadow=False)
            plt.savefig("Motionx" + str(n) +"time" + str(t) + ".png", dpi = 300, bbox_inches = "tight")
      
    # plot the motion for floors 1 and 10
    R, V = verlet(v_initial, x_initial, A(k_over_m, N[1]), time[1], dt)
    plt.figure()
    for i in [0, 9]:
        x_t = [] # position of a particular floor with respect to tim
        for k in range(len(R)):
            x_t.append(R[k][i])
        plt.plot(time[t], x_t, label="Floor " + str(i + 1))
    plt.xlabel("time (s)")
    plt.ylabel("Position from the vertical (meters)")
    plt.title("Motion of floor 1 and 10", y = 1.05)
    plt.grid()
    plt.xlim(0, time[1][-1])
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
    ncol=3, fancybox=True, shadow=False)
    plt.savefig("Motionx110.png", dpi = 300, bbox_inches = "tight")


# Solution for part b =======================================================================|
    # solution to eigenvalue problem
    # calculate the normal frequencies and 
    # normal modes from the eigenvalue problem
    A = A(k_over_m, N[0])
    eigenvalues = np.linalg.eigh(A)[0] # calculate raw eigenvalues
    normal_freq = np.sqrt(-eigenvalues)/(2*np.pi) # divide by 2pi to go from rad/s to Hz
    EV = []
    # this is here to extract the columns since they are
    # the eigenvectors
    for i in range(N[0]):
        eigenvectors = np.linalg.eigh(A)[1][:,i] 
        EV.append(eigenvectors)
    
    # print these solutions
    for i in range(len(EV)):
        print("Normal mode " + str(i + 1)+ ":")
        print(EV[i], " meters")
        print()
        print("With normal frequency: ")
        print(normal_freq[i], " Hz")
        print()

    # for each normal mode
    for j in range(len(EV)):

        # initial conditions
        v_initial = np.zeros(N[0]) # state of rest
        x_initial = EV[j] # normal mode

        # calculate R and V using the verlet function
        R, V = verlet(v_initial, x_initial, A, time[0], dt)
            
        plt.figure()
        # extract time series of each floor's position
        for i in range(len(R[0])):
            x_t = [] # position of a particular floor with respect to time
            for k in range(len(R)):
                x_t.append(R[k][i])

            # since there is sometimes overlap
            # in the position of different floors
            # this has been added for better legibility
            linestyle = "-"
            color = None
            if i == 2:
                linestyle ="--"
                color = 'r'
            # plot position of the floor against time
            plt.plot(time[0], x_t, linestyle = linestyle, color = color, label="Floor " + str(i + 1))

        # plottind details
        # this for loop plots the period derived from the calculated fre
        for d in range(10):
            plt.vlines(x = d/normal_freq[j], ymin = -0.8, ymax = 0.8, 
            colors = "black", linestyle="--")

        plt.xlabel("time (s)")
        plt.ylabel("Position from the vertical (meters)")
        plt.title("Normal mode " + str(j + 1) + " for 3 storey building", y = 1.05)
        plt.grid()
        plt.xlim(0, 0.8)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
        ncol=3, fancybox=True, shadow=False)
        plt.ylim(-0.8, 0.8)
        plt.savefig("normal" + str(j) + ".png", dpi = 300, bbox_inches = "tight")

    plt.show()