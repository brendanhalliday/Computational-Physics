from random import random,randrange
from math import exp,pi
from numpy import ones
import matplotlib.pyplot as plt
import os 

directory_name = os.path.dirname(__file__) # find directory of current running python file
os.chdir(directory_name) # change the working directory to the directory that includes the running file

N = 1000
T = [10.0, 40.0, 100.0, 400.0, 1200.0, 1600.0]
for t in T:
    steps = 250000

    # Create a 2D array to store the quantum numbers
    n = ones([N,3],int)

    # Main loop
    eplot = []
    E = 3*N*pi*pi/2
    for k in range(steps):

        # Choose the particle and the move
        i = randrange(N)
        j = randrange(3)
        if random()<0.5:
            dn = 1
            dE = (2*n[i,j]+1)*pi*pi/2
        else:
            dn = -1
            dE = (-2*n[i,j]+1)*pi*pi/2

        # Decide whether to accept the move
        if n[i,j]>1 or dn==1:
            if random() < exp(-dE/t):
                n[i,j] += dn
                E += dE

        eplot.append(E)

    # Make the graph
    plt.plot(eplot)
    plt.ylabel("Energy")
    plt.xlabel("Step Number")
    plt.savefig("Lab11Q1plot1.png", dpi = 300, bbox_inches = "tight")

    # this next part is from Nico Grisouard

    # This calculates the energy of each particle, neglecting constant factors
    energy_n = n[:, 0]**2 + n[:, 1]**2 + n[:, 2]**2 # n[:, 0]**2 => every row from zeroth column

    # This calculates the frequency distribution and creates a plot
    plt.figure(2)
    plt.clf()
    hist_output = plt.hist(energy_n, 50)
    plt.xlabel()
    plt.ylabel()

    # This is the frequency distribution
    energy_frequency = hist_output[0]

    # This is what the x-axis of the plot should look like
    # if we plot the energy distribution as a function of n
    # the 2nd axis of hist_output contains the boundaries of the array.
    # Instead, we want their central value, below.
    energy_vals = 0.5*(hist_output[1][:-1] + hist_output[1][1:])
    n_vals = energy_vals**0.5

    # Create the desired plot
    plt.figure(3)
    plt.clf()
    plt.bar(n_vals, energy_frequency, width=0.1)
    plt.xlabel(r"$n = \sqrt{e_{n}}$")
    plt.ylabel("Number of Particles")
    plt.title("Energy Frequency Histogram | T = " + str(int(t)))
    plt.savefig("Lab11Q1plot2_" + str(int(t)) + ".png", dpi = 300, bbox_inches = "tight")

plt.show()

