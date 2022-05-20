# This program was written by Brendan Halliday
# for PHY407 Lab11 Question 4
# This program is a modification of Nico Grisouard's 
# L11-Ising1D-start.py. 
# This program calculate the energy and magnetization
# of a 2D array of spin up and spin down particles 
# and it uses the Metroplois algorithm to perform a
# Monte Carlo simulation of magnetization 

# import modules
import numpy as np
import random
import matplotlib.pyplot as plt

# Part A
def energyfunction2D(J_, dipoles):
    """ 
    Function to calculate energy for 2D array
    This function was modified for a 2D array by Brendan
    from Nico Grisouard's 1D code.
    """
    energy = -J_*(np.sum(dipoles[:,0:-1]*dipoles[:,1:]) + np.sum(dipoles[0:-1,:]*dipoles[1:,:]))

    return energy

def generate_dipoles(N):
    """
    Generate an N x N matrix of spin up
    or spin down dipoles randomly
    """

    DIPOLES = np.zeros([N, N], int)

    for i in range(N):
        for j in range(N):
            z = random.randint(0,1)
            if z == 0:
                DIPOLES[i][j] = -1
            else:
                DIPOLES[i][j] = 1

    return DIPOLES


def acceptance(Enew, Eold, kB, T):
    """ 
    Function for acceptance probability; to be completed 
    """

    p = np.exp(-(Enew - Eold)/(kB*T)) # calculate the boltzmann value

    if Enew - Eold <= 0: # if the new energy is less than and or equal to the old
        result = True 

    elif Enew - Eold > 0 and p > random.random(): 
        # is the new energy is greater than the old and boltzmann value is greater 
        # than sum random float uniformly distributed between 0 and 1
        result = True

    else:# if neither condition holds, reject
        result = False
    
    return result  # result is True of False


if __name__ == "__main__":
    # Mech a grid of points for the animation
    x = np.arange(1, 21, 1)
    y = np.arange(1, 21, 1)
    X, Y = np.meshgrid(x,y)
    # define constants
    kB = 1.0
    T = [1.0, 2.0, 3.0]
    J = 1.0
    num_dipoles = 20
    N = 100000

    for t in T:
        # generate array of dipoles and initialize diagnostic quantities
        dipoles = generate_dipoles(num_dipoles)

        energy = []  # empty list; to add to it, use energy.append(value)
        magnet = []  # empty list; to add to it, use magnet.append(value)



        E = energyfunction2D(J, dipoles) # calculate the energy
        energy.append(E) # append this energy to rhe list
        magnet.append(np.sum(dipoles)) # calculate magnetization and append it to the list

        # This section is an implementation of the Metropolis algorithm
        # for N steps

        for i in range(N):

            E = energy[-1] # initialize the last energy value

            picked_i = random.randrange(num_dipoles)  # choose a victim
            picked_j = random.randrange(num_dipoles)
            dipoles[picked_i][picked_j] *= -1  # propose to flip the victim

            Enew = energyfunction2D(J, dipoles)  # compute Energy of proposed new state

            # calculate acceptance probability
            accepted = acceptance(Enew, E, kB, t)

            if accepted: # literally: if move is accpeted
                energy.append(Enew) # append the new energy and magnetization
                magnet.append(np.sum(dipoles))
            else: # if not...
                dipoles[picked_i][picked_j] *= -1 # flip state back to original value
                # recalculate old energy and magnetizations
                Eold = energyfunction2D(J, dipoles) 
                energy.append(Eold)
                magnet.append(np.sum(dipoles))

            # This next part is a simple animation using the python GUI 
            if i % 250 == 0:
                plt.clf() # clear the plot
                plt.scatter(X, Y, c=dipoles) # plot the current frame
                plt.title("Animation of Dipole array magnetization for T = " + str(t))
                plt.xlabel("X coordinate")
                plt.ylabel("Y coordinate")
                plt.xlim(0, 21)
                plt.ylim(0, 21)
                plt.colorbar()
                plt.draw()
                plt.pause(0.01) # pause to allow a smooth animation

        # plot energy, magnetization
        plt.figure()
        plt.plot(energy)
        plt.xlabel("Time in Steps")
        plt.ylabel("Energy (J)")
        plt.title("Energy Plot for N = " + str(N) + " steps, T = " + str(t))
        plt.grid()
        plt.savefig("Energy" + str(N) +"_"+str(t)+ ".png", dpi = 300, bbox_inches = "tight")

        plt.figure()
        plt.plot(magnet)
        plt.xlabel("Time in Steps")
        plt.ylabel("Magnetization")
        plt.title("Magnetization Plot for N = " + str(N) + " stpes, T = " + str(t))
        plt.grid()
        plt.savefig("Magnet" + str(N) +"_"+str(t)+ ".png", dpi = 300, bbox_inches = "tight")
