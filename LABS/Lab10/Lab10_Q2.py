# This program was written by Brendan Halliday
# for PHY407 Lab10 question 2



import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

directory_name = os.path.dirname(__file__) # find directory of current running python file
os.chdir(directory_name) # change the working directory to the directory that includes the running file

def model_f(x, a, b):
    """
    Define a model function for curve fitting
    """
    return a * x + b

def line_fit(x_data, y_data):
    """
    This function uses curve_fit from scipy.optimize
    to perform linear regreesion on an array of data y_data
    and x_array
    """

    parameter_optimized, p_cov = curve_fit(model_f, x_data, y_data)

    a = parameter_optimized[0]
    b = parameter_optimized[1]


    return a, b


def get_tau_step(): # written by NG
    """ Calculate how far a photon travels before it gets scattered.
    OUT: optical depth traveled """
    delta_tau = -np.log(np.random.random()) # positive
    return delta_tau

def emit_photon(tau_max): # written by NG
    """ Emit a photon from the stellar core.
    IN: tau max is max optical depth
    OUT:
    tau: optical depth at which the photon is created
    mu: directional cosine of the photon emitted """
    tau = tau_max
    delta_tau = get_tau_step()
    mu = np.random.random() # no chance of a emitted photon going directly into the core
    return tau - delta_tau*mu, mu

def scatter_photon(tau): # written by NG
    """ Scatter a photon.
    IN: tau, optical depth of the atmosphere
    OUT:
    tau: new optical depth
    mu: directional cosine of the photon scattered """
    delta_tau = get_tau_step()
    mu = 2*np.random.random()-1 # sample mu uniformly from -1 to 1
    return tau - delta_tau*mu, mu

def middle_bin(bins):
    """
    Calculates the value in the middle of the bin
    since plt.hist does not return this value
    """
    bin = []
    for i in range(0, len(bins)-1):
        b = (bins[i + 1] - bins[i])/2 + bins[i]
        bin.append(b)

    return bin

def luminosity_line(n, bins):
    """
    This value calculates the "experimental"
    values for AI(mu)/I_1 where a is some constant
    """
    # calculate I_1
    I_1 = n[-1]/bins[-1]

    n = np.array(n)
    bins = np.array(bins)

    I_mu = np.divide(n,bins)
    I = I_mu/I_1

    a, b = line_fit(bins, I)

    return a, b, I

# Part A
# Here, we want to write a function that akes tau_max
# as input and then emits a photon

def photon_path(tau_max):
    """

    This function tracks the path of a photon
    by updating the optical depth position for each scattering
    event. If the photon scatters back into the core, a new 
    photon is emitted and we start a new path. This is repeated
    until a photon escapes.

    It returns
    """

    tau, mu = emit_photon(tau_max) # emit first photon
    counter = 0 # initialize a counter for the number of times the photon scatters

    while tau >= 0:

        tau, mu = scatter_photon(tau) # update tau

         # check if it has been bounced back into the core or not
        if tau <= tau_max: # if not...
            counter = counter + 1 # count valid scattering 

        else: # if tau > tau_max...reinitialize the 1D path and start again
            tau, mu = emit_photon(tau_max) # reemit a photon
            counter = 0 # reset counter to zero


    return mu, counter
# End of Part A


if __name__ == "__main__":

    # Part B

    tau_max = 10 # intitialize tau_max
    PHOTONS = 100000 # number of photons to be emitted
    MU = [] # initialize final mu values
    for n in range(PHOTONS): # simulate for N photons
        mu, counter = photon_path(tau_max)
        MU.append(mu) # append to the final mu values list

    # this next bart plots a histogram as a function of mu 
    # for n_bins 
    n_bins = 20
    plt.figure(1)
    n, mu, patches = plt.hist(MU, bins = n_bins)
    plt.xlabel(r"$\mu$")
    plt.ylabel("Number of photons")
    plt.title(r"Histogram of photons as a fucntion os escape $\mu$")
    plt.xlim(0,1)
    plt.grid()
    plt.savefig("Photon_Histogram.png", dpi = 300)
    plt.show()


    mu = np.array(middle_bin(mu)) # calculate the middle value of the bins
    n = np.array(n) 
    a, b, I = luminosity_line(n, mu) # do linear regression to find specific luminosity

    # plot theoretical curve and simulated curve
    plt.figure(2)
    plt.plot(mu, model_f(mu, 0.6, 0.4), label="Theoretical Curve")
    plt.plot(mu, model_f(mu, a, b), label="Linear Regression")
    plt.scatter(mu, I, label="Simulated points")
    plt.xlabel(r"$\mu$")
    plt.ylabel(r"$\frac{I(\mu)}{I_1}$")
    plt.title("Linear Regression")
    plt.grid()
    plt.legend()
    plt.savefig("Linear_Regression.png", dpi = 300)
    plt.show()
