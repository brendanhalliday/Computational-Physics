# This program was written by Brendan Halliday
# for PHY407 lab10 question 1
# This program uses monte carlo integration to 
# apprimate area of land on earth. 

# load libraries
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.interpolate import RegularGridInterpolator
plt.rcParams["font.family"] = "Times New Roman"


N = [50,500,5000,50000] # number of points

A = 4 * np.pi # total area of are if radius is unity

# load data
loaded = np.load('Earth.npz')
data = loaded['data']
lon_array = loaded['lon']
lat_array = loaded['lat']

def theta(z):
    """
    Generates random theta from uniform distribution z
    """

    return np.arccos(-2*z + 1)

def phi(z):
    """
    Generates angle phi from uniform distribution z
    """
    return 2*np.pi*z

def convert_theta(lat):
    """
    Convert lat coordinate to theta coordinate
    """
    return (lat + 90)*(np.pi/180)

def convert_phi(lon):
    """
    Convert lon coordinate to phi coordinate
    """
    return (lon + 180)*(np.pi/180)

def convert_lat(theta):
    """
    Convert theta back to latitude cooridinate
    """
    return (theta - np.pi/2)*(180/np.pi)

def convert_lon(phi):
    """
    Convert phi back to longitudinal coordinate
    """
    return (phi -np.pi)*(180/np.pi)

def random_points(N):
    """
    This function generates random points on the
    surface of a sphere with an even distribution
    """
    PHI = []
    THETA = []


    for n in range(N):

        z = random.random()
        PHI.append(phi(z))

        z = random.random()
        THETA.append(theta(z))

    PHI = np.array(PHI)
    THETA = np.array(THETA)
    return PHI, THETA

def numerical_int(data, lon, lat):
    """
    This function sums all the points in data
    which gives a count for the number of land points.
    This is then divided by the total number of points to give
    a numerical approximation for the fraction of earth 
    that is land.
    """
    sum = 0
    for i in range(len(lon)):
        for j in range(len(lat)):
            sum = sum + data[i][j]
    fraction = sum / (len(lon)*len(lat))
    return fraction


if __name__ == "__main__":

    # Part B

    PHI, THETA = random_points(N[2])
    
    # plot random points
    plt.figure(1)
    plt.scatter(convert_lon(PHI),convert_lat(THETA), marker="+")
    plt.xlabel("Longitude (in degrees)")
    plt.ylabel("Latitude (in degrees)")
    plt.title("Random Points N = 5000")
    plt.xlim(-180, 180)
    plt.ylim(-90, 90)
    plt.grid()
    plt.savefig("Randompoints.png", dpi = 300, bbox_inches = "tight")

    # now for each point that was plotted, we want to see if it lands 
    # on land or in the ocean
    
    # convert latitude and longitude to theta and phi respectively
    lon = convert_phi(lon_array)
    lat = convert_theta(lat_array)

    # Part C
    print("This is the numerical value of the land fraction: ")
    numerical = numerical_int(data, lon_array, lat_array)
    print(numerical)
    print("This is the numerical value for the land area with r = 1:")
    print(numerical*A)

    # Part D

    # This litle block fixes the problem where points sometimes land
    # outside the bounds of the problem
    lon = np.concatenate((lon,2*np.pi), axis = None) # add the value 2pi to the end of lon
    sliver = np.zeros(len(lat)) # now set each latitude value for 2pi to zero
    data = np.vstack((data, sliver)) # concatonate this to the data set

    # define a interpolator to find nearest point if random point doesn't
    # land on grid point
    interp = RegularGridInterpolator((lon, lat), data, method='nearest')

    for n in range(len(N)):
        k = 0 # initialize an accumulator
        PHI, THETA = random_points(N[n]) # generate random points

        # initialize arrays for storing land and ocean coordinates
        PHI_green = [] 
        THETA_green = []
        PHI_blue = []
        THETA_blue = []

        for i in range(N[n]):
     
            # find the nearest grid point in (lon, lat)
            nearest_value = interp([PHI[i], THETA[i]])

            if int(nearest_value[0]) == 1: # if this is equal to 1, it is a land point
                k = k + 1 # add one to the accumulator
                PHI_green.append(PHI[i]) # append coordinates of land point
                THETA_green.append(THETA[i])
            
            else:# if not a land point
                k = k # do nothing to the accumulator
                PHI_blue.append(PHI[i]) # append coordinate of ocean point
                THETA_blue.append(THETA[i])


        PHI_green = np.array(PHI_green)
        THETA_green = np.array(THETA_green)
        PHI_blue = np.array(PHI_blue)
        THETA_blue = np.array(THETA_blue)

        PHI_green = convert_lon(PHI_green)
        THETA_green = convert_lat(THETA_green)
        PHI_blue = convert_lon(PHI_blue)
        THETA_blue = convert_lat(THETA_blue)

        fraction = k/N[n]

        I = round(k*A/N[n], 4)

        print("This is the approximate fraction for N = " + str(N[n]))
        print(fraction)
        print("This is the area that is land if the earth has radius 1: ")
        print(I)

        # plot the points as different colours
        plt.figure(n + 2)
        plt.scatter(PHI_green, THETA_green, marker=".", s = 1.9, color="green") # plot land points
        plt.scatter(PHI_blue, THETA_blue, marker=".", s = 1.9, color="blue") # plot ocean points
        plt.xlabel("Longitude (in degrees)")
        plt.ylabel("Latitude (in degrees)")
        plt.title("Monte Carlo Simulation for N = " + str(N[n]) + ", P = k/N = " + str(fraction))
        plt.xlim(-180, 180)
        plt.ylim(-90, 90)
        plt.savefig("MC_" + str(N[n]) + ".png", dpi = 300, bbox_inches = "tight")

    plt.show()