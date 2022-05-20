# This program was written by Brendan Halliday for
# Lab1 Q3b.
# This program computes the histogram
# using my histogram function and compares it
# to NumPy's np.histogram function in terms of shape
# and runtime.

#import libraries and helper module Lab1_Q3
import matplotlib.pyplot as plt
import numpy as np
from time import time


def my_bins(start, stop, num_bins):
    """This function returns an array where
    the entries define partitions or bins in
    a specified interval
    """
    return np.linspace(start, stop, num_bins + 1, endpoint=True)

def my_hist(R, bins):
    """This function takes in an array R
    of values and returns the histogram of those values
    """
    # initiallize empty list
    list = []

    # repeat the code below for each bin
    for i in range(len(bins) - 1):
        # initiallize accumulator to zero for each bin
        accumulator = 0

        # repeat the code below for each value in R
        for j in R:
            # for each j in R, check if value is within
            # the bin under question
            # if so, add one to the accumulator
            # if not, do not change accumulator value
            if bins[i] <= j < bins[i + 1] and i < len(bins) - 1:
                accumulator = accumulator + 1
            elif bins[i] <= j <= bins[i + 1] and i == len(bins) - 1:
                accumulator = accumulator + 1
            else:
                accumulator = accumulator

        # once each j has been checked, append
        # accumulator value to empty list
        list.append(accumulator)
    # repeat this whole process for each bin
    return np.array(list)

# store usleful constants
N = [10, 100, 1000, 10000, 100000, 1000000]
M = 1000

#use helper function to generate bins
bins = my_bins(-5, 5, 1000)

#initiallize empty lists for storing the time values
#for both np.histogram and my_hist()
my_TIME = []
np_TIME = []

for i in N:
    # R is an array of N random variables with
    R = np.random.randn(i)

    #time my.hist() which is my histogram program
    start_1 = time()
    hist = my_hist(R, bins)
    end_1 = time()
    diff_1 = end_1 - start_1
    #append run times to array my_TIME
    my_TIME.append(diff_1)

    # plot histogram for my_hist()
    plt.bar(x = bins[0:-1], width = 0.01 , height = hist, align = 'edge', color = 'blue')
    plt.xlabel('Bins')
    plt.ylabel('Datapoints')
    plt.title('Histogram for N = ' + str(i) + ' of my_hist()')
    plt.grid()
    plt.xlim(-5, 5)
    plt.savefig('myhistN_' + str(i) + '.jpg')
    plt.show()

    #time np.histogram() for the same R
    start_2 = time()
    nphist = np.histogram(a = R, bins = bins, density = False)[0]
    end_2 = time()
    diff_2 = end_2 - start_2
    #append times to array np.TIME
    np_TIME.append(diff_2)

    #plot histogram for np.histogram
    plt.bar(x = bins[0:-1], width = 0.01, height = nphist, align = 'edge', color = 'orange')
    plt.xlabel('Bins')
    plt.ylabel('Datapoints')
    plt.title('Histogram for N = ' + str(i) + ' of np.histogram()')
    plt.grid()
    plt.xlim(-5, 5)
    plt.savefig('nphistN_' + str(i) + '.jpg')
    plt.show()

#repeat process N times

# plot the times for my histogram calculator and np.histogram
plt.plot(N, my_TIME, 'r+', label='my_hist()')
plt.plot(N, np_TIME, 'bx', label='np.histogram()')
plt.yscale('log')
plt.xlabel('N')
plt.ylabel('Time(s)')
plt.title('Time comparison between my_hist() and np.histogram()')
plt.legend()
plt.savefig('timegraph.jpg')
plt.show()