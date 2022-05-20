import numpy as np
import matplotlib.pyplot as plt

def initial_psi(x, C, L, sigma, d):
    """
    This function defines the initial
    profile of the string's velocity.
    """
    return (C*x*(L - x)/(L**2)) * np.exp((-1/2)*((x - d)**2)/(sigma**2))

if __name__ == "__main__":
    h = 1e-6 # seconds
    v = 100 # meters per second
    C = 1.0 # meters per second
    L = 1.0 # meters
    sigma = 0.3 #meters
    d = 0.1 # meters
    M = 100
    a = L/M
    T = 100000
    X = np.arange(0, L + a, a)
    # initialize array to store spatial and time data for each point on string
    phi = np.zeros([T+1, M+1], float)
    # initialize array to store veolicty and time data for each point on string
    psi = np.zeros([T+1, M+1], float)
    psi[0,:] = initial_psi(X, C, L, sigma, d)

    for n in range(T): # for each time step
        for m in range(M): # for each x value
            phi[n + 1,m] = phi[n, m] + h * psi[n, m] #calculate the displacement for next time step
            psi[n + 1,m] = psi[n, m] + h * (phi[n, m - 1] + phi[n, m +1] - 2*phi[n, m]) * (v/a)**2 #calculate velocity for next time step

    # plt.figure(1)
    # plt.plot(X, phi[500,:], label = "t = " + str(500*h) + " s")
    # plt.plot(X, phi[1000,:], label = "t = " + str(1000*h)+ " s")
    # plt.plot(X, phi[2000,:], label = "t = " + str(2000*h)+ " s")
    # plt.plot(X, phi[3000,:], label = "t = " + str(3000*h)+ " s")
    # plt.plot(X, phi[4000,:], label = "t = " + str(4000*h)+ " s")

    # plt.title("Animation of a plucked string")
    # plt.xlabel("X (meters)")
    # plt.ylabel("Displacement (meters)")
    # plt.ylim(-0.0004, 0.0004)
    # plt.xlim(0, L)
    # plt.grid()
    # plt.legend()
    # plt.savefig("stable.png", dpi = 300, bbox_inches = "tight")

    # plt.figure(2)
    # plt.plot(X, phi[30000,:], label = "t = " + str(30000*h)+ " s")
    # plt.plot(X, phi[40000,:], label = "t = " + str(40000*h)+ " s")
    # plt.title("Animation of a plucked string")
    # plt.xlabel("X (meters)")
    # plt.ylabel("Displacement (meters)")
    # plt.ylim(-0.0004, 0.0004)
    # plt.xlim(0, L)
    # plt.grid()
    # plt.legend()
    # plt.savefig("unstable.png", dpi = 300, bbox_inches = "tight")

    # plt.figure(3)
    # plt.plot(X, phi[T,:], label = "t = " + str(40000*h)+ " s")
    # plt.title("Animation of a plucked string")
    # plt.xlabel("X (meters)")
    # plt.ylabel("Displacement (meters)")
    # plt.ylim(-0.0004, 0.0004)
    # plt.xlim(0, L)
    # plt.grid()
    # plt.legend()
    # plt.savefig("superunstable.png", dpi = 300, bbox_inches = "tight")
    
    # the following block generates an animation for the motion of a string 
    for t in range(T):
         if t % 50 == 0:
             plt.clf() # clear the plot
             plt.plot(X, phi[t,:]) # plot the current frame
             plt.title("Animation of a plucked string time = " + str(round(t*h,4)) + " s")
             plt.xlabel("X (meters)")
             plt.ylabel("Displacement (meters)")
             plt.grid()
             plt.xlim(0, L)
             plt.ylim(-0.0004, 0.0004)
             plt.draw()
             plt.pause(0.01) # pause to allow a smooth animation