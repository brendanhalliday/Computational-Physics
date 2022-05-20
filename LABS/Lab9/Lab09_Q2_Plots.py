import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as constant

N = 3000 
P = 1024 # number of segments in x
L = 1.0
L_0 = 1.0e-8 # in meters
X = np.linspace(-L/2, L/2, P + 1) # interval of x in units of L

M = constant.electron_mass # in kg
hbar = constant.hbar # in eV(attosecond)
E_0 = (hbar**2)/(M*(L_0**2))
TIME = np.arange(N+1) #time increment in attoseconds


T_0 = 0
T_1 = int(N/4)
T_2 = int(N/2)
T_3 = N

npzfile = np.load(r'plotting_data.npz')
psi = npzfile['psi']
ENERGY = npzfile['ENERGY']
NORM = npzfile['NORM']

npzfile1 = np.load(r'QHO_data.npz')
phi = npzfile1['psi']

plt.figure(1)
plt.plot(TIME,ENERGY/E_0)
plt.xlabel(r"$\frac{t}{t_{0}}$")
plt.ylabel(r"$\frac{E}{E_{0}}$")
plt.title("Energy as a function of time", y = 1.05)
plt.grid()
plt.savefig("Energytime.png", dpi = 300, bbox_inches = "tight")

    # Plot the noramlaization constant as a function of time.
plt.figure(2)
plt.plot(TIME,NORM)
plt.xlabel("Attoseconds")
plt.ylabel("Normalization constant")
plt.title("Normalization constant as a function of time", y = 1.05)
plt.grid()
plt.savefig("Normtime.png", dpi = 300, bbox_inches = "tight")

psi0 = psi[T_0,:]
psi1 = psi[T_1,:]
psi2 = psi[T_2,:]
psi3 = psi[T_3,:]


plt.figure(3)
plt.plot(X, np.real(psi0),label="t = 0")
plt.plot(X, np.real(psi1),label="t = T/4")
plt.plot(X, np.real(psi2),label="t = T/2")
plt.plot(X, np.real(psi3),label="t = T")
plt.xlabel("Position in units of L")
plt.ylabel(r"Re($\psi$)")
plt.title("Real part of Wavefunction", y = 1.05)
plt.grid()
plt.legend()
plt.savefig("PSIreal.png", dpi = 300, bbox_inches = "tight")


plt.figure(4)
plt.plot(X, np.imag(psi0), label="t = 0")
plt.plot(X, np.imag(psi1), label="t = T/4")
plt.plot(X, np.imag(psi2), label="t = T/2")
plt.plot(X, np.imag(psi3), label="t = T")
plt.xlabel("Position in units of L")
plt.ylabel(r"Im($\psi$)")
plt.title("Imaginary part of wavefunction", y = 1.05)
plt.grid()
plt.legend()
plt.savefig("PSIimage.png", dpi = 300, bbox_inches = "tight")


psi_square0 = abs(psi0)**2
psi_square1 = abs(psi1)**2
psi_square2 = abs(psi2)**2
psi_square3 = abs(psi3)**2


plt.figure(5)
plt.plot(X, psi_square0, label="t = 0")
plt.plot(X, psi_square1, label="t = T/4")
plt.plot(X, psi_square2, label="t = T/2")
plt.plot(X, psi_square3, label="t = T")
plt.xlabel("Position in units of L")
plt.ylabel("Probability Density")
plt.title("Probability Density vs. Position", y = 1.05)
plt.grid()
plt.legend()
plt.savefig("PSIsquared.png", dpi = 300, bbox_inches = "tight")


si01 = phi[T_0,:]
si11 = phi[T_1,:]
si21 = phi[T_2,:]
si31 = phi[T_3,:]

psi_square0 = abs(si01)**2
psi_square1 = abs(si11)**2
psi_square2 = abs(si21)**2
psi_square3 = abs(si31)**2

plt.figure(6)
plt.plot(X, psi_square0, label="t = 0")
plt.plot(X, psi_square1, label="t = T/4")
plt.plot(X, psi_square2, label="t = T/2")
plt.plot(X, psi_square3, label="t = T")
plt.xlabel("Position in units of L")
plt.ylabel("Probability Density")
plt.title("Probability Density vs. Position", y = 1.05)
plt.grid()
plt.legend()
plt.savefig("PSIQHO.png", dpi = 300, bbox_inches = "tight")



plt.show()