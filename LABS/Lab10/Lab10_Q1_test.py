import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

X = np.arange(0., 2*np.pi, 2*np.pi/5) # from 0 to 2pi (excluded), with 5 points
Y = np.linspace(0., np.pi, 3) # from 0 to pi (included) with 3 points
Yg, Xg = np.meshgrid(Y, X) # 2D grids to create data with correct shape
data = np.sin(Xg*Yg) # creating some data

print(X)
print(Y)

# create nearest interpolator
interp = RegularGridInterpolator((X, Y), data, method='nearest')
x, y = 1., 3. # these coordinates do not fall on a grid point above
nearest_value = interp([x, y])

print(nearest_value)

print(nearest_value, data[1, 2]) # in this example data[1, 2] is the nearest
plt.scatter(Xg, Yg, c=data)
plt.colorbar()
plt.show()
