import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define constants
c = 1  ## intercept
a = 2  ## coefficient for percentage of samples
b = 3  ## coefficient for random steps

# Create a grid of x and z values
x = np.linspace(-10, 10, 100)
z = np.linspace(-10, 10, 100)
X, Z = np.meshgrid(x, z)

# Calculate corresponding y values
Y = c + a * X - b * Z

# Plot the surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, alpha=0.7, rstride=100, cstride=100)

# Set labels
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

# Show plot
plt.show()

