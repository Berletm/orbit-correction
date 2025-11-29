import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('orbit_data.txt', skiprows=1)
time = data[:, 0] / (24 * 3600)
x = data[:, 1] / 1.496e11
y = data[:, 2] / 1.496e11

plt.figure(figsize=(8, 8))
plt.plot(x, y, 'b-', label='Earth orbit')
plt.plot(0, 0, 'yo', markersize=10, label='Sun')
plt.plot(x[0], y[0], 'go', label='Start')
plt.xlabel('X (AU)')
plt.ylabel('Y (AU)')
plt.title('Earth Orbit Simulation')
plt.axis('equal')
plt.legend()
plt.grid(True)
plt.show()
