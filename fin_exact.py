import matplotlib.pyplot as plt
import numpy as np

# Define parameters and geometry 
W=0.1
t=0.002
k=200
h=100
P=2*W
A=W*t
L=0.025

m = np.sqrt((h*P)/(k*A))
x = np.linspace(0, .025, 100)

# Calculate temperature 
y = np.cosh(m*(L-x)) / np.cosh(m*L)

# Plot results
plt.plot(x,y)
plt.xlabel('Distance from base')
plt.ylabel('Temperature excess')
plt.title('Fin Temp exact solution')
plt.show()
