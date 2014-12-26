import numpy as np
import matplotlib.pyplot as plt

# Define fin length and boundary conditions
x1 = 0; alpha = 50
x2 = 0.025; beta = 0

# Define nodal characteristics
npoints = 5
dx = (x2-x1) / npoints

# Define system parameters and geometry
k = 200; h = 100; W = 0.1; t = 0.002
m = np.sqrt((h*2*W)/(k*W*t))
Area = W*t
P = 2*W

# Define "b" matrix - the boundary condition matrix
b = np.zeros((npoints, 1));
b[0] = -alpha
b[-1] = 0

A = np.zeros((npoints, npoints));

# Calculate elements of the "A" matrix, with an adiabatic tip condition
for i in range(npoints):
    for j in range(npoints):
        if j == i: # the diagonal
            A[i,j] = -2-(m*dx)**2
        elif j == i - 1: # left of the diagonal
            A[i,j] = 1
        elif j == i + 1: # right of the diagonal
            A[i,j] = 1
A[npoints-1,npoints-2] = 2 

# solve for the temperature distribution
Y = np.linalg.solve(A,b)

x = np.linspace(x1, x2, npoints + 1)
y = np.hstack([alpha, Y[:,0]])

# numerically solve for the heat flow
Q1 = (y[0]-y[1])
Q2 = y[0]
Qb = (Area*k)*(Q1)*(1/dx)+(P*h*dx*Q2*0.5) 
q = 'The heat flow (approximate) at the base is' 
print q, Qb

# obtain the exact value of heat flow
Qexact = alpha*Area*k*m*np.tanh(m*x2)
qexact = 'The heat flow (exact) at the base is'
print qexact, Qexact

# plot results
plt.plot(x, y, label="Approximate")

z = np.linspace(x1, x2, 100)
plt.plot(z, (alpha*np.cosh(m*(x2-z))) / np.cosh(m*x2) , 'r--', label="Exact")

plt.xlim([0, x2])
plt.xlabel('Distance from base (mm)')
plt.ylabel('Temperature Excess (Celsius)')
legend = plt.legend( loc=1)
plt.show()
