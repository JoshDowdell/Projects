import numpy as np
import matplotlib.pyplot as plt

# Define fin length and boundary conditions 
x1 = 0; alpha = 50
x2 = 0.025; 

# Define finite difference grid
npoints = 100
dx = (x2-x1) / npoints

# Input parameters
k = 200; h = 100; W = 0.1; t = 0.002
m = np.sqrt((h*2*W)/(k*W*t))
Area = W*t
P = 2*W

# Calculate junction temperature
F = 1 / (np.cosh(m*x2)+np.sqrt(2)*np.tanh(m*x2)*np.sinh(m*x2))
theta1 = F*alpha
J = 'The junction temperature is' 
print J, theta1

# Define junction boundary condition
beta = theta1

# Populate required matrices for base fin
b = np.zeros((npoints, 1));
b[0] = -alpha
b[-1] = -beta

A = np.zeros((npoints, npoints));

for i in range(npoints):
    for j in range(npoints):
        if j == i: # the diagonal
            A[i,j] = -2-(m*dx)**2
        elif j == i - 1: # left of the diagonal
            A[i,j] = 1
        elif j == i + 1: # right of the diagonal
            A[i,j] = 1
            
#A[npoints-1,npoints-2] = 2 

# Solve system of equations for finite difference nodes for the base fin
Y = np.linalg.solve(A,b)
x = np.linspace(x1, x2, npoints + 2)
y = np.hstack([alpha, Y[:,0], beta])

# Populate required matrices for edge fins
b1 = np.zeros((npoints, 1));
b1[0] = -beta
b1[-1] = 0

A1 = np.zeros((npoints, npoints));

for i in range(npoints):
    for j in range(npoints):
        if j == i: # the diagonal
            A1[i,j] = -2-(m*dx)**2
        elif j == i - 1: # left of the diagonal
            A1[i,j] = 1
        elif j == i + 1: # right of the diagonal
            A1[i,j] = 1

A1[npoints-1,npoints-2] = 2 

# Solve system of equations for finite difference nodes for edge fins
Y1 = np.linalg.solve(A1,b1)
z = np.linspace(x1, x2, npoints + 1)
y1 = np.hstack([beta, Y1[:,0]])

# Use finite differencing to approximate heat flow in the base fin
Q1 = (y[0]-y[1])
Q2 = y[0]
Qb = (Area*k)*(Q1)*(1/dx)+(P*h*dx*Q2*0.5) 
q = 'The heat flow (approximate) at the base is' 
print q, Qb

Qexact = Area*k*m*(-beta+np.cosh(m*x2)*alpha) / np.sinh(m*x2)
qexact = 'The heat flow (exact) at the base is'
print qexact, Qexact

Q3 = (y1[0]-y1[1])
Q4 = y1[0]
Qb = (Area*k)*(Q3)*(1/dx)+(P*h*dx*Q4*0.5) 
q = 'The heat flow (approximate) into each edge fin is' 
print q, Qb

Qexact = beta*Area*k*m*np.tanh(m*x2)
qexact = 'The heat flow (exact) into each edge fin is'
print qexact, Qexact

# Plot the results
plt.plot(x, y, label="Approximate (Base Fin)")
plt.plot(z, y1, label="Approximate (Edge Fins)")

zz = np.linspace(x1, x2, 100)
plt.plot(zz, (beta*np.cosh(m*(x2-zz))) / np.cosh(m*x2) , 'b--', label="Exact (Edge Fins)")
plt.plot(zz, (beta*(np.sinh(m*zz))+alpha*np.sinh(m*(x2-zz))) / np.sinh(m*x2)  , 'r--', label="Exact (Base Fin)")

plt.xlim([0, x2])
plt.xlabel('Distance from base (mm)')
plt.ylabel('Temperature Excess (Celsius)')
legend = plt.legend( loc=1)
plt.show()
