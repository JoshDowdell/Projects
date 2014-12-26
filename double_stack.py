import numpy as np
import matplotlib.pyplot as plt

# Define fin length and boundary conditions 
x1 = 0; alpha = 50
x2 = 0.025; 

# Define finite difference grid
npoints = 5
dx = (x2-x1) / npoints

# Input parameters
k = 200; h = 100; W = 0.1; t = 0.002
m = np.sqrt((h*2*W)/(k*W*t))
Area = W*t
P = 2*W

# Calculate constant F2
F2 = 1 / (np.cosh(m*x2)+np.sqrt(2)*np.tanh(m*x2)*np.sinh(m*x2))
	
# Calculate constant K1
K1 = F2*(np.sinh(m*x2)+np.sqrt(2)*np.cosh(m*x2))+2*np.tanh(m*x2)

# Calculate junction temperature 1
F1 = 1 / (np.cosh(m*x2)+K1*np.sinh(m*x2))
theta1 = F1*alpha
J = 'The first junction temperature is' 
print J, theta1

# Calculate junction temperature 2
theta2 = F2*theta1
J = 'The second junction temperature is' 
print J, theta2

# Define first junction boundary condition
beta = theta1

# Define second junction boundary condition
gamma = theta2

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

# Populate required matrices for first set of edge fins
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

# Solve system of equations for finite difference nodes for first set edge fins
Y1 = np.linalg.solve(A1,b1)
z = np.linspace(x1, x2, npoints + 1)
y1 = np.hstack([beta, Y1[:,0]])

############################################################################
# second stack
  
# Populate required matrices for second base fin
b2 = np.zeros((npoints, 1));
b2[0] = -beta
b2[-1] = -gamma

A2 = np.zeros((npoints, npoints));

for i in range(npoints):
    for j in range(npoints):
        if j == i: # the diagonal
            A2[i,j] = -2-(m*dx)**2
        elif j == i - 1: # left of the diagonal
            A2[i,j] = 1
        elif j == i + 1: # right of the diagonal
            A2[i,j] = 1
            
#A[npoints-1,npoints-2] = 2 

# Solve system of equations for finite difference nodes for the base fin
Y2 = np.linalg.solve(A2,b2)
z2 = np.linspace(x1, x2, npoints + 2)
y2 = np.hstack([beta, Y2[:,0], gamma])

# Populate required matrices for second set of edge fins
b3 = np.zeros((npoints, 1));
b3[0] = -gamma
b3[-1] = 0

A3 = np.zeros((npoints, npoints));

for i in range(npoints):
    for j in range(npoints):
        if j == i: # the diagonal
            A3[i,j] = -2-(m*dx)**2
        elif j == i - 1: # left of the diagonal
            A3[i,j] = 1

        elif j == i + 1: # right of the diagonal
            A3[i,j] = 1

A3[npoints-1,npoints-2] = 2 

# Solve system of equations for finite difference nodes for first set edge fins
Y3 = np.linalg.solve(A3,b3)
z3 = np.linspace(x1, x2, npoints + 1)
y3 = np.hstack([gamma, Y3[:,0]])

############################################################################
# heat flows for first stack

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


#####################################################################
# heat flows for second stack

# Use finite differencing to approximate heat flow in the second base fin
Q5 = (y2[0]-y2[1])
Q6 = y2[0]
Qb = (Area*k)*(Q5)*(1/dx)+(P*h*dx*Q6*0.5) 
q = 'The heat flow (approximate) at the second base of the is' 
print q, Qb

Qexact = Area*k*m*(-gamma+np.cosh(m*x2)*beta) / np.sinh(m*x2)
qexact = 'The heat flow (exact) at the second base is'
print qexact, Qexact

Q7 = (y3[0]-y3[1])
Q8 = y3[0]
Qb = (Area*k)*(Q7)*(1/dx)+(P*h*dx*Q6*0.5) 
q = 'The heat flow (approximate) into each second edge fin is' 
print q, Qb

Qexact = gamma*Area*k*m*np.tanh(m*x2)
qexact = 'The heat flow (exact) into each second edge fin is'
print qexact, Qexact

#####################################################################
# plot results 

plt.plot(x, y, label="Approx. (Base 1)")
plt.plot(z, y1, label="Approx. (Edge 1)")

plt.plot(z2, y2, label="Approx. (Base 2)")
plt.plot(z3, y3, label="Approx. (Edge 2)")

zz = np.linspace(x1, x2, 100)

plt.plot(zz, (beta*np.cosh(m*(x2-zz))) / np.cosh(m*x2) , 'b--', label="Exact (Edge 1)")
plt.plot(zz, (beta*(np.sinh(m*zz))+alpha*np.sinh(m*(x2-zz))) / np.sinh(m*x2)  , 'r:', label="Exact (Base 2)")

plt.plot(zz, (gamma*np.cosh(m*(x2-zz))) / np.cosh(m*x2) , 'g-.', label="Exact (Edge 1)")
plt.plot(zz, (gamma*(np.sinh(m*zz))+beta*np.sinh(m*(x2-zz))) / np.sinh(m*x2)  , 'b:', label="Exact (Base 1)")

plt.xlim([0, x2])
plt.xlabel('Distance from base (mm)')
plt.ylabel('Temperature Excess (Celsius)')
plt.legend(loc=1)
plt.show()
