# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 17:15:11 2023

@author: guill
"""

import numpy as np
import matplotlib.pyplot as plt


### constantes definition: ###

#G =
#M =

#constantes of the potential miyamoto-nagai
a = 1
b = 1

#Celerity initial conditions:

vx = 0
vy = 0.2   #modif
vz = 0.2   #modif
#v0=np.sqrt(vx*vx+vy*vy+vz*vz)

#Position initial conditions:

x = 1    #modif
y = 0
z = 0



### Functions definition ###

def potential(x,y,z):
    R=np.sqrt(x*x+y*y)
    return -1/(np.sqrt(R*R+(a+np.sqrt(z*z+b*b))**2))

def acceleration_x(x,y,z):
    return -x/pow(x*x+y*y+(a+np.sqrt(z*z+b*b))**2,1.5)

def acceleration_y(x,y,z):
    return -y/pow(x*x+y*y+(a+np.sqrt(z*z+b*b))**2,1.5)

def acceleration_z(x,y,z):
    return -z*(1+a/np.sqrt(z*z+b*b))/pow(x*x+y*y+(a+np.sqrt(z*z+b*b))**2,3/2)


### Runge Kutta ordre 4 ###

#step definition
h=1e-2

#empty lists
X = []
Y = []
Z = []
Vx = []
Vy = []
Vz = []
E = []

N = 1e3 #number of points

for k in np.arange (0,N,h):

    #ordre1
    kx1 = vx*h
    ky1 = vy*h
    kz1 = vz*h
    ku1 = h*acceleration_x( x, y, z)
    kv1 = h*acceleration_y( x, y, z)
    kw1 = h*acceleration_z( x, y, z)

    #ordre2
    kx2 = (vx+0.5*ku1)*h
    ky2 = (vy+0.5*kv1)*h
    kz2 = (vz+0.5*kw1)*h
    ku2 = h*acceleration_x( x+kx1/2, y+ky1/2, z+kz1/2)
    kv2 = h*acceleration_y( x+kx1/2, y+ky1/2, z+kz1/2)
    kw2 = h*acceleration_z( x+kx1/2, y+ky1/2, z+kz1/2)

    #ordre3
    kx3 = (vx+0.5*ku2)*h
    ky3 = (vy+0.5*kv2)*h
    kz3 = (vz+0.5*kw2)*h
    ku3 = h*acceleration_x( x+kx2/2, y+ky2/2, z+kz2/2)
    kv3 = h*acceleration_y( x+kx2/2, y+ky2/2, z+kz2/2)
    kw3 = h*acceleration_z( x+kx2/2, y+ky2/2, z+kz2/2)

    #ordre4
    kx4 = (vx+ku3)*h
    ky4 = (vy+kv3)*h
    kz4 = (vz+kw3)*h
    ku4 = h*acceleration_x( x+kx3, y+ky3, z+kz3)
    kv4 = h*acceleration_y( x+kx3, y+ky3, z+kz3)
    kw4 = h*acceleration_z( x+kx3, y+ky3, z+kz3)
    
    #computation of the positions 
    x += (kx1+2*kx2+2*kx3+kx4)/6.
    y += (ky1+2*ky2+2*ky3+ky4)/6.
    z += (kz1+2*kz2+2*kz3+kz4)/6.
    
    #computations of the celerities
    vx += (ku1+2*ku2+2*ku3+ku4)/6.
    vy += (kv1+2*kv2+2*kv3+kv4)/6.
    vz += (kw1+2*kw2+2*kw3+kw4)/6.
    
    #adding to the lists
    X.append(x)
    Y.append(y)
    Z.append(z)
    Vx.append(vx)
    Vy.append(vy)
    Vz.append(vz)
    
    #computation of the energy
    E.append(0.5*(vx*vx+vy*vy+vz*vz)+potential(x,y,z))
 
    
 
### Plot 2D trajectoire ###    

plt.figure()
plt.plot(X,Y)
plt.xlabel('x')
plt.ylabel('y')
plt.show()

### Plot 3D trajectoire ###

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(X, Y, Z)
