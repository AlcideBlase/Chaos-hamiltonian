import numpy as np
import math
import matplotlib.pyplot as plt
from playsound import playsound

print('début du programme')
#step definition
h=0.01
#number of points
N=100000

#define constant
for b0 in np.arange(0.29,0.32,0.001): #We define what are the boundary of the thing we will study
    b = b0
    rc = 50
    c = 0.1
    vx = 0
    vy0 = 0.2
    vz0 = 0.2
    vy = vy0
    vz = vz0
    v0=np.sqrt(vx*vx+vy*vy+vz*vz)
    x = 1
    y = 0
    z = 0


#define functions
    def potential(x,y,z):
        return (v0*v0/2)*math.log((x*x+(y/b)*(y/b)+(z/c)*(z/c))/(rc*rc))

    def acceleration_x(x,y,z):
        return -((v0*v0)/(x*x+(y/b)*(y/b)+(z/c)*(z/c)))*x

    def acceleration_y(x,y,z):
        return -((v0*v0)/(x*x+(y/b)*(y/b)+(z/c)*(z/c)))*(y/(b*b))

    def acceleration_z(x,y,z):
        return -((v0*v0)/(x*x+(y/b)*(y/b)+(z/c)*(z/c)))*(z/(c*c))
                    #Runge Kutta order 4 

             #Definition of the localisation of the section
    C_x = 0
    C_y = 0
    C_z = 0

                    #empty list
    X = []
    Y = []
    Z = []
    Vx = []
    Vy = []
    Vz = []
    Coupe_Xy = []
    Coupe_Xz = []
    Coupe_Yx = []
    Coupe_Yz = []
    Coupe_Zx = []
    Coupe_Zy = []

    for k in np.arange (0,N,h):

                        #order1
        kx1 = vx*h
        ky1 = vy*h
        kz1 = vz*h
        ku1 = h*acceleration_x( x, y, z)
        kv1 = h*acceleration_y( x, y, z)
        kw1 = h*acceleration_z( x, y, z)

                        #order2
        kx2 = (vx+0.5*ku1)*h
        ky2 = (vy+0.5*kv1)*h
        kz2 = (vz+0.5*kw1)*h
        ku2 = h*acceleration_x( x+kx1/2, y+ky1/2, z+kz1/2)
        kv2 = h*acceleration_y( x+kx1/2, y+ky1/2, z+kz1/2)
        kw2 = h*acceleration_z( x+kx1/2, y+ky1/2, z+kz1/2)

                        #order3
        kx3 = (vx+0.5*ku2)*h
        ky3 = (vy+0.5*kv2)*h
        kz3 = (vz+0.5*kw2)*h
        ku3 = h*acceleration_x( x+kx2/2, y+ky2/2, z+kz2/2)
        kv3 = h*acceleration_y( x+kx2/2, y+ky2/2, z+kz2/2)
        kw3 = h*acceleration_z( x+kx2/2, y+ky2/2, z+kz2/2)

                        #order4
        kx4 = (vx+ku3)*h
        ky4 = (vy+kv3)*h
        kz4 = (vz+kw3)*h
        ku4 = h*acceleration_x( x+kx3, y+ky3, z+kz3)
        kv4 = h*acceleration_y( x+kx3, y+ky3, z+kz3)
        kw4 = h*acceleration_z( x+kx3, y+ky3, z+kz3)

                        #calcul of the positions
        x += (kx1+2*kx2+2*kx3+kx4)/6.
        y += (ky1+2*ky2+2*ky3+ky4)/6.
        z += (kz1+2*kz2+2*kz3+kz4)/6.

                        #calcul of the celerities
        vx += (ku1+2*ku2+2*ku3+ku4)/6.
        vy += (kv1+2*kv2+2*kv3+kv4)/6.
        vz += (kw1+2*kw2+2*kw3+kw4)/6.

                        #adding to the list to get the informations
        X.append(x)
        Y.append(y)
        Z.append(z)
        Vx.append(vx)
        Vy.append(vy)
        Vz.append(vz)

                        #Poincaré section :
        n = int(k/h)
        if n > 0 :
            verif_X = (X[n-1]-C_x)*(X[n]-C_x)
            verif_Y = (Y[n-1]-C_y)*(Y[n]-C_y)
            verif_Z = (Z[n-1]-C_z)*(Z[n]-C_z)
            if verif_X <=0 :
                Coupe_Xy.append(y)
                Coupe_Xz.append(vy)
            if verif_Y <=0 :
                Coupe_Yz.append(x)
                Coupe_Yx.append(vx)
            if verif_Z <=0 :
                Coupe_Zx.append(x)
                Coupe_Zy.append(vx)

#We create the name of the picture such that we can easly have the information needed when we have a lot of pictures
    name = "b = " + str(round(b,3)) + ", c = " + str(c) + ", rc = " + str(rc) + ", vy0 = " + str(vy0) + ", vz0 = " + str(vz0) + ".png"
    name_save_X = "Images_b_X/Coupe X " + name
    name_save_Y = "Images_b_Y/Coupe Y " + name
    name_save_Z = "Images_b_Z/Coupe Z " + name

#Creation of the pictures
    figX=plt.figure()
    plt.scatter(Coupe_Xy,Coupe_Xz,s=1)
    plt.title(name)
    plt.xlabel("y")
    plt.ylabel("vy")
    plt.savefig(name_save_X, dpi=300, bbox_inches='tight')
    plt.close(figX)

    figY=plt.figure()
    plt.scatter(Coupe_Yz,Coupe_Yx,s=1)
    plt.title(name)
    plt.xlabel("x")
    plt.ylabel("vx")
    plt.savefig(name_save_Y, dpi=300, bbox_inches='tight')
    plt.close(figY)

    figZ=plt.figure()
    plt.scatter(Coupe_Zx,Coupe_Zy,s=1)
    plt.title(name)
    plt.xlabel("x")
    plt.ylabel("vx")
    plt.savefig(name_save_Z, dpi=300, bbox_inches='tight')
    plt.close(figZ)


playsound('pouet.mp3') #The code can be long so it made a sound to let us now when its ending
print('Fin du programme')
