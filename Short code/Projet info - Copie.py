import numpy as np
import math
import matplotlib.pyplot as plt
from playsound import playsound

print('début du programme')
#step definition
h=0.01
vy0 = 0.2
vz0 = 0.2


b = 0.39
rc = 50
c = 0.1
vx = 0
vy = vy0
vz = vz0
v0=np.sqrt(vx*vx+vy*vy+vz*vz)
x = 1
y = 0
z = 0


def potential(x,y,z):
    return (v0*v0/2)*math.log((x*x+(y/b)*(y/b)+(z/c)*(z/c))/(rc*rc))

def acceleration_x(x,y,z):
    return -((v0*v0)/(x*x+(y/b)*(y/b)+(z/c)*(z/c)))*x

def acceleration_y(x,y,z):
    return -((v0*v0)/(x*x+(y/b)*(y/b)+(z/c)*(z/c)))*(y/(b*b))

def acceleration_z(x,y,z):
    return -((v0*v0)/(x*x+(y/b)*(y/b)+(z/c)*(z/c)))*(z/(c*c))
                #Runge Kutta ordre 4 

         #Definition de l'endroit de la coupe
C_x = 0
C_y = 0
C_z = 0

                #listes vides
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
Coupe_Xvy = []
Coupe_Xvz = []
Coupe_Yvx = []
Coupe_Yvz = []
Coupe_Zvx = []
Coupe_Zvy = []

for k in np.arange (0,100000,h):

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

                    #calcul des positions 
    x += (kx1+2*kx2+2*kx3+kx4)/6.
    y += (ky1+2*ky2+2*ky3+ky4)/6.
    z += (kz1+2*kz2+2*kz3+kz4)/6.

                    #calcul des vitesses
    vx += (ku1+2*ku2+2*ku3+ku4)/6.
    vy += (kv1+2*kv2+2*kv3+kv4)/6.
    vz += (kw1+2*kw2+2*kw3+kw4)/6.

                    #ajout aux listes
    X.append(x)
    Y.append(y)
    Z.append(z)
    Vx.append(vx)
    Vy.append(vy)
    Vz.append(vz)

#Coupe de poincaré :
    n = int(k/h)
    if n > 0 :
        verif_X = (X[n-1]-C_x)*(X[n]-C_x)
        verif_Y = (Y[n-1]-C_y)*(Y[n]-C_y)
        verif_Z = (Z[n-1]-C_z)*(Z[n]-C_z)
        if verif_X <=0 :
            Coupe_Xy.append(y)
            Coupe_Xvy.append(vy)
            Coupe_Xz.append(z)
            Coupe_Xvz.append(vz)
        if verif_Y <=0 :
            Coupe_Yx.append(x)
            Coupe_Yvx.append(vx)
            Coupe_Yz.append(z)
            Coupe_Yvz.append(vz)
        if verif_Z <=0 :
            Coupe_Zx.append(x)
            Coupe_Zvx.append(vx)
            Coupe_Zy.append(y)
            Coupe_Zvy.append(vy)
            
name_X = "b = " + str(round(b,3)) + ", c = " + str(c) + ", rc = " + str(rc) + ", vy0 = " + str(vy0) + ", vz0 = " + str(vz0) + ".png"
name_save_X = "Images_b_X/Coupe X " + name_X
name_Y = "b = " + str(round(b,3)) + ", c = " + str(c) + ", rc = " + str(rc) + ", vy0 = " + str(vy0) + ", vz0 = " + str(vz0) + ".png"
name_save_Y = "Images_b_Y/Coupe Y " + name_Y
name_Z = "b = " + str(round(b,3)) + ", c = " + str(c) + ", rc = " + str(rc) + ", vy0 = " + str(vy0) + ", vz0 = " + str(vz0) + ".png"
name_save_Z = "Images_b_Z/Coupe Z " + name_Z

figX=plt.figure()
plt.scatter(Coupe_Xy,Coupe_Xvy,s=0.05)
plt.title(name_X)
plt.xlabel("y")
plt.ylabel("vy")
plt.savefig("Images_b_X/_y_vy", dpi=300, bbox_inches='tight')
plt.close(figX)

figY=plt.figure()
plt.scatter(Coupe_Xz,Coupe_Xvz,s=0.05)
plt.title(name_Y)
plt.xlabel("z")
plt.ylabel("vz")
plt.savefig("Images_b_X/_z_vz", dpi=300, bbox_inches='tight')
plt.close(figY)

figZ=plt.figure()
plt.scatter(Coupe_Yx,Coupe_Yvx,s=0.05)
plt.title(name_Z)
plt.xlabel("y")
plt.ylabel("vy")
plt.savefig("Images_b_Y/_x_vx", dpi=300, bbox_inches='tight')
plt.close(figZ)

figX=plt.figure()
plt.scatter(Coupe_Yz,Coupe_Yvz,s=0.05)
plt.title(name_X)
plt.xlabel("z")
plt.ylabel("vz")
plt.savefig("Images_b_Y/_z_vz", dpi=300, bbox_inches='tight')
plt.close(figX)

figY=plt.figure()
plt.scatter(Coupe_Zx,Coupe_Zvx,s=0.05)
plt.title(name_Y)
plt.xlabel("x")
plt.ylabel("vx")
plt.savefig("Images_b_Z/_x_vx", dpi=300, bbox_inches='tight')
plt.close(figY)

figZ=plt.figure()
plt.scatter(Coupe_Zy,Coupe_Zvy,s=0.05)
plt.title(name_Z)
plt.xlabel("y")
plt.ylabel("vy")
plt.savefig("Images_b_Z/_y_vy", dpi=300, bbox_inches='tight')
plt.close(figZ)

print('Fin du programme')