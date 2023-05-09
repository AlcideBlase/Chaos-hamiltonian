import matplotlib.pyplot as plt
import pickle


print('début du programme')


#Definition of the location of the Poincaré section
C_x = 0
C_y = 0
C_z = 0

#List definition

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

#Unpickling the list of data points
data = open('Data/myPoints.pickle', 'rb')
h = pickle.load(data)
X = pickle.load(data)
Y = pickle.load(data)
Z = pickle.load(data)
Vx = pickle.load(data)
Vy = pickle.load(data)
Vz = pickle.load(data)
data.close()

#Unpickling the list of constantes and initial conditions 
data2 = open('Data/constantes.pickle', 'rb')
a = pickle.load(data2)
b = pickle.load(data2)
c = pickle.load(data2)
vx = pickle.load(data2)
vy = pickle.load(data2)
vz = pickle.load(data2)
rc = pickle.load(data2)
data2.close()

#Poincaré section :
for n in range(len(X)):
    if n > 0 :
        verif_X = (X[n-1]-C_x)*(X[n]-C_x)
        verif_Y = (Y[n-1]-C_y)*(Y[n]-C_y)
        verif_Z = (Z[n-1]-C_z)*(Z[n]-C_z)
        if verif_X <=0 :
            Coupe_Xy.append(Y[n])
            Coupe_Xvy.append(Vy[n])
            Coupe_Xz.append(Z[n])
            Coupe_Xvz.append(Vz[n])
        if verif_Y <=0 :
            Coupe_Yx.append(X[n])
            Coupe_Yvx.append(Vx[n])
            Coupe_Yz.append(Z[n])
            Coupe_Yvz.append(Vz[n])
        if verif_Z <=0 :
            Coupe_Zx.append(X[n])
            Coupe_Zvx.append(Vx[n])
            Coupe_Zy.append(Y[n])
            Coupe_Zvy.append(Vy[n])

#We create the name of the picture such that we can easly have the information needed when we have a lot of pictures
name = "a = " + str(a) +", b = " + str(b) + ", c = " + str(c) + ", rc = " + str(rc) + ", vy0 = " + str(vy) + ", vz0 = " + str(vz) 
name_save_X = "Images_b_X/Coupe X " + name
name_save_Y = "Images_b_Y/Coupe Y " + name
name_save_Z = "Images_b_Z/Coupe Z " + name

#Creation of the pictures
figX=plt.figure()
plt.scatter(Coupe_Xy,Coupe_Xvy,s=0.05)
plt.title(name)
plt.xlabel("y")
plt.ylabel("vy")
plt.savefig("Images_b_X/_y_vy", dpi=300, bbox_inches='tight')
plt.close(figX)

figY=plt.figure()
plt.scatter(Coupe_Xz,Coupe_Xvz,s=0.05)
plt.title(name)
plt.xlabel("z")
plt.ylabel("vz")
plt.savefig("Images_b_X/_z_vz", dpi=300, bbox_inches='tight')
plt.close(figY)

figZ=plt.figure()
plt.scatter(Coupe_Yx,Coupe_Yvx,s=0.05)
plt.title(name)
plt.xlabel("x")
plt.ylabel("vx")
plt.savefig("Images_b_Y/_x_vx", dpi=300, bbox_inches='tight')
plt.close(figZ)

figX=plt.figure()
plt.scatter(Coupe_Yz,Coupe_Yvz,s=0.05)
plt.title(name)
plt.xlabel("z")
plt.ylabel("vz")
plt.savefig("Images_b_Y/_z_vz", dpi=300, bbox_inches='tight')
plt.close(figX)

figY=plt.figure()
plt.scatter(Coupe_Zx,Coupe_Zvx,s=0.05)
plt.title(name)
plt.xlabel("x")
plt.ylabel("vx")
plt.savefig("Images_b_Z/_x_vx", dpi=300, bbox_inches='tight')
plt.close(figY)

figZ=plt.figure()
plt.scatter(Coupe_Zy,Coupe_Zvy,s=0.05)
plt.title(name)
plt.xlabel("y")
plt.ylabel("vy")
plt.savefig("Images_b_Z/_y_vy", dpi=300, bbox_inches='tight')
plt.close(figZ)


print('Fin du programme')