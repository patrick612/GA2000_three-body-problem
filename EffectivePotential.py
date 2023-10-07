import numpy as np
import matplotlib.pyplot as plt

#Set these
M1 = 10
M2 = 0.5
G = 1
a = 10

R1=M2/(M1+M2)*a
R2=M1/(M1+M2)*a
Omega = (G*(M1+M2)/a**3)**0.5

dims=((-20.0,20.0),(-20.0,20.0))
N=401

xs=np.linspace(*dims[0],N)
ys=np.linspace(*dims[1],N)

x,y = np.meshgrid(xs,ys)

#Position from array indices
pos = lambda x,y:(x*(dims[0][1]-dims[0][0])/xs.shape[0]+dims[0][0],y*(dims[1][1]-dims[1][0])/ys.shape[0]+dims[1][0])

#Potential
U = lambda x,y:-G*M1/((x+R1)**2+y**2)**0.5-G*M2/((x-R2)**2+y**2)**0.5-Omega**2*(x**2+y**2)/2
arr = U(x,y)

#Calculate L1-3
L_points=[]
for i in range(1,N-1):
    if U(xs[i],0)>=U(xs[i-1],0) and U(xs[i],0)>U(xs[i+1],0):
        L_points+=[(i*(dims[0][1]-dims[0][0])/xs.shape[0]+dims[0][0],0)]
#Reorder
L_points = [L_points[i] for i in [1,2,0]]
#Add L4 and L5
L4=pos(*(np.unravel_index(np.argmax(arr),arr.shape)[::-1]))
L_points+=[L4, (L4[0],-L4[1])]

#This is for finding the best possible contour division
L1L3potdiff=U(*L_points[2])-U(*L_points[0])
contourn = 0
eps = 0.001
besteps=1e10
for i in range(1,15):
    L2=(U(*L_points[1])-U(*L_points[0]))/L1L3potdiff
    if abs(L2-int(L2*i)/i)<besteps:
        contourn=i+1
        besteps=abs(L2-int(L2*i)/i)
        if besteps<eps: break

print(f"Lagrange points: {L_points}")

plt.imshow(arr, vmin=-20.0, vmax=-1.0, extent=[*dims[0],*dims[1]])
plt.contour(x,y,arr,np.linspace(U(*L_points[0])-L1L3potdiff*3, U(*L_points[2])+L1L3potdiff*2,contourn*6-5), cmap='binary')
plt.colorbar()
plt.show()
