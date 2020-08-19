import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
L1=8
L2=12
L3=16
Tc=2.269

data1=np.loadtxt('./data/data1.txt')
data2=np.loadtxt('./data/data2.txt')
data3=np.loadtxt('./data/data3.txt')
T_list=data1[:,0]
Q_list_8=data1[:,1]
Q_list_12=data2[:,1]
Q_list_16=data3[:,1]

def generate_xLyL(kappa,nu):
    N=len(T_list)
    xL=np.zeros((3*N,))
    yL=np.zeros((3*N,))
    const=1
    for index in range(len(T_list)):
        t=T_list[index]
        reduce_T=(t-Tc)/Tc
        x1=reduce_T*L1**(1/nu)
        x2=reduce_T*L2**(1/nu)
        x3=reduce_T*L3**(1/nu)
        if abs(x1)<=const:   
            xL[index]=x1
            c1=Q_list_8[index]
            yL[index]=c1*L1**(-kappa/nu)
        if abs(x2)<=const:   
            xL[index+1]=x2
            c2=Q_list_12[index]
            yL[index+1]=c2*L2**(-kappa/nu)
        if abs(x3)<=const:   
            xL[index+2]=x3
            c3=Q_list_16[index]
            yL[index+2]=c3*L3**(-kappa/nu)
    return xL,yL


error_list=[]
exponents=[]
kappa_list=[]
nu_list=[]
for kappa in np.arange(-0.1,0.8,0.01):
    for nu in np.arange(0.8,1.2,0.01):
        tL_list,QL_list=generate_xLyL(kappa,nu)
        tL_array=np.array(tL_list)
        QL_array=np.array(QL_list)
        poly=np.polyfit(tL_array,QL_array,3)
        length=len(QL_list)
        error=0
        for indexnew in range(len(tL_list)):
            x=tL_list[indexnew]
            y1=np.polyval(poly,x)
            y=QL_list[indexnew]
            error += sqrt((y-y1)**2/length)
        exponents.append((kappa,nu))
        kappa_list.append(kappa)
        nu_list.append(nu)
        error_list.append(error)
index_min=error_list.index(min(error_list))
kappa_nu=exponents[index_min]
print(min(error_list))
print(kappa_nu)

#plot 3d
kappa_list_new=[]
nu_list_new=[]
error_list_new=[]
for index in range(len(error_list)):
    error2=error_list[index]
    if error2<0.002:
        kappa_list_new.append(kappa_list[index])
        nu_list_new.append(nu_list[index])
        error_list_new.append(error2)
kappa_array=np.array(kappa_list_new)
nu_array=np.array(nu_list_new)
error_array=np.array(error_list_new)

fig=plt.figure()
ax=Axes3D(fig)
ax.scatter(kappa_array,nu_array,error_array)
plt.show()
