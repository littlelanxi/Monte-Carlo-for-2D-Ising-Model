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
Q_list_L1=data1[:,1]
Q_list_L2=data2[:,1]
Q_list_L3=data3[:,1]

def generate_xLyL(kappa,nu):
    N=len(T_list)
    xL=[]
    yL=[]
    limit=1
    for index in range(N):
        T=T_list[index]
        reduced_T=(T-Tc)/Tc
        x1=reduced_T*L1**(1/nu)
        x2=reduced_T*L2**(1/nu)
        x3=reduced_T*L3**(1/nu)
        if abs(x1)<=limit:   
            xL.append(x1)
            c1=Q_list_L1[index]
            yL.append(c1*L1**(-kappa/nu))
        if abs(x2)<=limit:   
            xL.append(x2)
            c2=Q_list_L2[index]
            yL.append(c2*L2**(-kappa/nu))
        if abs(x3)<=limit:   
            xL.append(x3)
            c3=Q_list_L3[index]
            yL.append(c3*L3**(-kappa/nu))
    return xL,yL


error_list=[]
exponents=[]
kappa_list=[]
nu_list=[]
for kappa in np.arange(1,3,0.1):
    for nu in np.arange(0.8,2,0.1):
        tL_list,QL_list=generate_xLyL(kappa,nu)
        tL_array=np.array(tL_list)
        QL_array=np.array(QL_list)
        poly=np.polyfit(tL_array,QL_array,3)
        length=len(QL_array)
        error=0
        for indexnew in range(length):
            x=tL_array[indexnew]
            y1=np.polyval(poly,x)
            y=QL_array[indexnew]
            error += sqrt((y-y1)**2)
        exponents.append((kappa,nu))
        kappa_list.append(kappa)
        nu_list.append(nu)
        error_list.append(error/length)
index_min=error_list.index(min(error_list))
kappa_nu=exponents[index_min]
print(error_list)
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
