import numpy as np
import matplotlib.pyplot as plt
from math import sqrt


def get_probability(delta_energy, beta):
    return np.exp(-delta_energy*beta)

def get_energy():
    return -np.sum(
    interaction * spins * np.roll(spins, 1, axis=0) +
    interaction * spins * np.roll(spins, -1, axis=0) +
    interaction * spins * np.roll(spins, 1, axis=1) +
    interaction * spins * np.roll(spins, -1, axis=1)
  )/2 - moment * np.sum(field * spins)


def update(beta):
    global spins
    global current_energy
    global spins_sum
    i = np.random.randint(spins.shape[0])
    j = np.random.randint(spins.shape[1])
    delta_energy=2*spins[i,j]*(spins[(i+1)%L,j]+spins[(i-1)%L,j]+spins[i,(j+1)%L]+spins[i,(j-1)%L])
    if get_probability(delta_energy, beta) > np.random.random():
        spins[i,j] *= -1
        current_energy += delta_energy
        spins_sum += spins[i,j]*2


def isweep(N,beta):
    for i in range(N):
        update(beta)

def bin_sweep(N_per_bin,beta,mag_sum,mag2_sum,energy1_sum,energy2_sum):
    for i in range(N_per_bin):
        isweep(N,beta)
        mag_sum += abs(spins_sum)/N  
        mag2_sum += (spins_sum/N)**2
        energy1_sum += current_energy  
        energy2_sum += current_energy**2 
    return mag_sum,mag2_sum,energy1_sum,energy2_sum


if __name__ == '__main__':
    L = 8
    shape = (L, L)
    N = L*L
    moment = 1
    delta=0.05
    # External magnetic field
    field = np.full(shape, 0)

    # Interaction (ferromagnetic if positive, antiferromagnetic if negative)
    interaction = 1

    # Spin configuration
    global spins
    spins = np.full(shape,1)
    global spins_sum
    global current_energy
    current_energy=get_energy() 
    T_list=[]
    C_list=[]
    m2_list=[]
    M_susceptibility_list=[]
    m2_std=[]
    C_std=[]
    M_suscep_std=[]
    N_bin=50
    N_per_bin=5000
    N_total=N_bin*N_per_bin
    for temperature in np.arange(1.5,4,delta):
        beta=1/temperature
        T_list.append(temperature)
        m2_bin=np.zeros((N_bin))
        C_bin=np.zeros((N_bin))
        Chi_bin=np.zeros((N_bin))
        spins_sum=np.sum(spins)
        isweep(int(N_total/10),beta)
        for i in range(N_bin):
            magnetization_sum=0
            magnetization2_sum=0
            energy_sum=0
            energy_2_sum=0
            magnetization_sum,magnetization2_sum,energy_sum,energy_2_sum=bin_sweep(N_per_bin,beta,magnetization_sum,magnetization2_sum,energy_sum,energy_2_sum)
            magnetization_b=magnetization_sum/N_per_bin
            magnetization2_b=magnetization2_sum/N_per_bin
            energy_b=energy_sum/N_per_bin
            energy2_b=energy_2_sum/N_per_bin
            heat_capacity_b=(energy2_b-energy_b**2)*(beta**2)/N
            magnetic_susceptibility_b=(magnetization2_b-magnetization_b**2)*N*beta
            m2_bin[i]=magnetization2_b
            C_bin[i]=heat_capacity_b
            Chi_bin[i]=magnetic_susceptibility_b

        magnetization2=np.mean(m2_bin)
        stantard_deviation_m2=np.std(m2_bin)/sqrt(N_bin-1)

        heat_capacity=np.mean(C_bin)
        stantard_deviation_C=np.std(C_bin)/sqrt(N_bin-1)

        magnetic_susceptibility=np.mean(Chi_bin)
        stantard_deviation_Msus=np.std(Chi_bin)/sqrt(N_bin-1)
        
        m2_list.append(magnetization2)
        m2_std.append(stantard_deviation_m2)
        C_list.append(heat_capacity)
        C_std.append(stantard_deviation_C)
        M_susceptibility_list.append(magnetic_susceptibility)
        M_suscep_std.append(stantard_deviation_Msus)

    T_m2_array[:,0]=np.array(T_list)
    T_m2_array[:,1]=np.array(m2_list)
    T_C_array[:,0]=np.array(T_list)
    T_C_array[:,1]=np.array(C_list)
    T_Chi_array[:,0]=np.array(T_list)
    T_Chi_array[:,1]=np.array(M_susceptibility_list)
    f1='data1_L={}.txt'.format
    np.savetxt(f1,T_m2_array)
    f2='data2_L={}.txt'.format
    np.savetxt(f2,T_C_array)
    f3='data3_L={}.txt'.format
    np.savetxt(f3,T_Chi_array)
    
    plt.title('graph of m^2-T')
    plt.errorbar(T_list,m2_list,yerr=m2_std,xerr=None,fmt='.-',ecolor='r',elinewidth=1,capsize=2)
    plt.xlabel('T/t0')
    plt.ylabel('m^2')
    # plt.plot(T_list,m2_list,'bo-',label='N={}*{} M^2--T'.format(n,n))
    # plt.legend()
    plt.savefig('L={} m^2-T_new1'.format(L))
    plt.show()

    plt.title('graph of C-T')
    plt.xlabel('T/t0')
    plt.ylabel('C')
    plt.errorbar(T_list,C_list,yerr=C_std,fmt='.-',ecolor='r',elinewidth=1,capsize=2)
    # plt.legend()
    plt.savefig('L={} C-T_new1'.format(L))
    plt.show()

    plt.title('graph of magnetic_susceptibility-T')
    plt.xlabel('T/T0')
    plt.ylabel('magnetic_susceptibility')
    plt.errorbar(T_list,M_susceptibility_list,yerr=M_suscep_std,fmt='.-',ecolor='r',elinewidth=1,capsize=2)
    # plt.plot(T_list,M_susceptibility_list,'bo-',label='N={}*{} M_S--T'.format(n,n))
    # plt.legend()
    plt.savefig('L={} M_S^2-T_new1'.format(L))
    plt.show()
