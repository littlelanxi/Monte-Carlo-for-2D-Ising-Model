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
        # print("isweep",i)

def bin_sweep(N_per_bin,beta,mag_sum,mag2_sum,energy1_sum,energy2_sum):
    for i in range(N_per_bin):
        isweep(int(2*N),beta)
        mag_sum += abs(spins_sum)/N  
        mag2_sum += (spins_sum/N)**2
        energy1_sum += current_energy  
        energy2_sum += current_energy**2 
    return mag_sum,mag2_sum,energy1_sum,energy2_sum


if __name__ == '__main__':
    L = 4
    shape = (L, L)
    N = L*L
    moment = 1
    delta=0.05
    # External magnetic field
    field = np.full(shape, 0)

    # Temperature (in units of energy)
    # temperature = 0                   取倒数

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
    M2_list=[]
    M_susceptibility_list=[]
    M2_std=[]
    C_std=[]
    M_suscep_std=[]
    N_per_bin = 4096
    error_n=[]
    N_total_list=[]
    temperature=4.0
    for N_bin in [666,888]:
        N_total=N_bin*N_per_bin
        beta=1/temperature
        # T_list.append(temperature)
        M2_bin=np.zeros((N_bin))
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
            M2_bin[i]=magnetization2_b
            C_bin[i]=heat_capacity_b
            Chi_bin[i]=magnetic_susceptibility_b

        # magnetization2=np.mean(M2_bin)
        stantard_deviation_M2=np.std(M2_bin)/sqrt(N_bin-1)

        N_total_list.append(N_total)
        error_n.append(stantard_deviation_M2*sqrt(N_total))
    
    with open('error-n={}2.txt'.format(N_total),'w') as f:
        f.write(str(N_total_list)+'\n'+str(error_n))

    plt.title('graph of errorbar-N')
    plt.xlabel('N')
    plt.ylabel('errorbar')
    plt.plot(N_total_list,error_n,'bo-')
    # plt.legend()
    plt.savefig('errorbar-N')
    plt.show()
