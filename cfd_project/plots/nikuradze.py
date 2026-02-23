########################################################################
########################################################################
################### riportare parametri utilizzati #####################
##### turbulent region, k-omeg solution for epsilon/D=0.005 

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

def reynolds(ro,mi,v,Dh):
    
    Re=(ro*v*Dh)/mi
    return Re

def velocity(ro,mi,Re,Dh):

    v=mi*Re/(ro*Dh)
    return v

def p_drop(f,ro,Dh,v,L):

    DeltaP=(f*L*ro*v**2)/(2*Dh)
    return DeltaP

def friction_factors(ro,Dh,L,v,delta_p):

    f = delta_p*2*Dh/(L*ro*v**2)
    return f

Dh=2*0.01875#2*0.0375 #m
ro=997.561 #kg*m^-3
mi= 8.8871E-4 #Pa-s
Re=np.array([4e+3,1.6e+4,2e+4,4.5e+4,6.3e+4,2.5e+5,1e+6])
L=0.2 #m qdelta_L
rel_roughness=5e-3
v=velocity(ro,mi,Re,Dh)
print(v)
p_1,p_2 =  np.loadtxt("/home/f_z/code_project/Starccm/data_nik.txt", unpack=True)

delta_p_sccm= p_1 - p_2 

np.set_printoptions(formatter={'float': '{: 0.2e}'.format})

fig, ax1 = plt.subplots(figsize=(8,6))

### Computing the rough pipe
f_rough=np.zeros(5) #theoretical friction factor

for j in range(len(Re)):
    #f_rough[i,j]=(-1.8*(np.log10(((rel_roughness[i])/(3.7))**(1.11)+ 6.9/Re[j])))**(-2) Haaland
    if j>0:
        f_rough[i,j]=0.25*(np.log10(((rel_roughness)/(3.7))+ 5.74/Re[j]))**(-2) #Jain Swanee >5000
    else:
    #Churcil per ogni Re:
        argomento_log = 1 / ( (7/Re[j])**0.9 + 0.27 * rel_roughness )
        A = (2.457 * np.log(argomento_log))**16
        B = (37530 / Re[j])**16
        term1 = (8 / Re[j])**12
        term2 = 1 / (A + B)**1.5
        f_rough[j] = 8 * (term1 + term2)**(1/12)
    
f_f=friction_factors(ro,Dh,L,v,delta_p_sccm) #sccm friction factor
np.set_printoptions(formatter={'float': '{: 0.2e}'.format})

print(f"f{f_rough}")
print(f"reynolds{Re}")
print("diff",np.abs(f_f-f_rough)/f_rough)
###plotting differences 

ax1.plot(Re, f_f, marker= "s", label="friction factors")
plt.show()

