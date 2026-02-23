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
#print(v)
p_1,p_2 =  np.loadtxt("/home/f_z/code_project/Starccm/data_omega.txt", unpack=True)

delta_p_omega= p_1[:7] - p_2[:7]
delta_p_epsilon= p_1[7:] - p_2[7:]


np.set_printoptions(formatter={'float': '{: 0.2e}'.format})

fig, ax1 = plt.subplots(figsize=(8,6))

### Computing the rough pipe
f_rough=np.zeros(7) #theoretical friction factor

for j in range(len(Re)):
    #f_rough[i,j]=(-1.8*(np.log10(((rel_roughness[i])/(3.7))**(1.11)+ 6.9/Re[j])))**(-2) Haaland
    if j!=0:
        f_rough[j]=0.25*(np.log10(((rel_roughness)/(3.7))+ 5.74/(Re[j]**0.9)))**(-2) #Jain Swanee >5000
    else:
    #Churcil per ogni Re:
        argomento_log = 1 / ( (7/Re[j])**0.9 + 0.27 * rel_roughness )
        A = (2.457 * np.log(argomento_log))**16
        B = (37530 / Re[j])**16
        term1 = (8 / Re[j])**12
        term2 = 1 / (A + B)**1.5
        f_rough[j] = 8 * (term1 + term2)**(1/12)
#print("delta_omega",delta_p_omega.shape)
f_omega=friction_factors(ro,Dh,L,v,delta_p_omega) #k-omega friction factor
f_epsilon=friction_factors(ro,Dh,L,v,delta_p_epsilon) #k-omega friction factor
gamma_omega=np.abs(f_omega-f_rough)/f_rough
gamma_epsilon=np.abs(f_epsilon-f_rough)/f_rough
np.set_printoptions(formatter={'float': '{: 0.2e}'.format})

print(f"f{f_rough}")
print(f"reynolds{Re}")
print("diff k-omega",gamma_omega)
print("diff k-eps",gamma_epsilon)
###plotting differences 

ax1.plot(Re, f_omega, marker= "s", label="k-omega")
ax1.plot(Re, f_epsilon, marker= "^", label="k-epsilon")
ax1.plot(Re, f_rough, linestyle= "--", label="Relationships")
ax1.set_xscale("log")
ax1.set_xlabel("Reynolds")
ax1.set_ylabel("f")
plt.grid()
ax1.legend(title="Models, $\epsilon /D=0.005$", title_fontsize='11', 
          fontsize='10', loc='best')

#Error histgrams
fig1, ax = plt.subplots(figsize=(8,6))
bar_width = 0.15
diff_values=np.concatenate((gamma_omega,gamma_epsilon))
indices = np.arange(7) #len(Re)

for i in range(2):
    start = i * 7
    end = (i + 1) * 7   
    
    current_errors = diff_values[start:end]
    if i == 0 :
            model = "omega"  
    else: 
        model= "epsilon"
    offset = (i - 2 // 2) * bar_width
    x_positions = indices + offset

    ax.bar(x_positions, current_errors, 
           width=bar_width,
           edgecolor='black', 
           linewidth=0.5,
           label=f'k-{model}') 
    
ax.set_xlabel('Reynolds', fontsize=12, fontweight='bold')
ax.set_ylabel(' $|f_{CFD} - f_{Corr}|/f_{Corr}$', fontsize=12, fontweight='bold')

unique_reynolds = Re[:7] 
ax.set_xticks(indices)
ax.set_xticklabels([f'{val:.1e}' for val in unique_reynolds])
ax.legend(title="Models", title_fontsize='11', 
          fontsize='10', loc='best')

plt.grid()
plt.show()

