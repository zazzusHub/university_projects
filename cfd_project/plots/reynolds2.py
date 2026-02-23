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

'''
Dh=2*0.01875#2*0.0375 #m
ro=997.561 #kg*m^-3
mi= 8.8871E-4 #Pa-s
Re=np.array([4e+3,1.6e+4,6.3e+4,2.5e+5,1e+6])
#v=0.5 #m*s^-1
L=0.2#2.13-0.89
rel_roughness=np.array([1e-4,1e-3,2e-3,5e-3,1e-2])

#Re=reynolds(ro,mi,v,Dh)
#print("Reynolds:",Re)

v=velocity(ro,mi,Re,Dh)
print("Velocity",v)
epsilon=rel_roughness*Dh
print("Rel_rough",epsilon)
f=0.316/(Re**0.25) #lisci blausius<3e+4
#f=0.184/(Re**0.2) #lisci mcadams >3e+4, <1e+6
#f=64/Re

#f =  (-1.8*np.log10(((rel_roughness)/(3.7))**(1.11)+ 6.9/Re))**(-2)
delta_p=p_drop(f,ro,Dh,v,L)
np.set_printoptions(formatter={'float': '{: 0.2e}'.format})
#print("deltaP",np.round(delta_p,2))

a= 6.989193e+00-  5.760494e+00
print(a, 1.5/25.47)

p_1,p_2 =  np.loadtxt("/home/f_z/code_project/Starccm/data.txt", unpack=True)

v=np.array([0.095,0.38,1.5,5.94,23.76]) 
delta_p_sccm= p_1 - p_2 

####Calcolo relativo liscio
f=np.array([0.316/(Re[0]**0.25),0.316/(Re[1]**0.25),0.184/(Re[2]**0.2),0.184/(Re[3]**0.2),0.184/(Re[4]**0.2)]) 
#delta_p=p_drop(f,ro,Dh,v,L)
#diff = delta_p - delta_p_sccm
#rel_error = diff / delta_p_sccm

f_f=friction_factors(ro,Dh,L,v,delta_p_sccm[:5]) 
np.set_printoptions(formatter={'float': '{: 0.2e}'.format})

#print(f"delta_p{delta_p}")
#print(f"delta_p_sccm{delta_p_sccm}")
#print(f"diff{rel_error}")

print(f"f{f}")
print(f"reynolds{Re}")

fig, ax1 = plt.subplots(figsize=(8,6))

ax1.errorbar(Re, f_f, color='black', marker="P", label="friction factors")#np.abs(f-f_f),
# Plot keff values
ax1.set_xlabel("Reynolds")
ax1.set_ylabel("f", color="blue")
ax1.tick_params(axis="y", labelcolor="blue")
ax1.set_xscale("log")
ax1.set_yscale("log")

###color_map
norm = mcolors.LogNorm(vmin=1e-4, vmax=1e-2)
cmap = plt.get_cmap('viridis')
for i in rel_roughness:
    line_color = cmap(norm(i)) #rgba values
    print(i)
    if i==1e-4: 
        epsilon1=line_color
    elif i==1e-3:
        epsilon2=line_color
    elif i==2e-3:
        epsilon3=line_color
    elif i==5e-3:
        epsilon4=line_color
    elif i==1e-2:
        epsilon5=line_color

print("line color",epsilon1,epsilon2,epsilon3)
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([]) 
cbar = plt.colorbar(sm, ax=ax1)
cbar.set_label('Relative Roughness ($\epsilon/D$)')

# Plot relative error on second y-axis

ax2 = ax1.twinx()
tick_vals = [r for r in rel_roughness if r > 0] # Filtriamo eventuali 0
ax2.set_yticks(tick_vals) # Posiziona le tacche
ax2.set_yticklabels([f"{r}" for r in tick_vals]) # Mette le etichette numeriche
ax2.set_ylabel("Relative Roughness $\epsilon/D$", color="gray")
ax2.tick_params(axis="y", labelcolor="gray")

### Computing the rough pipe
f_rough=np.zeros((5,5)) #theoretical friction factor
v=np.array([0.095,0.38,1.5,5.94,23.76,
            0.095,0.38,1.5,5.94,23.76,
            0.095,0.38,1.5,5.94,23.76,
            0.095,0.38,1.5,5.94,23.76,
            0.095,0.38,1.5,5.94,23.76,]) 

for i in range(len(rel_roughness)):
      for j in range(len(Re)):
        f_rough[i,j]=(-1.8*(np.log10(((rel_roughness[i])/(3.7))**(1.11)+ 6.9/Re[j])))**(-2)

f_f=friction_factors(ro,Dh,L,v,delta_p_sccm[5:]) #sccm friction factor
np.set_printoptions(formatter={'float': '{: 0.2e}'.format})

print(f"f{f_rough}")
print(f"reynolds{Re}")
print("diff",np.abs(f_f[5:10]-f_rough[1,:])/f_rough[1,:])
###plotting differences 

#ax1.errorbar(Re, f_f[:5] ,np.abs(f_f[:5]-f_rough[0,:]), color='blue', label="friction factors")
#ax1.errorbar(Re,f_f[5:10] ,   np.abs(f_f[5:10]-f_rough[1,:]), color='red', label="friction factors")
#ax1.errorbar(Re,f_f[10:15] ,  np.abs(f_f[10:15]-f_rough[2,:]), color='purple', label="friction factors")
#ax1.errorbar(Re, f_f[15:20] , np.abs(f_f[15:20]-f_rough[3,:]), color='yellow', label="friction factors")
#ax1.errorbar(Re, f_f[20:] ,   np.abs(f_f[20:]-f_rough[4,:]), color='green', label="friction factors")
ax1.plot(Re, f_f[:5] , marker= "s",color=epsilon1, label="friction factors")
ax1.plot(Re,f_f[5:10] ,   marker= "o",color=epsilon2, label="friction factors")
ax1.plot(Re,f_f[10:15] ,   marker= "^",color=epsilon3, label="friction factors")
ax1.plot(Re, f_f[15:20] ,  marker= "X",color=epsilon4, label="friction factors")
ax1.plot(Re, f_f[20:] ,    marker= "D",color=epsilon5, label="friction factors")

fig1, ax = plt.subplots(figsize=(8,6))

bar_width = 0.15
f_rough_flat = f_rough.flatten() #flatting the 2darray array
diff_values = np.abs(f_f - f_rough_flat)/f_rough_flat
indices = np.arange(5)

for i in range(5):
    # a. Estraiamo i 5 errori corrispondenti alla roughness 'i'
    #    I dati sono ordinati: primi 5 sono Roughness 1, secondi 5 sono Roughness 2...
    start = i * 5
    end = (i + 1) * 5
    
    current_errors = diff_values[start:end]
    
    # b. Calcoliamo il colore
    r_val = rel_roughness[i]
    color = 'black' if r_val <= 0 else cmap(norm(r_val))
    
    # c. Calcoliamo la posizione X delle barre (offset)
    #    (i - 2) serve a centrare il gruppo: -2, -1, 0, +1, +2
    offset = (i - 5 // 2) * bar_width
    x_positions = indices + offset
    
    # d. Plottiamo le barre
    ax.bar(x_positions, current_errors, 
           width=bar_width, 
           color=color, 
           edgecolor='black', 
           linewidth=0.5,
           label=f'$\epsilon/D={r_val:.1e}$') 
    
ax.set_xlabel('Reynolds Re', fontsize=12, fontweight='bold')
ax.set_ylabel(' $|f_{CFD} - f_{Corr}|/f_{Corr}$', fontsize=12, fontweight='bold')

# Impostazione Etichette asse X (i valori di Reynolds)
# Assumiamo che 'Re' contenga i 5 valori unici. Se 'Re' è lungo 25, prendiamo solo i primi 5.
unique_reynolds = Re[:5] 
ax.set_xticks(indices)
ax.set_xticklabels([f'{val:.1e}' for val in unique_reynolds])
ax.legend(title="Relative Roughness", title_fontsize='11', 
          fontsize='10', loc='upper left', bbox_to_anchor=(0.80, 0.90))
plt.grid()
plt.show()

for i in range(len(dp_eval)):
    print(f"{dp_eval[i]:8.2f} | {dp_star[i]:8.2f} | {differenza[i]:8.2f} | {errore_relativo[i]:8.4f}")
###### y+ evaluation
ni=mi/ro#kinematic viscosity 
f= 0.079 * Re**(-0.25)
tau=0.055#=0.5*ro*f*v**2
ustar=(tau/ro)**0.5
delta=5*ni/ustar #sublayer thickenss
ustar= 5*ni/delta
yplus=30 #m
y=yplus*ni/ustar #ydistance from the wall , in this case the centroid of the oitermost cell==> if i need to search for y+=30 then just invert it!

S=5 #mm total thickness
r=1.2 #stretching ratio
n=5


y1=S*(r-1)/(r**n-1)
print(y,y1)

##################
##################
##################
#Re=2.3e+5;Turb_int=0.01;L=2.25 m

delta_p1=2.102494e+04 - 1.963530e+04

#Re=2.3e+5;Turb_int=0.1;L=2.25 m

delta_p2=2.227650e+04 - 2.037559e+04

#Re=2.3e+5;Turb_int=0.01;L=4 m (3. e 3.2)

delta_p3=2.337055e+04 - 2.474837e+04
#                               (1.8 e 2)
delta_p4=2.102662e+04 - 1.963677e+04

print("Re=2.3e+5;Turb_int=0.01;L=2.25 m:", delta_p1)
print("Re=2.3e+5;Turb_int=0.1;L=2.25 m:", delta_p2)
print("Re=2.3e+5;Turb_int=0.01;L=4 m, sezioni 1:", delta_p3)
print("Re=2.3e+5;Turb_int=0.01;L=4 m, sezioni 2:", delta_p4)


#mesh1
delta_p1= 9.700491e-01 - 6.250381e-01
#                               mesh2
delta_p2= 9.774695e-01 - 6.316912e-01
#mesh2 + doubling the number of prism
delta_p3=9.774267e-01 -  6.315597e-01
#mesh1 + doubling the number of prism 
delta_p4=9.705980e-01 - 6.249090e-01
print("mesh1:", delta_p1)
print("mesh2:", delta_p2)
print("mesh2*doubling prismes:", delta_p3)
print("mesh1*doubling prismes:", delta_p3)
'''
    

Re=np.array([1000,1500,2000])
Dh=2*0.01875#2*0.0375 #m
ro=997.561 #kg*m^-3
mi= 8.8871E-4 #Pa*s
L=0.2#2.13-0.89 


p_1,p_2 =  np.loadtxt("/home/f_z/code_project/Starccm/data_lam.txt", unpack=True)
delta_p_CFD=p_1-p_2

v=velocity(ro,mi,Re,Dh)
f=friction_factors(ro,Dh,L,v,delta_p_CFD)

f_f=64/Re
print("velocity",v)
print("f_CFD:",f)
print("f_pouiselle",f_f)
print("Relative:",np.abs(f-f_f)/f_f)



