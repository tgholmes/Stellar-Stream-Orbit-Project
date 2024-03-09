#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 15:44:52 2023

@author: tomholmes
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"

#Keplerian Orbit Integrator: Using 2 Methods (Euler Method and KDK Method)

M = 9 #Units of 10^10 Solar masses
G = 43007.10573 #Units of solar masses, kpc and kms^-1
kpc_to_km = 3.086E+16 #Few times where units need converting from kpc to km
location_0 = np.array([8,0,0]) #Defining x coordinate at t=0 
r_0 = np.sqrt((location_0[0]**2)+(location_0[1]**2)+(location_0[2]**2)) #Defining r
v_0 = np.array([0,np.sqrt((G*M)/r_0),0])
delta_t = 0.001
euler_delta_t = 0.00005 #time in kpc/kms^-1, roughly 1 Gyr
 #Number of iterations the for loops will go through to iterate the Euler and Leapfrog methods
 #Total time that has passed in Gyr

#----------------------------------------------------------------------------------------------------------------

#Defining Euler Step function
def Euler_step(location_0, v_0, delta_t):
    
    a = (-G*M/(np.linalg.norm(location_0)**3))*location_0  
    location_new = np.add(location_0,(v_0*delta_t)) #Creating next position
    v_new = np.add(v_0, (a*delta_t)) #Creating next velocity
        
    return location_new, v_new

#----------------------------------------------------------------------------------------------------------------

a_0 = (-G*M/(np.linalg.norm(location_0)**3))*location_0

#Defining Leapfrog function
def Leapfrog_step(location_0, v_0, a_0, delta_t):
    
    
    v_halfstep = np.add(v_0,(a_0*(delta_t/2)))
    location_new = np.add(location_0,(v_halfstep*delta_t))
    a_new = (-G*M/(np.linalg.norm(location_new)**3))*location_new
    v_new = np.add(v_halfstep,(a_new*(delta_t/2)))
    
    return location_new, v_new, a_new
                        
#----------------------------------------------------------------------------------------------------------------

# Defining Leapfrog integrator with variable timesteps
def Leapfrog_step_variable_t(location_0, v_0, a_0):
    

    eta = 0.1
    vdelta_t = eta*np.sqrt(np.linalg.norm(location_0)/np.linalg.norm(a_0))
    v_halfstep = np.add(v_0,(a_0*(vdelta_t/2)))
    location_new = np.add(location_0,(v_halfstep*vdelta_t))
    a_new = (-G*M/(np.linalg.norm(location_new)**3))*location_new
    v_new = np.add(v_halfstep,(a_new*(vdelta_t/2)))
    
    return location_new, v_new, a_new,eta, vdelta_t

#----------------------------------------------------------------------------------------------------------------

euler_t = 0
euler_t_end = 3
euler_t_array = [0]

#Setting empty lists that can take the outputs of the iterated function
euler_coordinates = [location_0]
euler_velocities = [v_0]

while euler_t < euler_t_end:
    euler_location_new, euler_v_new = Euler_step(location_0,v_0,euler_delta_t)
    euler_coordinates.append(euler_location_new)
    euler_velocities.append(euler_v_new)
    location_0 = euler_location_new
    v_0 = euler_v_new
    euler_t += euler_delta_t
    euler_t_array.append(euler_t)
    
#Resetting these parameters following the Euler function
location_0 = np.array([8,0,0])
v_0 = np.array([0,np.sqrt((G*M)/r_0),0])  

#Converting coordinate and velocity lists to arrays     
np_euler_coordinates = np.array(euler_coordinates)
np_euler_velocities = np.array(euler_velocities)

lf_coordinates_static = [location_0]
lf_velocities_static = [v_0]

lf_coordinates_variable = [location_0]
lf_velocities_variable= [v_0]

# Defining time arrays for Leapfrog method, keeping seperate from Euler time arrays
t_static = 0
t_static_end = 3
t_static_array = [0]

# Iterating static timesteps Leapfrog
while t_static < t_static_end:
    lf_location_new_static, lf_v_new_static, lf_a_new_static = Leapfrog_step(location_0,v_0,a_0,delta_t)
    lf_coordinates_static.append(lf_location_new_static)
    lf_velocities_static.append(lf_v_new_static)
    location_0 = lf_location_new_static
    v_0 = lf_v_new_static
    a_0 = lf_a_new_static
    t_static += delta_t
    t_static_array.append(t_static)

t = 0
t_end = 1
t_array = [0]

# Resetting starting parameters for next integrator
location_0 = np.array([8,0,0])
v_0 = np.array([0,np.sqrt((G*M)/r_0),0])  
apocentre = np.array([4.979649, -14.694955, 0])
# apocentre = np.array([8,0,0])
pericentre = np.array ([-1,0,0])
r_apo = np.sqrt((apocentre[0]**2)+(apocentre[1]**2)+(apocentre[2]**2))
r_peri = np.sqrt((pericentre[0]**2)+(pericentre[1]**2)+(pericentre[2]**2))
ecc_v_0 = np.array(([0,np.sqrt((2*((G*M)/r_peri-(G*M)/r_apo)/((r_apo/r_peri)**2))-1),0]))
# ecc_v_0 = np.array([165.772981, -45.175504, 0])

lf_coordinates_variable = [apocentre]
lf_velocities_variable = [ecc_v_0]

# Iterating Leapfrog variable timesteps method (M=9)
a_0 = (-G*M/(np.linalg.norm(apocentre)**3))*apocentre
while t < t_end:
    lf_location_new, lf_v_new, lf_a_new, eta, vdelta_t = Leapfrog_step_variable_t(location_0,ecc_v_0,a_0)
    lf_coordinates_variable.append(lf_location_new)
    lf_velocities_variable.append(lf_v_new)
    apocentre = lf_location_new
    ecc_v_0 = lf_v_new
    a_0 = lf_a_new
    t += vdelta_t
    t_array.append(t)    

t = 0
t_end = 1
t_array_heavy = [0]

np_lf_coordinates_static = np.array(lf_coordinates_static)
np_lf_velocities_static = np.array(lf_velocities_static)
np_lf_coordinates_variable = np.array(lf_coordinates_variable)
np_lf_velocities_variable = np.array(lf_velocities_variable)

# Heavy M calculations
apocentre = np.array([4.979649, -14.694955, 0])
pericentre = np.array ([-1,0,0])
ecc_v_0 = np.array([165.772981, -45.175504, 0])
lf_coordinates_variable_HEAVY = [apocentre]
lf_velocities_variable_HEAVY = [ecc_v_0]
M=9*1.05
a_0 = (-G*M/(np.linalg.norm(apocentre)**3))*apocentre
while t < t_end:
    lf_location_new_HEAVY, lf_v_new_HEAVY, lf_a_new_HEAVY, eta, vdelta_t = Leapfrog_step_variable_t(apocentre,ecc_v_0,a_0)
    lf_coordinates_variable_HEAVY.append(lf_location_new_HEAVY)
    lf_velocities_variable_HEAVY.append(lf_v_new_HEAVY)
    apocentre = lf_location_new_HEAVY
    ecc_v_0 = lf_v_new_HEAVY
    a_0 = lf_a_new_HEAVY
    t += vdelta_t
    t_array_heavy.append(t)

t = 0
t_end = 1
t_array_light = [0]

np_lf_coordinates_variable_HEAVY = np.array(lf_coordinates_variable_HEAVY)
np_lf_velocities_variable_HEAVY = np.array(lf_velocities_variable_HEAVY)

# LIGHT M CALCULATIONS
apocentre = np.array([4.979649, -14.694955, 0])
pericentre = np.array ([-1,0,0])
ecc_v_0 = np.array([165.772981, -45.175504, 0])
lf_coordinates_variable_LIGHT = [apocentre]
lf_velocities_variable_LIGHT = [ecc_v_0]
M=9*0.97
a_0 = (-G*M/(np.linalg.norm(apocentre)**3))*apocentre
while t < t_end:
    lf_location_new_LIGHT, lf_v_new_LIGHT, lf_a_new_LIGHT, eta, vdelta_t = Leapfrog_step_variable_t(apocentre,ecc_v_0,a_0)
    lf_coordinates_variable_LIGHT.append(lf_location_new_LIGHT)
    lf_velocities_variable_LIGHT.append(lf_v_new_LIGHT)
    apocentre = lf_location_new_LIGHT
    ecc_v_0 = lf_v_new_LIGHT
    a_0 = lf_a_new_LIGHT
    t += vdelta_t
    t_array_light.append(t)

np_lf_coordinates_variable_LIGHT = np.array(lf_coordinates_variable_LIGHT)
np_lf_velocities_variable_LIGHT = np.array(lf_velocities_variable_LIGHT)

t_array_np = np.array(t_array)
t_array_heavy_np = np.array(t_array_heavy)
t_array_light_np = np.array(t_array_light)

  # Plot the Euler orbit
plt.plot(np_euler_coordinates[:,0],np_euler_coordinates[:,1], label = "Euler")
# plt.plot(np_lf_coordinates_static[:,0],np_lf_coordinates_static[:,1], label = "KDK",color='red')
plt.legend(prop={'size': 15},loc='center left', bbox_to_anchor=(1, 0.5))
# plt.title('Orbit Simulation',weight='bold')
plt.xlabel("x (kpc)",weight='bold',fontsize=17)
plt.xticks(fontsize=16)
plt.ylabel("y (kpc)",weight='bold',fontsize=17)
plt.yticks(fontsize=16)
# plt.xlim(-15,15)
# plt.ylim(-15,15)
plt.tight_layout()
plt.gca().set_aspect('equal')
plt.savefig("Euler Orbit Simulation.png",dpi=250)
plt.show()


#Saving values to text files for further analysis if needed
# np.savetxt('Euler Coordinates.txt', euler_coordinates)
# np.savetxt('Euler Velocities.txt', euler_velocities)
# np.savetxt('Leapfrog Coordinates Static.txt', lf_coordinates_static)
# np.savetxt('Leapfrog Velocities Static.txt', lf_velocities_static)

plt.plot(np_lf_coordinates_variable[:,0],np_lf_coordinates_variable[:,1], label = "KDK Simulation (Variable)",color='g')
# plt.legend()
# plt.title(('KDK Variable Orbit with Eta = 0.003'))
plt.xlabel("x (kpc)",weight='bold',fontsize=17)
plt.xticks(fontsize=16)
plt.ylabel("y (kpc)",weight='bold',fontsize=17)
plt.yticks(fontsize=16)
plt.axhline(0, color='black',linewidth=0.75)
plt.axvline(0, color='black',linewidth=0.75)
plt.tight_layout()
plt.gca().set_aspect('equal')
plt.savefig("KDK Variable Orbit Simulation.png",dpi=250)
plt.show()

#Creating arrays that will hold the magnitudes of the velocity vectors and the position vectors
euler_v_mag = np.zeros(len(euler_t_array))
euler_r_mag = np.zeros(len(euler_t_array))

M=9
# EULER CALCULATIONS
#Using for loop to add to the magnitude arrays from the vector arrays
for i in range(len(euler_velocities)):
    
    euler_v_mag[i] = np.linalg.norm(euler_velocities[i])
    euler_r_mag[i] = np.linalg.norm(euler_coordinates[i])
    
#Creating empty array for the euler energy values at each time slot
euler_energy = np.zeros(len(euler_t_array))

#Iterating with for loop to calculate the energy for each v and r magnitudes in the array
for i in range(len(euler_energy)):
    euler_energy[i] = (0.5*euler_v_mag[i]**2)-(G*M/euler_r_mag[i])

# STATIC LEAPFROG CALCULATIONS
#Doing same procedure for Leapfrog method
lf_v_mag_static = np.zeros(len(np_lf_velocities_static))
lf_r_mag_static = np.zeros(len(np_lf_coordinates_static))

for i in range(len(lf_velocities_static)):
    
    lf_v_mag_static[i] = np.linalg.norm(lf_velocities_static[i])
    lf_r_mag_static[i] = np.linalg.norm(lf_coordinates_static[i])

#Creating empty array for the euler energy values at each time slot
lf_energy_static = np.zeros(len(t_static_array))

#Iterating with for loop to calculate the energy for each v and r magnitudes in the array
for i in range(len(lf_energy_static)):
    lf_energy_static[i] = (0.5*lf_v_mag_static[i]**2)-(G*M/lf_r_mag_static[i])
    
# VARIABLE LEAPFROG CALCULATIONS
#Doing same procedure for Leapfrog method
lf_v_mag_variable = np.zeros(len(np_lf_velocities_variable_LIGHT))
lf_r_mag_variable = np.zeros(len(np_lf_coordinates_variable_LIGHT))
                            
for i in range(len(lf_velocities_variable_LIGHT)):
    lf_v_mag_variable[i] = np.linalg.norm(lf_velocities_variable_LIGHT[i])
    lf_r_mag_variable[i] = np.linalg.norm(lf_coordinates_variable_LIGHT[i])
    
#Creating empty array for the euler energy values at each time slot
lf_energy_variable = np.zeros(len(t_array_light_np))


#Iterating with for loop to calculate the energy for each v and r magnitudes in the array
for i in range(len(lf_energy_variable)):
    
    lf_energy_variable[i] = (0.5*lf_v_mag_variable[i]**2)-(G*M/lf_r_mag_variable[i])
    
# Plotting energy against time to see conservation
plt.plot(euler_t_array,euler_energy, label='Euler')
plt.plot(t_static_array, lf_energy_static, label='KDK',color='red')
# plt.title('KDK vs. Euler: Conservation of Energy',weight='bold')
plt.xlabel("Time (~Gyr)",weight='bold',fontsize=17)
plt.xticks(fontsize=16)
plt.ylabel("Energy $\mathregular{(km/s)^{-2}}$",weight='bold',fontsize=17)
plt.yticks(fontsize=16)
plt.ylim(-30000,0)
plt.legend(prop={'size': 15},loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig("Euler vs. KDK Energy.png",dpi=250)
plt.show()

# plt.plot(t_array, lf_r_mag_variable, label='KDK',color='red')
# # plt.title('KDK vs. Euler: Conservation of Energy',weight='bold')
# plt.xlabel("Time (~Gyr)",weight='bold',fontsize=17)
# plt.xticks(fontsize=16)
# plt.ylabel("Radius (kpc)",weight='bold',fontsize=17)
# plt.yticks(fontsize=16)
# # plt.ylim(7.995,8.005)
# # plt.legend()
# plt.tight_layout()
# plt.savefig("KDK radius conservation",dpi=250)
# plt.show()

fig, ax1 = plt.subplots()

ax1.plot(t_array_light_np,lf_energy_variable, label='Energy',color='orange')
# ax1.plot(t_array_np,(((lf_energy_variable-(lf_energy_variable[0]))/np.abs(lf_energy_variable[0]))*100), label='Energy',color='orange')
# ax2.plot(t_array_np,(((lf_r_mag_variable-np.average(lf_r_mag_variable))/np.average(lf_r_mag_variable))*100), label='Radius',color='b')

ax1.set_xlabel('Time (~Gyr)',weight='bold',fontsize=17)
plt.xticks(fontsize=16) #Actual units: kpc/(km/s), but roughly 0.97 Gyr
ax1.set_ylabel('Energy $\mathregular{(km/s)^{-2}}$ ',weight='bold',fontsize=17)
plt.yticks(fontsize=16)
# ax2.set_ylabel('Radius Deviations % change', color='b')
# ax2.set_ylim(-0.01, 0.01)
# ax1.set_ylim(-0.01, 0.01)
# ax1.legend(loc='upper right')
# ax2.legend(loc='upper left')
plt.title('KDK Variable Timesteps: Energy Conservation',weight='bold')
plt.tight_layout()
plt.savefig("Variable Timesteps Conservation.png",dpi=250)
plt.show()

#Calculating percentage error in energy
euler_energy_error = np.abs((np.std(euler_energy)/np.mean(euler_energy))*100)
lf_energy_static_error = np.abs((np.std(lf_energy_static)/np.mean(lf_energy_static))*100)
lf_energy_variable_error = np.abs((np.std(lf_energy_variable)/np.mean(lf_energy_variable))*100)

print("This is the percentage error of the energy conservation from the EULER SIMULATION: ",euler_energy_error, "%")
print("This is the percentage error of the energy conservation from the KDK SIMULATION (WITH CONSTANT TIMESTEPS): ",lf_energy_static_error, "%")
print("AVERAGE ENERGY FOR HEAVY ORBIT: ", np.mean(lf_energy_variable))
print("This is the percentage error of the energy conservation from the KDK SIMULATION (WITH VARIABLE TIMESTEPS): ",lf_energy_variable_error, "%")
print("Value of Eta: ",eta)

# Comparing model orbit to real orbit data
# x_data, y_data, x_vel, y_vel = np.loadtxt("stream_model_GC[11].txt", skiprows=4, usecols = (0,1,3,4),unpack=True)

plt.plot(np_lf_coordinates_variable[:,0],np_lf_coordinates_variable[:,1], label = "KDK (M=9)",color='purple')
# plt.plot(np_lf_coordinates_variable_HEAVY[:,0],np_lf_coordinates_variable_HEAVY[:,1], label = "KDK (M=9.45)",color='grey')
# plt.plot(np_lf_coordinates_variable_LIGHT[:,0],np_lf_coordinates_variable_LIGHT[:,1], label = "KDK (M=8.73)",color='black')
# plt.scatter(x_data,y_data,label='Mock data',color='deepskyblue',s=0.3)
plt.legend(prop={'size': 15},loc='center left', bbox_to_anchor=(1, 0.5))
# plt.title(('KDK model: Comparison with Mock stellar data'))
plt.xlabel("x (kpc)",weight='bold',fontsize=17)
plt.xticks(fontsize=16)
plt.ylabel("y (kpc)",weight='bold',fontsize=17)
plt.yticks(fontsize=16)
plt.tight_layout()
plt.gca().set_aspect('equal')
# plt.axhline(0, color='black',linewidth=0.75)
# plt.axvline(0, color='black',linewidth=0.75)
plt.savefig("KDK Variable vs. Mock Data.png",dpi=250)
plt.show()




