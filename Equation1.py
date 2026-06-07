import numpy as np
import matplotlib.pyplot as plt
import EField
from Equation3 import WET, ion_ratio, Ea, omega_a, alpha, beta
from sawtooth import sawwave
from scipy.fft import fft,fftfreq

m = 9.1093837e-31
q = 1.60217663e-19
epsilon_0 = 8.854e-12 #F/m
hbar = 1.0545718e-34 #J/s

r_h = ion_ratio(15.76)   #Argon ionisation energy 

E_a = Ea(m,q,epsilon_0,hbar) # Calculate E_a used in ionisation rates

omegaa = omega_a(m,q,epsilon_0, hbar) # Calculate Omega_a used in ionisation rate

alph = alpha(omegaa,r_h,E_a) # Calculate Alpha used in Ionisation rate
bet = beta(r_h,E_a) # Calculate Beta used in Ionisation rate

t_max = 100e-15  # Length of simulation
dt = 1e-18  # Time step
num_step = int(2*t_max/dt)
t = np.linspace(-t_max,t_max,num_step)

rho_0 = 1e25    #Atomic density
rho_n = [[0],[0],[0]]  #Set up lists to store free electron density at each time step for each Field
jn = [[0],[0],[0]] #Set up lists to store current density at each time step for each Field

tau = 30e-15   # Envelope duration gives pulse of 10s of femtoseconds
omega = 2*np.pi * 1.875e14# 1600nm Laser
tc = 3000000    #Collision rate or electron decay rate
EStrengthsing = 3.78e10 # Does this give an ionisation yield of 10-20%
EStrengthDoub = 3e10
EStrengthSaw = 4.18e10
EFields= [[0],[0],[0]] #Set up lists to store Electric fields at each time step 
W = [[0],[0],[0]] # Set up lists to store instantaneous ionisation rate at each time step for each field
## Compute individual fields for each time step
EFields[0] = EField.ESingle(EStrengthsing,t,tau,omega)
EFields[1] = EField.EDouble(EStrengthDoub,t,tau,omega)
EFields[2] = sawwave(EStrengthSaw,4,omega,t,tau)

# Main forward Euler compute loop.
for j in range(len(EFields)): # Compute for each Field
    for i in range(num_step-1): # Compute values for each time step (-1 as starts at 0)
        E = EFields[j][i]

        if abs(E) < 1e-30:  # Catch if the electric field is close to 0 as this may induce divison by 0 errors.
            WeT = 0   
        else:
            WeT = WET(alph, E, bet)
        W[j].append(WeT)

        ## Free electron density
        rh = rho_n[j][-1] + dt * WeT * (rho_0 - rho_n[j][-1])
        #print(rh)
        rho_n[j].append(rh)

        #Current density
        jn_1 = jn[j][-1] + dt * ((((q**2)/m) * rho_n[j][-1] * E) - (jn[j][-1] /tc))  
        jn[j].append(jn_1)


#Plotting routines
plt.figure(0)
plt.plot(t,EFields[0],label = "Single Colour")
plt.plot(t,EFields[1],label = "Two Colour")
plt.plot(t,EFields[2],label = "SawTooth")
plt.xlabel("Time (s)")
plt.ylabel("Electric Field (V/m)")
plt.legend()
plt.show()

plt.figure(1)
plt.plot(t,W[0],label = "Single Colour")
plt.plot(t,W[1],label = "Two Colour")
plt.plot(t,W[2],label = "SawTooth")
plt.xlabel("Time (s)")
plt.ylabel("Ionisation rate (arb units)")
plt.legend()
plt.show()

#print(len(rho_n))
plt.figure(2)
plt.plot(t,np.array(rho_n[0])/rho_0,label = "Single Colour")
plt.plot(t,np.array(rho_n[1])/rho_0,label = "Two Colour")
plt.plot(t,np.array(rho_n[2])/rho_0,label = "SawTooth")
plt.xlabel("Time (s)")
plt.ylabel("Free electron density (Normalised against rho_0)")
plt.legend()
plt.show()

plt.figure(3)
plt.plot(t,np.array(jn[0])/np.max(jn[0]),label = "Single Colour")
plt.plot(t,np.array(jn[1])/np.max(jn[1]),label = "Two Colour")
plt.plot(t,np.array(jn[2])/np.max(jn[2]),label = "SawTooth")
plt.xlabel("Time (s)")
plt.ylabel("Current Density (Normalised)")
plt.legend()
plt.show()

# Fourier transform function
def ft(currentdense, time,dt):
    frequencies = fft((currentdense))
    bins = fftfreq(time,dt)
    return (frequencies, bins)

#Compute fourier transform of the current denisites
singfft = ft(jn[0],len(t),dt)
doubfft = ft(jn[1],len(t),dt)
sawfft = ft(jn[2],len(t),dt)

#Plot Fourier transform of current densities to find low frequency current component
plt.figure(4)
plt.plot(doubfft[1],np.abs(doubfft[0])/np.max(np.abs(singfft[0])), label ="Double")
plt.plot(singfft[1],np.abs(singfft[0])/np.max(np.abs(singfft[0])), label ="Single")
plt.plot(sawfft[1],np.abs(sawfft[0])/np.max(np.abs(singfft[0])), label = "Saw")
plt.xlim(0,1e14)
plt.legend()



