import numpy as np
import matplotlib.pyplot as plt
import EField

def WET(alpha,E,beta):
    prefac = alpha/(np.abs(E))
    exp = np.exp(-beta/(np.abs(E)))
    return (prefac*exp)

def alpha(omega_a, r_h,E_a):
    a = 4 * omega_a * ((r_h)**(5/2)) *E_a
    return(a) 

def beta(r_h,E_a):
    b = 2 * ((r_h)**(3/2)) * (E_a/3)
    return(b)

def Ea(m,q,epsilon_0,hbar):
    num = (m**2) * (q**5)
    denom = (4*np.pi*epsilon_0)**3 * (hbar**4)
    return num/denom

def omega_a(m,q,epsilon_0,hbar):
    num = m * (q**4)
    denom = (4*np.pi*epsilon_0)**2 * (hbar**3)
    return num/denom

def ion_ratio(U_i):
    U_H = 	13.6 #eV should this be in different units?
    return U_i/U_H

t = np.linspace(-100,100,1000)
tau = 20
omega = 1000e-2
E = (EField.ESingle(10e10,t,tau,omega)) #1e10 #V/m
print(len(E))
# plt.figure(1)
# plt.plot(t,E)
# plt.show

m = 9.1093837e-31
q = 1.6e-19
epsilon_0 = 8.85e-12 #F/m
hbar = 1.0545e-34 #J/s

r_h = ion_ratio(13.6)

E_a = Ea(m,q,epsilon_0,hbar)

omegaa = omega_a(m,q,epsilon_0, hbar)

alph = alpha(omegaa,r_h,E_a)
bet = beta(r_h,E_a)


W = WET(alph, E, bet)

# plt.figure(2)
# plt.plot(t,E/np.max(E),color = 'red')
# plt.plot(t,W/np.max(W))
# plt.show()

print("r_h =", r_h)
print("E_a =", E_a)
print("omega_a =", omegaa)
print("alpha =", alph)
print("beta =", bet)
print("W =", np.max(W))