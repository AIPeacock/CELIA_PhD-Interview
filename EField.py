import numpy as np
import matplotlib.pyplot as plt

def ESingle(E_0,t,tau,omega):
    E = E_0 * np.exp(-t**2/tau**2)*np.cos(omega*t)
    return (E)

def EDouble(E_0,t,tau,omega):
    E = E_0 * np.exp(-t**2/tau**2) * (np.cos(omega*t) + 0.5 * np.cos(2*omega*t - np.pi/2))
    return (E)
t = np.linspace(-100,100,1000)
# Efield = ESingle(10,t,25,4)
# plt.figure(1)
# plt.plot(t,Efield)
# plt.show()

# Efield2 = EDouble(10,t,25,4)
# plt.figure(2)
# plt.plot(t,Efield2)
# plt.show()
