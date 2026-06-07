import numpy as np
import matplotlib.pyplot as plt 

def sawwave(a,n,f,t,tau):
    sawwave = 0 
    for n in range(1,n+1):
        y = -(2*a/np.pi) * ((-1)**n) * (np.sin(n*f*t)/n) * np.exp(-t**2/tau**2)
        sawwave = sawwave + y
    return(sawwave)

# a = 10
# f = 1
# t = np.linspace(-5,5,1000)
# n=6
# plt.figure(n)
# plt.title(f"Waveform")
# plt.plot(t,sawtooth(a,n,f,t))
# plt.show()
