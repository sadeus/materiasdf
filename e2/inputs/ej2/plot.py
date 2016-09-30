import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import integrate

A = np.loadtxt("Alk1579.pdos_tot")

print(sp.integrate.simps(A[:,1],A[:,0]))

plt.plot(A[:,0],A[:,1],"ro-", A[:,0],A[:,2],"go-")

plt.savefig("Alk1579_d_E.png")
plt.show()
