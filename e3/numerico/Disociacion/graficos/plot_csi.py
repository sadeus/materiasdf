import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("cis.csv")

fsize = 7
plt.figure(figsize=(fsize*(1+np.sqrt(5))/2, fsize))

markers = "o", "x", "d", "p", "+"
labels = "RHF","MP2", "CIS"
max = 0
min = -1
for i in range(2, data.shape[1]):
    plt.plot(data[max:min,1],data[max:min,i],"-", marker=markers[i-2], label = labels[i - 2])

plt.tick_params(labelsize=12)
plt.xlabel("r[A]", fontsize=14)
plt.ylabel("E[ua]", fontsize=14)
plt.grid(axis="y")
plt.legend(loc=0, fontsize=14)
plt.show()
