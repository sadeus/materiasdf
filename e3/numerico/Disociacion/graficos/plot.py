import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("rhf_rohf(triple)_uhf.csv", delimiter=";")

fsize = 7
plt.figure(figsize=(fsize*(1+np.sqrt(5))/2, fsize))

markers = "o", "x", "d", "p", "+"
labels = "ROHF(Triplete)", "RHF (Singlete)", "UHF (Singlete)", "B3LYP (Singlete)", "UB3LYP (Singlete)"
max = 4
min = -3
for i in range(2, data.shape[1]):
    plt.plot(data[max:min,0],data[max:min,i],"-", marker=markers[i-1], label = labels[i - 1])

plt.tick_params(labelsize=12)
plt.xlabel("r[A]", fontsize=14)
plt.ylabel("E[ua]", fontsize=14)
plt.grid(axis="y")
plt.legend(loc=0, fontsize=14)
plt.show()
