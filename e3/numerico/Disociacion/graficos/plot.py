import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("plot.csv", delimiter=",",skip_header=True, missing_values="NaN", usecols=(0,1,3, 5, 7)) #3, 5, 7))

fsize = 7
plt.figure(figsize=(fsize*(1+np.sqrt(5))/2, fsize))

markers = "o", "x", "d", "p", "+"

#labels = "RHF(6-31g)", "UHF(6-31g)"
#labels = "B3LYP(cc-pVDZ)", "UB3LYP(cc-pVDZ)"
labels =  "RHF(6-31g)","RHF(cc-pVDZ)","B3LYP(cc-pVDZ)","CISD(cc-pVDZ)"

for i in range(1, data.shape[1]):
    plt.plot(data[:,0],data[:,i],"-", marker=markers[i-1], label = labels[i - 1])

plt.tick_params(labelsize=12)
plt.xlabel("r[A]", fontsize=14)
plt.ylabel("E[ua]", fontsize=14)
plt.grid(axis="y")
plt.legend(loc=0, fontsize=14)

rEq = 1.336
E = -399.1012
#plt.axhline(E, c='k', lw=2, ls="--")
plt.axvline(rEq, c="k", lw=2, ls="--")
#plt.text(3.5, E+0.007,"E$_{S+2H}$ = " + "{} ua".format(E), fontsize=17)
plt.text(rEq + 0.05, max(data[:,1]) - 0.01, r"r$_{eq}$ = " + "{}Å".format(rEq), fontsize=17)

plt.savefig("plot.png", bbox_inches='tight')
