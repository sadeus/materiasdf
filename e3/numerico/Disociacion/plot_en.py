import glob
import numpy as np
import matplotlib.pyplot as plt


#t = np.arange(2.0, 0.76, -.04)


def read_file(path):
    with open(path, "r") as f:
        s = f.read()
    return s

def parse_variable(s, var):
    init = s.find("|" + var + "=")
    end = s.find("|", init+1)
    return s[init:end]

def parse_energy(s):
    init = s.find("HF=")
    end = s.find("|", init)
    array = s[init+3 : end].strip()
    return np.fromstring(array, dtype=float, sep=",")

files = glob.glob("*/*.dat")
for f in files:
    #path = f.replace("H2S-", "").replace("-SCAN-GEO-SYM.out", "")
    #path = path[path.find("\\") + 1:]
    E = np.loadtxt(f)
    #s = read_file(f).replace("\n\r ", "")
    #E = parse_energy(s)
    #print(E)
    #print(E.size, t.size)
    plt.plot(E[:,1], E[:,2], "o-", label=f)

plt.legend(loc=0)
plt.show()
