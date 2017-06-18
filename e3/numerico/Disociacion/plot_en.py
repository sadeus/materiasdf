import glob
import numpy as np
import matplotlib.pyplot as plt


t = np.arange(2.0, 0.76, -.04)


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

files = glob.glob("*/*.out")
for f in files:
    path = f.replace("H2S-", "").replace("-SCAN-GEO-SYM.out", "")
    path = path[path.find("\\") + 1:]

    s = read_file(f).replace("\n ", "")
    E = parse_energy(s)
    if E.size > 0:
        plt.plot(t, E, "o-", label=path)

plt.legend(loc=0)
plt.show()
#with open(name,"r") as f:
#    s = ""
#    for i in f:
#        s += i
#    fInd = s.find("HF=")
#    lInd = s.find("|", fInd)
#    A = np.array(s[fInd+3:lInd].replace("\n ", "").split(","), dtype=float)
#    plt.plot(t, A, 'go-')
#    plt.show()
