import glob
import numpy as np
import matplotlib.pyplot as plt
import fileinput

A = np.linspace(6.8,8.5,100)
for a in A:
	with open('ej1.inp', 'r') as in_file:
		for line in in_file:
			if "celldm(1) =" in line:
				print(line)

	

out_files = glob.glob("./out_data/*")
rexp = "!"
E = []
for g in out_files:
	a = g.split("a_")[1]
	with open(g,"r") as f:
		for line in f:
			if rexp in line:
				E.append([a, line.split("=")[1].strip().split(" ")[0]])

E = np.array(E)
E = E[np.argsort(E[:,0])]

plt.grid()
plt.xlabel("a[ua]")
plt.ylabel("Et[Ry]")
plt.plot(E[:,0],E[:,1],"ro-")
plt.savefig("res.png")
plt.show()

