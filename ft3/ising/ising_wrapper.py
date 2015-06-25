#!/usr/bin/python
# -*- coding: utf-8 -*-
import subprocess as sub
import numpy as np
import os
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (5 * (1 + np.sqrt(5)) / 2, 5)
plt.rcParams["lines.linewidth"] = 2.5
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["axes.labelsize"] = 20

L = 64
path = '.'
temps = np.linspace(0.1,1,20)
fs = 100 #Cantidad de pasos entre mediciones.
n_samp = 2000 #Cantidad de sampleos
warm = 500 #Cantidad de pasos antes de termalizar
seed = 672791038.0
sub.check_output(["gcc","-o","ising","ising.c","-lm"])
out = "med_L_{}".format(L)
filePath = os.path.join(path, out)
mag = []
e = []
c = []
for t in temps:
    cmd = ['./ising']
    cmd += ['-T', str(t)]
    cmd += ['-L', str(L)]
    cmd += ['-n',str(n_samp)]
    cmd += ['-nT',str(warm)]
    cmd += ['-fs',str(fs)]
    cmd += ['-s', str(seed)]
    out = "med_L_{}".format(L)
    #N = warm + fs * Nsamp #Cantidad de iteraciones finales
    filePath = os.path.join(path, out)
    open(filePath,'w+').close() #Lo crea de nuevo el archivo
    with open(filePath, "a+") as file:
        file.write(str(sub.check_output(cmd),"utf-8"))
    data = np.loadtxt(filePath)
    mag.append(np.abs(np.mean(data[:,1])))
    e.append(np.mean(data[:,2]))
    c.append(np.var(data[:,2])**2)

#Tamaño de figuras, siguiendo la regla de oro
plt.figure(1)
plt.xlabel(r"$\beta J$")
plt.ylabel(r'$\langle m \rangle$')           
plt.plot(temps, mag, 'bo-')
plt.grid()
#plt.axvline(np.log(1+ np.sqrt(2))/2, ls = "--", c = "k")
plt.savefig(os.path.join(path, 'mag_L_{}'.format(L)), bbox_inches = 'tight')

plt.figure(2)
plt.xlabel(r"$\beta J$")
plt.ylabel(r'$\langle e \rangle$')   
plt.plot(temps, e,'ro-')
plt.grid()
#plt.axvline(np.log(1+ np.sqrt(2))/2, ls = "--", c = "k")
plt.savefig(os.path.join(path, 'e_L_{}'.format(L)), bbox_inches = 'tight')


plt.figure(3)
plt.xlabel(r"$T' = \frac{k}{J} T$")
plt.ylabel(r"$c = \frac{<\Delta E>^2}{N T'^2}$")  
plt.plot(temps, c,'go-') 
plt.grid()
plt.savefig(os.path.join(path, 'c_L_{}'.format(L)) , bbox_inches = 'tight')

#plt.show()

#
