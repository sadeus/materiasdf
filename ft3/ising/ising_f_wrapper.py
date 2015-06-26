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

L = 128
nterm = 1000 # pasos de termalizacion
ngrupo = 100 # cantidad de subcadenas (muestras)
nfrec = 200 # cada cuantos pasos se toman mediciones dentro de una subcadena
nsize = 100 # cuantas mediciones se hacen dentro de una subcadena, sobre las cuales se promedia
temps = np.concatenate((np.linspace(0.01, 0.43, 10), np.linspace(0.44, 0.46, 10), np.linspace(0.47, 1, 10)))

seed = 672791038.0
#N = warm + fs * Nsamp #Cantidad de iteraciones finales
path = os.path.abspath(".")
dataPath = path + "/data{}".format(L)
if not os.path.isdir(dataPath):
    os.mkdir(dataPath)
os.chdir(dataPath)
np.savetxt("temps", temps)
if os.path.isfile("./Ising"):
    os.remove("./Ising")
sub.call("gfortran {}/ising.f -o Ising".format(path), shell = True)
mag = []
e = []
c = []
chi = []
for i in range(temps.shape[0]):
    filePath = "data{}".format(i)
    open("isingdat.dat", "w+").close()
    s = "{}, {} ! x y dimensiones\n".format(L,L)
    s += "{} , {} , {} ! nterm, ngrupo, nfrec\n".format(nterm, ngrupo, nfrec)
    s += "1, {} ! iflag, nsize\n".format(nsize)
    s += "{} ! acomplamiento\n".format(temps[i])
    s += "0.0 ! b\n"
    s += filePath + "\n"
    s += "{}\n".format(seed)
    with open("isingdat.dat", "a+") as file:
            file.write(s)
    sub.call("./Ising") #Ejecuto
    print("{}% finalizado".format(i/temps.shape[0]*100))
    datos = np.loadtxt(filePath)

    mag_total = datos[:,1]
    mag_cadenas = datos[:,2]
    energia = datos[:,3]
    cal_esp = datos[:,4]

    e.append(np.mean(energia)) # Energia
    mag.append(np.abs(np.mean(mag_cadenas))) # Magnetizacion
    chi.append(np.var(mag_cadenas)) #*(Acopl[i-1]**2) # Suscptibilidad
    c.append(np.mean(cal_esp)) # Calor especifico

mag = np.array(mag)
e = np.array(e)
chi = np.array(chi)
c = np.array(c)

os.chdir(path)
#Tamaño de figuras, siguiendo la regla de oro
plt.figure(1)
plt.grid()
plt.xlabel(r"$\beta J$")
plt.ylabel(r'$m = \frac{M}{N}$')
plt.plot(temps, mag,'bo-')
#plt.axvline(np.log(1+ np.sqrt(2)) / 2, ls = "--", c = "k")
plt.savefig(os.path.join(path, 'mag_L_{}'.format(L)), bbox_inches = 'tight')

plt.figure(2)
plt.grid()
plt.xlabel(r"$\beta J$")
plt.ylabel(r'$e = \frac{E}{N}$')
plt.plot(temps, e,'ro-')
#plt.axvline(np.log(1+ np.sqrt(2)) / 2, ls = "--", c = "k")
plt.savefig(os.path.join(path, 'e_L_{}'.format(L)), bbox_inches = 'tight')


plt.figure(3)
plt.grid()
plt.xlabel(r"$\beta J $")
plt.ylabel(r"$c = (\Delta e)^2$")
plt.plot(temps, c,'go-')
#plt.axvline(np.log(1+ np.sqrt(2)) / 2, ls = "--", c = "k")
plt.savefig(os.path.join(path, 'c_L_{}'.format(L)), bbox_inches = 'tight')

