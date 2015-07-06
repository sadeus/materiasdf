# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 19:31:36 2015

@author: sadeus
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats #Contiene distribuciones útiles
import scipy.optimize as opt

plt.rcParams["figure.figsize"] = (6 * (1 + np.sqrt(5)) / 2, 6)
plt.rcParams["lines.markersize"] = 10
plt.rcParams["lines.linewidth"] = 2.5
plt.rcParams["ytick.labelsize"] = 15
plt.rcParams["xtick.labelsize"] = 15
plt.rcParams["axes.labelsize"] = 22

L = [8, 32, 256]
colors = ["b","g","r","m","c","k"]
for i, l in enumerate(L): 
    datos = np.loadtxt("./2014/med_L_{}.dat".format(l), delimiter=',')
    datos = datos[1/datos[:,0] <= 1]
    k = 1/datos[:,0]
    mag = datos[:,1]
    chi = datos[:,2]
    e = datos[:,3] # Energia
    c = datos[:,4]
    
    plt.figure(5)
    plt.grid()
    plt.xlabel(r"$\beta J$")
    plt.ylabel(r'$m = \frac{M}{N}$')
    plt.plot(k, mag, '{}o-'.format(colors[i]), label="L = {}".format(l))
    plt.legend(loc = 0)
    #plt.axvline(np.log(1+ np.sqrt(2)) / 2, ls = "--", c = "k")
    plt.savefig('mag_L', bbox_inches = 'tight')
    
    #Tamaño de figuras, siguiendo la regla de oro
    plt.figure(1)
    plt.clf()
    plt.grid()
    plt.xlabel(r"$\beta J$")
    plt.ylabel(r'$m = \frac{M}{N}$')
    plt.plot(k, mag,'bo-')
    #plt.axvline(np.log(1+ np.sqrt(2)) / 2, ls = "--", c = "k")
    plt.savefig('mag_L_{}'.format(l), bbox_inches = 'tight')
    
    plt.figure(2)
    plt.clf()
    plt.grid()
    plt.xlabel(r"$\beta J$")
    plt.ylabel(r'$e = \frac{E}{N}$')
    plt.plot(k, e,'ro-')
    #plt.axvline(np.log(1+ np.sqrt(2)) / 2, ls = "--", c = "k")
    plt.savefig('e_L_{}'.format(l), bbox_inches = 'tight')
    
    
    plt.figure(3)
    plt.clf()
    plt.grid()
    plt.xlabel(r"$\beta J $")
    plt.ylabel(r"$\chi = (\Delta m)^2$")
    plt.plot(k, chi,'ko-')
    #plt.axvline(np.log(1+ np.sqrt(2)) / 2, ls = "--", c = "k")
    plt.savefig('chi_L_{}'.format(l), bbox_inches = 'tight')
    
    plt.figure(4)
    plt.clf()
    plt.grid()
    plt.xlabel(r"$\beta J $")
    plt.ylabel(r"$c = (\Delta e)^2$")
    plt.plot(k, c,'go-')
    #plt.axvline(np.log(1+ np.sqrt(2)) / 2, ls = "--", c = "k")
    plt.savefig('c_L_{}'.format(l), bbox_inches = 'tight')