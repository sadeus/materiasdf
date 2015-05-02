#!/usr/bin/python
# -*- coding: utf-8 -*-
import subprocess as sub
import numpy as np
import os
import matplotlib.pyplot as plt


def c_wrapper():
    path = '.'
    L = [128]#[8, 16, 64, 128, 256]
    temps = np.linspace(0.1,10,100)
    for l in L:
        fs = 100 #Cantidad de pasos entre mediciones.
        n_samp = 100 #Cantidad de sampleos
        warm = 500 #Cantidad de pasos antes de termalizar
        #N = warm + fs * Nsamp #Cantidad de iteraciones finales
        filePath = os.path.join(path, "med_L_{}".format(l))
        open(filePath,'w+').close() #Lo crea de nuevo el archivo
        for t in temps:
            cmd = ['./ising']
            cmd += ['-T', str(t)]
            cmd += ['-L', str(l)]
            cmd += ['-n',str(n_samp)]
            cmd += ['-nT',str(warm)]
            cmd += ['-fs',str(fs)]
            with open(filePath, "a+") as file:
                file.write(str(sub.check_output(cmd),"utf-8"))
                                
def plotTerm(termPath):
    data = np.loadtxt(termPath, delimiter=',')
    path = os.path.dirname(termPath)
    fileName = os.path.splitext(termPath)[0] + '.png'

    fig = _makeFig(r"Pasos",r"$<m> = \frac{1}{\mu_B} \frac{<M>}{N}$")
    plt.plot(data[:,0],'bo-',figure = fig)
    plt.savefig(os.path.join(path, fileName),
                              bbox_inches = 'tight')
    
def plotData(dataPath):
    data = np.loadtxt(dataPath, delimiter=',')
    path = os.path.dirname(dataPath)
    l = os.path.splitext(dataPath)[0].split("_")[2] #Archivo tiene que ser med_L_{L}.EXT
    
    #Tamaño de figuras, siguiendo la regla de oro

    fig = _makeFig(r"$T' = \frac{k}{J} T$",
                   r'$<m> = \frac{1}{\mu_B} \frac{<M>}{N}$')
    plt.plot(data[:,0],data[:,1],'bo-',figure=fig)
    plt.savefig(os.path.join(path, 'mag_L_{}'.format(l)),
                              bbox_inches = 'tight')
    
    fig = _makeFig(r"$T' = \frac{k}{J} T$",r'$e = \frac{E}{N J}$')
    plt.plot(data[:,0], -data[:,3],'ro-', figure = fig) 
    plt.savefig(os.path.join(path, 'e_L_{}'.format(l)) , bbox_inches = 'tight')
    
    fig = _makeFig(r"$T' = \frac{k}{J} T$",r"$c = \frac{<\Delta E>^2}{N T'^2}$")
    plt.plot(data[:,0], data[:,4]**2/(data[:,0])**2,'go-',figure = fig) 
    plt.savefig(os.path.join(path, 'c_L_{}'.format(l)) , bbox_inches = 'tight')


def _makeFig(xLabel, yLabel):
    fig = plt.figure(figsize =(10, 10 * (np.sqrt(5) - 1)/2))
    plt.tick_params(labelsize= 15)
    plt.ticklabel_format(style = 'sci', scilimits = (-2,2),axis = 'y')
    plt.grid()
    plt.xlabel(xLabel, fontsize=20)
    plt.ylabel(yLabel, fontsize=20)
    return fig

def c_v(): 
    fig = _makeFig(r"$T' = \frac{k}{J} T$",r"$\frac{c}{\max(c)}$")
    plt.plot(T,c_8,'bo-',label='L = 16')
    plt.plot(T,c_16,'r^-',label='L = 64')
    plt.plot(T,c_256,'gs-',label='L = 256')
    plt.legend(loc=0)
