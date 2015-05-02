# -*- coding: utf-8 -*-

import numpy as np
import sympy as symb
import matplotlib.pyplot as plt

z, T, V = symb.symbols('z T V')
#para resolver el problema 2 de la guia de Bose, 2C 2014
#defino algunas funciones 

#logaritmo de la funcion Gran Canonica
logZ = lambda z, T, V : symb.log(1-z) + 2*V* T**(3/2) * symb.log(2/(2-z))

#de la derivacion analitica <N> = z d logZGC / dz
Nm = lambda z, T, V : z * symb.diff(logZ(z,T,V), z) 
#Nm = z/(1-z) + 2* V * np.power(T,3/2)* z / (2-z)
#separo explicitamente N0 y Nex

N0= lambda z: z / (1-z)
Nex = lambda z,T,V: 2 * V* T**(3/2) * z /(2-z)

#N = N0 + Nex

#analiticamente encuentro que v T^(3/2) = 1/2 o Tc=1/(2v)^(2/3) me define la
#condicion critica, cuando el numero de particulas es igual al de
#excitados.
Tc = lambda v: (1/(2 * v))**(2/3)

#%La solucion analitica de <N> = N me determina z(N,T,V), esta es la solucion de la cuadratica

z = lambda N, T, V: (2 + 3*N + 2 * T**(3/2) * V - symb.sqrt((N + 2)**2 - 4 * (N-2) * T**(3/2) * V + 4 * T**3 * V**2)) / (2* (1+N+ 2 * T**(3/2) *V))



#Tomo kB=1
#Presion, directamente del Gran Canonico
Pv= lambda z,T,V: T * logZ(z,T,V) / Nm(z,T,V)


#Energia Interna, derivando analiticamente logZGC)
U = lambda z, T, V: 3 * T**(5/2) * V * symb.log(2 / (2 - z))



#EJEMPLOS:
#plots de la fraccion condensada, N0/N

v = 1; # analizo solo este valor de v=V/N
T = np.linspace(0, 3 * Tc(v), 100) #preparo un rango de temperaturas

plt.figure(1)

#clf;
for logN in range(9):
    NN=10**logN;
    plt.subplot(2,1,1)
    plt.plot(T / Tc(v), N0(z(NN,T,v*NN))/ NN,'b-', T/Tc(v), Nex(z(NN, T, v*NN), T, v*NN) / NN, 'r-')
    plt.subplot(2,1,2)
    plt.plot(T / Tc(v),z(NN,T,v*NN))
plt.subplot(2,1,1)
plt.xlabel('T/Tc')
plt.ylabel('Ni/N');
#legend('N0/N','N_{exc}/N')
plt.title('Fraccion condensada, excitaciones y fugacidad para N=1,10,100,..10^8')
plt.subplot(2,1,2)
plt.xlabel('T/Tc')
plt.ylabel('z');


#%NB: aca estoy recalculando en cada plot z, pero deberia
#calcularlo al principio y re-usarlo, para hacerlo mas
#eficiente

#repito para la presion
plt.figure(2)
for logN in range(9):
    NN=10**logN
    plt.plot(T / Tc(v), Pv(z(NN,T,v*NN),T,v*NN))
    
plt.xlabel('T/Tc')
plt.ylabel('Pv');
plt.title('Pv para  N=1,10,100,..10^8')


#Ahora la energia interna por particula
plt.figure(3)
for logN in range(9):
    NN = 10**logN
    plt.plot(T/Tc(v), U(z(NN,T,v*NN),T,v*NN) / NN)
plt.xlabel('T/Tc')
plt.ylabel('U/N')
plt.title('U/N para  N=1,10,100,..10^8')

#Ahora comparo energia interna por particula con 3/2 Pv, son iguales?
plt.figure(4)
for logN in range(0,8,2):
    NN=10**logN
    plt.subplot(4,1,logN/2+1)
    plt.plot(T / Tc(v), U(z(NN,T,v*NN),T,v*NN) / NN,'r-',T / Tc(v), 3/2 * Pv(z(NN,T,v*NN),T,v*NN),'bo')
    plt.title('N = {0}'.format(NN))
#legend('U/N','3/2Pv'
    plt.ylabel('U/N y 3/2 Pv')
plt.xlabel('T/Tc')



#Ahora el Cv (por particula) 
plt.figure(5)
#hold off;
for logN in range(9):
    NN= 10**logN
    #hago el calculo numericamente, la expresion analitica es un poco larga...
    #NB: Cv es a N (no z) constante.
    UiN= U(z(NN,T,v*NN),T,v*NN) / NN
    Cvi = symb.diff(UiN,T)
    plt.plot(T / Tc(v),Cvi)
plt.xlabel('T/Tc')
plt.ylabel('Cv')
plt.title('Cv para N=1,10,100,..10^8')


