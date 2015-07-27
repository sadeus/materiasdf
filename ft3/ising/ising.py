# -*- coding: utf-8 -*-
"""

"""
import numpy.random as random
import numpy as np
import matplotlib.pyplot as plt

        
def energia(latt):
    resultado = np.roll(latt, -1, axis = 1) 
    resultado += np.roll(latt, 1, axis = 1)
    resultado += np.roll(latt, -1, axis = 0)
    resultado += np.roll(latt, 1, axis = 0)
    resultado *= latt
    print(resultado)
    return np.sum(resultado)

    
def sampleo(latt, L):    
    N = L*L    
    m = 0.0
    e = 0.0
    it = np.nditer(latt, flags=['f_index'])
    while not it.finished:
        m += it[0];
        #Multiplica con los vecinos posteriores.
        i = it.index
        i_der = i + 1
        i_down = i + L
        if (i + 1) % L == 0:
            i_der = i - L + 1        
        if i + L >= N:
            i_down = i + L - N
        e += it[0] * (latt[i_der] + latt[i_down])
        it.iternext()
    m /= N
    e /= N
    return m, e
    
def isingMetropolis(L, T, f_med = 100, n_samp = 1000, n_term = 500):
    random.seed()
    n_iter = n_samp * f_med + n_term
    #Tiene que se grande porque en cada iteración se flipea *un* solo spin
    ### Inicialización de la red    
    N = L*L
    
    latt = np.ones(N)
    ### Temperatura ###
    beta = 1.0/T
    ##########
    ##### Magnetización #####
    mag = np.ndarray(shape = (n_samp))
    e = np.ndarray(shape = (n_samp))
    ##########
    #ims = [latt.copy()] #En ims están las observaciones de la red
    p = [np.exp(-beta * 4), np.exp(-beta * 8)]
    n_samp = 0
    
    def f(x, i):
        '''Calcula la energia asociada a una posición de la red'''
        #print(i)
        dE = latt[(i - 1) % L]
        dE += latt[(i + 1) % L]
        dE += latt[(i - L) % N]
        dE += latt[(i + L) % N]
        dE *= 2 * x
        if dE != 0: 
            if dE < 0 or p[1 if dE == 8 else 0] > np.random.rand():
                x *= 1
    
    for step in range(n_iter):
        #print(latt)
        it = np.nditer(latt, flags=['f_index'], op_flags = [["readwrite"]])
        while not it.finished:
            f(it[0], it.index)
            it.iternext()
        if step % f_med == 0 and step > n_term:
            #ims.append(latt.copy())
            mag[n_samp], e[n_samp] = sampleo(latt, L)
            mag[n_samp] = np.abs(mag[n_samp])
            n_samp += 1
    return T, np.mean(mag), np.var(mag), np.mean(e), np.var(e)
    #plt.savefig('mag_en.png')

if __name__ == "__main__":
    temps = np.concatenate((np.linspace(0.5,1,5), np.arange(1,4,0.1), np.linspace(4, 10, 5)))
    L = 2
    data = []
    for t in temps:
        data.append(isingMetropolis(L,t))
    data = np.array(data)
    plt.figure(1)
    plt.subplot(2,1,1)
    plt.plot(data[:,0],data[:,1],'go-')
    plt.subplot(2,1,2)
    plt.plot(data[:,0],data[:,3],'bo-')
    plt.show()