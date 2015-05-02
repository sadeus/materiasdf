# -*- coding: utf-8 -*-
"""

"""
import numpy.random as random
import numpy as np
import matplotlib.pyplot as plt

    
def deltaE(latt, i):
    '''Calcula la energia asociada a una posición de la red'''
    #def bc(i,N):
    #    '''Determina las condiciones de contorno periodicas'''
        #if i+1 > N - 1:
    #        return 0
    #    elif i-1 < 0:
    #        return N - 1
    #    else:
    #        return i
    L = latt.shape[0]
    temp = latt.ravel()
    N = temp.shape[0]
    result = temp[(i - 1) % L]
    result += temp[(i+1) % L]
    result += temp[(i - L) % N]
    result += temp[(i + L) % N]
    result *= 2 * temp[i] 
    return result
    

def initLatt(L, isRandom = True):
    if isRandom:
        return np.sign(2 * random.rand(L, L) - 1)
    else:
        return np.ones((L,L))
    #return randomLatt(L)
        
def energia(latt):
    resultado = np.roll(latt, -1, axis = 1) 
    resultado += np.roll(latt, 1, axis = 1)
    resultado += np.roll(latt, -1, axis = 0)
    resultado += np.roll(latt, 1, axis = 0)
    resultado *= latt
    print(resultado)
    return np.sum(resultado)

    
def sampleo(latt):    
    temp = latt.ravel()
    L = latt.shape[0]    
    N = temp.shape[0]
    m = 0.0
    e = 0.0
    for i in range(N):
        m += temp[i];
        #Multiplica con los vecinos posteriores.
        i_der = i + 1
        i_down = i + L
        if (i + 1) % L == 0:
            i_der = i - L + 1        
        if i + L >= N:
            i_down = i + L - N
        e += temp[i] * (temp[i_der] + temp[i_down])
    m /= N
    e /= N
    return m, e
    
def isingMetropolis(L, T, f_med = 100, n_samp = 1000, n_term = 500):
    random.seed()
    n_iter = n_samp * f_med + n_term
    #Tiene que se grande porque en cada iteración se flipea *un* solo spin
    ### Inicialización de la red    
    L = 10  #Red de NxN
    N = L*L   
    latt = initLatt(L, isRandom = False)
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
    for step in range(n_iter):
        for i in range(N):
            dE = deltaE(latt, i)
            if dE != 0: #No hacer nada si dE = 0
                if dE < 0 or p[1 if dE == 8 else 0] > np.random.rand():
                    temp = latt.ravel()
                    temp[i] *= -1
                    latt = temp.reshape((L,L))                
        
        if step % f_med == 0 and step > n_term:
            #ims.append(latt.copy())
            mag[n_samp], e[n_samp] = sampleo(latt)
            mag[n_samp] = np.abs(mag[n_samp])
            n_samp += 1
    return T, np.mean(mag), np.var(mag), np.mean(e), np.var(e)
    #plt.savefig('mag_en.png')

if __name__ == "__main__":
    temps = np.concatenate((np.linspace(0.5,1,5), np.arange(1,4,0.1), np.linspace(4, 10, 5)))
    L = 100
    data = []
    for t in temps:
        data.append(isingMetropolis(L,t))
    data = np.array(data)
    plt.figure(1)
    plt.subplot(2,1,1)
    plt.plot(data[:,0],data[:,1],'go-')
    plt.subplot(2,1,2)
    plt.plot(data[:,0],data[:,4]**2/data[:,0]**2,'bo-')
    plt.show()