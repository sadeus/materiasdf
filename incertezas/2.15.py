import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt


def bernoulli(n, p):
    if (p > 1.0): #Probabilidades menores a 1 no existen
        return 0
    count = 0
    for i in range(n):
        r = random.uniform(0, 1)
        if r < p:
            count += 1
    return count
    

def detectaFotones():
    N, n, p = 1000, 15,0.75 #Datos de la experiencia
    bins = np.arange(4, 18, 1)
    A = [] #Array de conteo
    #Experimento
    for i in range(N):
        A.append(bernoulli(n, p))
    #Gráficos
    hist, bins = np.histogram(A, bins = bins)
    histTeo, bins = np.histogram(random.binomial(n, p, N), bins)
    plt.figure()
    N = np.sum(hist)
    yerr = np.sqrt(hist)
    plt.bar(bins[:-1], hist/N, width = 1, yerr=yerr/N, ecolor="r", label="Experiencia")
    plt.plot(bins[:-1] + 0.5, histTeo/np.sum(histTeo), 'go', label="Distribución esperada")
    plt.legend(loc=0)
    plt.savefig("detecta_fotones.png")
    
def fuenteFotones():
    bins = np.arange(0,35,1)
    A = [] #Array de conteo
    N, I, m = 1000, 15, 1000 #Datos de la experiencia
    pEm = I/m  #Probabilidad de emitir un fotón en dt
    #Experimento
    for j in range(N):
        count = 0
        for i in range(m):
            count += bernoulli(1, pEm) #Un proceso de Bernoulli de medir un fotón en un intervalo dt 
        A.append(count) #Después de cada Dt agrego a la lista de experimentos la cuenta de fotones    
    #Histograma
    hist, bins = np.histogram(A, bins = bins)
    #Error de cada bin, con una poissoneana
    yerr = np.sqrt(hist)
    #Grafico
    plt.figure()    
    plt.bar(bins[:-1], hist/np.sum(hist), width = 1, yerr=yerr/np.sum(hist), ecolor="r", label="Experiencia")
    #Gráfico teórico. Poisson con mu = I = 15
    histTeo, bins = np.histogram(random.poisson(lam = I, size = N), bins = bins)
    plt.plot(bins[:-1]+0.5, histTeo/N, 'go', label="Distribución teórica")
    plt.legend(loc=0)
    plt.savefig("fuente_fotones.png")
    
def detectaFuenteFotones():
    bins = np.arange(0,35,1)
    A = []
    p = 0.75
    N = 1000
    I = 15
    m = 1000
    #Experimento
    for j in range(N):
        count = 0
        for i in range(m):
            count += bernoulli(1, I/m) #Generación de 1 fotón por cada dt = Dt/m con p = I/m
        A.append(bernoulli(count,p)) #Dado count fotones, al azar mido con probabilidad p.
    #Histograma
    hist, bins = np.histogram(A, bins = bins)
    #Error de cada bin, con una poissoneana
    N = np.sum(hist)
    yerr = np.sqrt(hist)
    #Figuras    
    plt.figure()    
    plt.bar(bins[:-1], hist/N, yerr=yerr/N, width=1, ecolor="r", label="Experiencia")
    #Distribución teórica
    histTeo, bins = np.histogram(random.poisson(lam = I*p, size = N), bins = bins)
    plt.plot(bins[:-1] + 0.5, histTeo/np.sum(histTeo), 'go', label="Distribución esperada") #Bins corridos para tener el dato teórico al centro
    plt.legend(loc=0)
    plt.savefig("detecta_fuente.png")
    
def probaEfectiva():
    bins = np.arange(0,35,1)
    A = []
    pDec, N, m = 0.75,  1000, 1000
    pEm = 15/m #Probabilidad de fotones por dt
    #Experimento
    for j in range(N):
        count = 0
        for i in range(m):
            count += bernoulli(1,pDec)*bernoulli(1, pEm) #Cuenta de exitos de dos eventos independientes de Bernoulli.
        A.append(count) #La cuenta de exitos en cada iteración la obtengo usando un proceso de Bernoulli
    #Histograma
    hist, bins = np.histogram(A, bins = bins)
    #Error de cada bin, con una multinomial
    N = np.sum(hist)
    yerr = np.sqrt(hist)
    #Figuras
    plt.figure()    
    plt.bar(bins[:-1], hist/N, yerr= yerr/N, ecolor="r", width = 1, label="Experiencia")    
    histTeo, bins = np.histogram(random.poisson(lam = pEm*m*pDec, size = N), bins, density = True)
    plt.plot(bins[:-1], histTeo, 'go', label="Distribución esperada")
    plt.legend(loc=0)
    plt.savefig("proba_efectiva.png")

    
if __name__ == "__main__":
    #Formato de los gráficos consistente
    plt.rcParams["figure.figsize"] = (15/2.54*(np.sqrt(5) + 1)/2,15/2.54)
    plt.rcParams["xtick.labelsize"] = 16
    plt.rcParams["ytick.labelsize"] = 16
    op = input("Ingrese una opción: ")
    if op == "b":
        detectaFotones()
    if op == "c":
         fuenteFotones()
    if op == "d":
        detectaFuenteFotones()
    if op == "e":
        probaEfectiva()
    if op == "all":
        detectaFotones()
        fuenteFotones()
        detectaFuenteFotones()
        probaEfectiva()
    plt.show()
        
    
    
    
    
