#!/usr/bin/julia

import PyPlot #ploteo de datos
using LaTeXStrings #Permite usar LaTeX en las cadenas

#Funciones auxiliares. 
#Nota: Los arrays son pasados por referencia -> no es necesario retornarlos

function deltaE(latt, i)
    #Calcula la energia asociada a una posición de la red, 
    #en caso de cambiar el spin i
    L = size(latt)[1]
    N = L * L
    #Condiciones de contorno
    i_down = i % L != 0 ? i + 1 : i - L + 1
    i_up = (i - 1) % L != 0 ? i - 1 : i + L - 1
    i_izq = (i - L) > 0 ? i - L : i - L + N
    i_der = (i + L) <= N ? i + L : i + L - N
    # delta E =  2 * i * (abajo + arriba + izquierda + derecha)
    result = latt[i_up] 
    result += latt[i_down]
    result += latt[i_izq] 
    result += latt[i_der]
    result *= 2 * latt[i] 
    return result
end


function transicion(latt, T, isRandom = false)
    #Lógica para crear configuraciones de spines
    srand(time_ns())
    L = size(latt)[1] #Tamaño de la red
    N = L * L 
    beta = 1.0 / T
    p = [exp(-beta * 4), exp(-beta * 8)]
    for s = 1:N
            #Permite elegir un spin al azar o ordenadamente
            if isRandom
                i = rand(1:L)
            else
                i = s
            end
            dE = deltaE(latt, i)
            if dE != 0 #No hacer nada si dE = 0
                if dE < 0 || p[dE == 4 ? 1 : 2] > rand()
                    latt[i] *= -1                
                end
            end
        end
end

function initLatt(L, isRandom = true)
    #Inicializa la red
     if isRandom
        return sign(rand(Int32, L, L))
        #return rand(-1:2:1, L, L) #Devuelve un número aleatorio de -1:2:1 que es -1,1
    else
        #Devuelve red ordenada con unos.
        return ones(L, L)
    end
end


function sampleo(latt)
    #Funcion utilizada para medir observables
    L = size(latt)[1]
    N = L * L
    m = 0.0
    e = 0.0
    m = sum(latt) #Con sumar es suficiente
    for i = 1:N
        #Multiplica con los vecinos posteriores.
        i_der = (i + L) <= N ? i + L : i + L - N
        i_down = i % L != 0 ? i + 1 : i - L + 1
        e += latt[i] * (latt[i_der] + latt[i_down])
    end
    m /= N
    e /= N
    return (abs(m), e) #Tomo sólo la magnitud de la magnetización
end


function isingMetropolis(L, T, fMed = 100, nSamp = 1000, nTerm = 500, lattRandom = false, randTrans = false) 
    #Función donde se implementa el recorrido
    nIter = nSamp * fMed #Cantidad de iteraciones necesarias para medir lo pedido
    N = L*L   
    latt = initLatt(L, lattRandom) #Inicializa la red
    #Observables
    e = zeros(nSamp + 1)
    mag = zeros(nSamp + 1)
    nSamp = 1 #Resteo el nSamp así lo uso como un indice para los observables
    #Termalizado de la red
    for s = 1:nTerm
        transicion(latt, T, randTrans)
    end
    #Primera medición
    mag[nSamp], e[nSamp] = sampleo(latt)
    nSamp += 1
    for s = 1:nIter
        transicion(latt,T, randTrans)
        if s % fMed == 0
            mag[nSamp], e[nSamp] = sampleo(latt)
            nSamp += 1
        end
    end
    return (T, mag, e)
end


function adqTermalizado(L = 32, T = 2.2, N = 2000, initLattRandom = false, randTrans = false, path = ".")
    #Función que plotea el termalizado de la red
    fMed = 1 #Mide paso a paso de MC
    nTerm = 0 #No termaliza ningún paso
    filePath = "$path/term_L_$L\_$T.csv"
    data = isingMetropolis(L, T, fMed, N, nTerm, initLattRandom, randTrans)
    writedlm(filePath, [data[2] data[3]], ',') #Escribe en el archivo en filePath con coma como delimitador (CSV)
    return filePath
end

function adqObservablesT(L, nTerm = 500, nSamp = 100, fMed = 100, initLattRandom = false, randTrans = false, path = ".")
    temps = [linspace(0.5,1,5),1:0.1:4,linspace(4, 10, 10)] #Temperaturas para graficar observables f(T)
    filePath = path * "/med_L_$L.csv"
    f = open(filePath, "w+")
    for t in temps
        data = isingMetropolis(L, t, fMed, nSamp, nTerm, initLattRandom, randTrans)
        printData = [data[1] mean(data[2]) var(data[2]) mean(data[3]) var(data[3])]
        writedlm(f, printData, ',') #Escribe en el archivo en filePath con coma como delimitador (CSV)
    end
    return filePath
end


function plotStyle(width, xLabel, yLabel)
    #Función para plotear los datos
    figSize = (width, (sqrt(5) - 1)/2 * width) #En pulgadas, sigue la relación aurea
    isGrid, gridTick, gridStyle = true, 1, "--"
    labelSize = 25
    fig = PyPlot.figure(figsize= figSize) 
    PyPlot.tick_params(labelsize = labelSize - 5)
    PyPlot.grid(isGrid, linestyle = gridStyle, linewidth = gridTick)
    PyPlot.xlabel(xLabel, fontsize = labelSize)
    PyPlot.ylabel(yLabel, fontsize = 1.5*labelSize)
    return fig
end
    

#Ejecución del script
L = 64
T = 10
nSamp = 1000
nTerm = 5000
fSamp = 100
initLattRandom = false
randTrans = true

filePath = adqObservablesT(L, nTerm, nSamp, fSamp , initLattRandom, randTrans)
#filePath = adqTermalizado(L, T, nSamp, initLattRandom, randTrans)
figPath = replace(filePath,"csv", "png")
data = readcsv(filePath)

plotStyle(10, "T", L"<m> = \frac{1}{\mu_B} \frac{<M>}{N}")
PyPlot.plot(data[:,1], data[:,2], "bo-", markersize = 8)
PyPlot.savefig(replace(figPath, "med", "mag"), bbox_inches = "tight")

plotStyle(10, "T", L"<e> = \frac{k}{J} \frac{<E>}{N}")
PyPlot.plot(data[:,1], data[:,3], "ro-", markersize = 8)
PyPlot.savefig(replace(figPath, "med", "e"), bbox_inches = "tight")

plotStyle(10, "T", L"c = \frac{<\Delta E>^2}{N T'^2}")
PyPlot.plot(data[:,1], data[:,5].^2 .* data[:,1].^(-2), "go-")
PyPlot.savefig(replace(figPath, "med", "c"), bbox_inches = "tight")

