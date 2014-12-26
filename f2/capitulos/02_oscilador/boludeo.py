import scipy
import numpy
import matplotlib
from scipy.integrate import odeint
from numpy import sin,pi,linspace
from matplotlib.pyplot import plot,show


def f(y,t):
    x = y[0]
    xdot = y[1]
    return [xdot, -sin(x)]

t = linspace(0,10,100)
y0 = [pi-0.01,0]
y = odeint(f, y0, t)

plot(t,y[:,0])
show()
