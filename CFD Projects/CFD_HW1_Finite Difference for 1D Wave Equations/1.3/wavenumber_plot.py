

import math

import numpy as np


import matplotlib.pyplot as plt
import matplotlib as mpl



def main():

    x_max = 2*math.pi


    t_max = 60

    dx = 2*math.pi/20
    dt = 0.001

    m = 2

    def sin(x):
        return np.sin(x)
    
    def cos(x):
        return np.cos(x)


    x = np.arange(0,math.pi+dx, dx)
    #k = 2*math.pi/dx
    y_exact = x
    
    y_1Re = sin(x)/dx
    y_1Im = (cos(x)-1)/dx

    y_2Re = (4*sin(x)-sin(2*x))/(2*dx)
    y_2Im = (4*cos(x)-cos(2*x)-3)/(2*dx)

    y_3Re = (3*sin(x)-3/2*(sin(2*x))+1/3*sin(3*x))/dx
    y_3Im = (3*cos(x)-3/2*(cos(2*x))+1/3*cos(3*x)-11/6)/dx


    y_3bRe = (6*sin(x)-sin(2*x)+2*sin(x))/(6*dx)
    y_3bIm = (6*cos(x)-cos(2*x)-2*cos(x)-3)/(6*dx)

    
    plt.plot(x, y_exact, color = 'blue', label = 'Exact k')

    plt.plot(x, y_1Re, color = (190/255,184/255,220/255), label = '1st Re k',linewidth = 2)
    plt.plot(x, y_1Im, color = (190/255,184/255,220/255), label = '1st Im k',linewidth = 2)

    plt.plot(x, y_2Re, color = (255/255,190/255,222/255), label = '2st Re k',linewidth = 2)
    plt.plot(x, y_2Im, color = (255/255,190/255,222/255), label = '2st Im k',linewidth = 2)

    plt.plot(x, y_3Re, color = (250/255,127/255,111/255), label = '3st Re k',linewidth = 2)
    plt.plot(x, y_3Im, color = (250/255,127/255,111/255), label = '3st Im k',linewidth = 2)

    plt.plot(x, y_3bRe, color = (130/255,176/255,210/255), label = '3st bias Re k',linewidth = 2)
    plt.plot(x, y_3bIm, color = (130/255,176/255,210/255), label = '3st bias Im k',linewidth = 2)


    #plt.annotate('1st Re k', xy=(3, y_1Re[3]), xytext=(3 +0.1, y_1Re[3]+0.1))

    plt.xlabel('kdx')
    plt.ylabel('Modified k')

    plt.legend()

    plt.show()

    












if __name__ == '__main__':
    main()













