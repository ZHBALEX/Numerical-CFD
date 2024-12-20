

import math

import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
import copy

class FTCS_Solver:
    def __init__(self, dx, dt, x_max, t_max):
        self.dx = dx
        self.dt = dt
        self.x_max = x_max
        self.t_max = t_max

    def grid_generate(self):
        self.i_max = int(self.x_max/self.dx)
        self.n_max = int(self.t_max/self.dt)

        self.u = np.zeros(( self.i_max +1 ))

    def IC(self):
        for i in range(0, self.i_max+1):
            self.u[i] = 0

    def Iteration_Formula(self, u_W, u_E, u):
        self.Pe = 50
        Pe = self.Pe
        C = self.dt/self.dx
        r = ((1/Pe)*self.dt)/(self.dx**2)

        u_new = (r-C/2)*u_E + (1-2*r)*u + (r+C/2)*u_W

        return u_new

    def Iterative(self):
        u_Next = copy.deepcopy(self.u)
        
        self.Pe = 50

        for n in range(0, self.n_max):
            u_Next[0] = 0
            for i in range(1,self.i_max):
                u_Next[i] = self.Iteration_Formula(self.u[i-1], self.u[i+1], self.u[i])
            u_Next[self.i_max] = 1
            self.u = u_Next



            u_e = copy.deepcopy(self.u)
            for i in range(0,self.i_max+1):
                u_e[i] = u_exact(i*self.dx)
            if np.sum(np.abs(self.u-u_e))/(np.sum(u_e))<0.001:
                print('aaaaaaa')


        return self.u
    
    def plot_result(self):
        x = np.arange(0,self.i_max+1)
        plt.plot(x,self.u)
        plt.show()



    


def Runner(dx, dt, x_max, t_max):

    u_Final = 0

    Run = FTCS_Solver(dx, dt, x_max, t_max)

    Run.grid_generate()

    Run.IC()

    u_Final = Run.Iterative()

    return u_Final



def plot_result(x_max, t_max, dx_20, dx_50, dx_100, u_20, u_50, u_100):
    x_20 = np.linspace(0, x_max, len(u_20))
    x_50 = np.linspace(0, x_max, len(u_50))
    x_100 = np.linspace(0, x_max, len(u_100))

    x_e = np.linspace(0,1,1000)
    u_e = u_exact(x_e)

    plt.figure(figsize=(10, 6))

    plt.plot(x_20, u_20, label='dx=1/20')
    plt.plot(x_50, u_50, label='dx=1/50')
    plt.plot(x_100, u_100, label='dx=1/100')
    plt.plot(x_e, u_e, linestyle='--', label='Exact Solution')

    plt.legend()
    plt.title('Numerical vs. Exact Final Solution ')


    plt.show()



def u_exact(x):
    Pe = 50
    return (np.exp(x*Pe)-1)/(np.exp(Pe)-1)



    


    
    

    


    

class Post_op(FTCS_Solver): #UNFINISH!!!!! DO NOT FORGRT
    

    pass

        





def main():

    x_max = 1

    t_max = 100

    dx = 1/20
    dt = 0.0001

    dx_20 = 1/20
    u_20 = Runner(dx_20, dt, x_max, t_max)

    dx_50 = 1/50
    u_50 = Runner(dx_50, dt, x_max, t_max)

    dx_100 = 1/100
    u_100 = Runner(dx_100, dt, x_max, t_max)

    plot_result(x_max, t_max, dx_20, dx_50, dx_100, u_20, u_50, u_100)





    







if __name__ == '__main__':
    main()













