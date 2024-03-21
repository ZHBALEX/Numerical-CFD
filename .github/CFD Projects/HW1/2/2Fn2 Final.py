

import math

import numpy as np


import matplotlib.pyplot as plt
import matplotlib as mpl


from scipy.fftpack import fft, fftfreq

import copy




class UPwind1st_Solver:
    def __init__(self, dx, dt, x_max, t_max,m):
        self.dx = dx
        self.dt = dt
        self.x_max = x_max
        self.t_max = t_max


        self.m = m  #added condition input

    def grid_generate(self):
        self.i_max = int(self.x_max/self.dx)
        self.n_max = int(self.t_max/self.dt)

        self.u = np.zeros(( self.i_max  ))

    
    def IC(self,m):
        for i in range(0, self.i_max):
            self.u[i] = math.sin(m*i*self.dx)

    def Iteration_Formula(self, u_W, u):
        u_new = u + (self.dt/(self.dx))*(u_W-u)
        return u_new

    def Iterative(self):
        u_next = copy.deepcopy(self.u)

        for n in range(1, self.n_max+1):
            for i in range(0,self.i_max):
                u_next[i] = self.Iteration_Formula(self.u[i-1], self.u[i])
            self.u[:] = u_next[:]
        self.u_Full = self.getTHElastBACK(self.u)
        return self.u_Full
    

    def getTHElastBACK(self,u_lost):
        u_Fullget = np.copy(u_lost)
        u_Fullget = np.append(u_Fullget, [u_lost[0]])
        return u_Fullget
    
    def plot_result(self):
        x = np.arange(0,self.i_max+1)
        print(self.u_Full) #TESTING ##############DELETE AFT TEST
        
        plt.plot(x*(2*math.pi/20),self.u_Full)
        plt.show()

    
    def exactsoln(self):
        u_e = np.zeros(self.i_max+1)
        for i in range(self.i_max+1):
            x = i * self.dx
            u_e[i] = math.sin(self.m * (x - self.t_max))
        return u_e

    def plot_result_compare(self):
        x = np.arange(0,self.i_max+1)

        u_e =  self.exactsoln()

        
        plt.plot(x*(2*math.pi/20),self.u_Full)
        plt.plot(x*(2*math.pi/20),u_e)

        print(u_e - self.u_Full)

        plt.show()


    



class Fn1(UPwind1st_Solver):

    def IC(self):
        for i in range(0, self.i_max):
            self.u[i] = math.sin(i*self.dx) +0.5*math.sin(4*i*self.dx)

    def Iteration_Formula(self, u):
        u_new = u - self.dt*(u**2)
        return u_new
    
    def Iterative(self):
        u_next = copy.deepcopy(self.u)

        for n in range(1, self.n_max+1):
            for i in range(0,self.i_max):
                u_next[i] = self.Iteration_Formula(self.u[i])
            self.u[:] = u_next[:]
        self.u_Full = self.getTHElastBACK(self.u)
        return self.u_Full



class Fn2(UPwind1st_Solver):

    def IC(self):
        for i in range(0, self.i_max):
            self.u[i] = math.sin(i*self.dx) +0.5*math.sin(4*i*self.dx)

    def Iteration_Formula(self, u,i):
        u_new = u - self.dt * math.sin(3*i*self.dx)*u
        return u_new

    def Iterative(self):
        u_next = copy.deepcopy(self.u)

        for n in range(1, self.n_max+1):
            for i in range(0,self.i_max):
                u_next[i] = self.Iteration_Formula(self.u[i],i)
            self.u[:] = u_next[:]
        self.u_Full = self.getTHElastBACK(self.u)
        return self.u_Full










class Post_op(): 
    def __init__(self, Function):
        self.Function = Function
    def append(self):
        return 0


    def FFT(self, k, u_analysis, delta_x):
        dx = copy.deepcopy(delta_x)
        u = copy.deepcopy(u_analysis)
        N = self.Function.i_max
        yf = fft(u[:-1])  
        kx = fftfreq(N, dx)[:N // 2]

        y_plot = 2.0 / N * np.abs(yf[:N // 2])

        return kx, y_plot
    




def Aliasing_grid_rearrange(dx, dt):
    t = copy.deepcopy(dt)
    x = copy.deepcopy(dx)
    t = t
    x= x/100
    return x, t



def resultcompare(k, x_max, dx, dxa, y1, y1a):
    i_max = int(x_max/dx)
    i_maxa = int(x_max/dxa)

    x = np.arange(0,i_max+1)

    xa = np.arange(0,i_maxa+1)


    plt.plot(x*dx, y1, label = 'aliasied')  
    plt.plot(xa*dxa, y1a, label = 'unaliasied')  
    
    plt.grid()
    

    plt.xlabel('Result on k=%d'%k)
    plt.ylabel('Amplitude')
    plt.legend()
    plt.show()

    

    

        





def main():

    x_max = 2*math.pi

    dx = 2*math.pi/20
    dt = 0.1


    k = 100                     ## Contorl Iteration counts

    t_max = dt*k

    m = 2


    Function1 = copy.deepcopy(Fn2(dx, dt, x_max, t_max, m))

    Function1.grid_generate()

    Function1.IC()

    u1 = Function1.Iterative()


    Post1 = Post_op(Function1)
    xf, y1  = Post1.FFT(k, u1, dx)


    
    dxa, dta = Aliasing_grid_rearrange(dx, dt)
    



    Function1a = copy.deepcopy(Fn2(dxa, dta, x_max, t_max, m))

    Function1a.grid_generate()

    Function1a.IC()

    u1a = Function1a.Iterative()


    Post1a = Post_op(Function1a)
    xfa, y1a  = Post1a.FFT(k, u1a, dxa)


    #resultcompare(k, x_max, dx, dxa, u1, u1a) # Plot compare in real space

###################### Plot compare in frequence space #############
    plt.plot(xf, y1, label = 'aliasied')  
    plt.plot(xfa, y1a, label = 'unaliasied')  
    
    plt.grid()
    plt.xlim(0,2)

    plt.xlabel('wavenumber')
    plt.ylabel('Amplitude')
    plt.legend()
    plt.title('k = %d'%k)
    plt.show()

######################  ###################### ######################


    




    










if __name__ == '__main__':
    main()













