

import math

import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl

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

    # def BC(self):

    #     self.u[0] = 0

    #     self.u[self.i_max] = self.u[self.i_max]
        
    
    def IC(self):
        for i in range(0, self.i_max+1):
            self.u[i] = 0

    def Iteration_Formula(self, u_W, u_E, u):
        Pe = 50
        C = self.dt/self.dx
        r = (1/Pe)*self.dt/(self.dx**2)

        u_new = (r-C)*u_E + (1-2*r)*u + (r+C)*u_W
        return u_new

    def Iterative(self):

        for n in range(0, self.n_max):
            self.u[0] = 0
            for i in range(1,self.i_max):
                self.u[i] = self.Iteration_Formula(self.u[i-1], self.u[i+1], self.u[i])
            self.u[self.i_max] = 1

        return self.u
    
    def plot_result(self):
        x = np.arange(0,self.i_max+1)
        plt.plot(x,self.u)
        plt.show()

    


def Runner(dx, dt, x_max, t_max):
    return 0

    


    
    

    


    

class Post_op(): #UNFINISH!!!!! DO NOT FORGRT
    

    pass

        





def main():

    x_max = 1

    t_max = 0.02

    dx = 1/50
    dt = 0.001

    m = 2


    upwind1st = FTCS_Solver(dx, dt, x_max, t_max)

    upwind1st.grid_generate()

    upwind1st.IC()

    a = upwind1st.Iterative()

    print(a)

    upwind1st.Iterative()
    upwind1st.plot_result()

    







if __name__ == '__main__':
    main()













