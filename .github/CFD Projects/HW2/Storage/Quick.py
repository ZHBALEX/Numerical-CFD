

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

        u_new = (r-C)*u_E + (1-2*r)*u + (r+C)*u_W
        return u_new

    def Iterative(self):
        u_Next = copy.deepcopy(self.u)

        for n in range(0, self.n_max):
            u_Next[0] = 0
            for i in range(1,self.i_max):
                u_Next[i] = self.Iteration_Formula(self.u[i-1], self.u[i+1], self.u[i])
            u_Next[self.i_max] = 1
            self.u = u_Next

        return self.u
    
    def plot_result(self):
        x = np.arange(0,self.i_max+1)
        plt.plot(x,self.u)
        plt.show()


###############################################################################################
###############################################################################################
###############################################################################################


class QUICK_Solver(FTCS_Solver):
    
    # def Iteration_Formula(self, u_WW, u_W, u_E, u):
    #     self.Pe = 50
    #     Pe = self.Pe
    #     C = self.dt/self.dx
    #     r = ((1/Pe)*self.dt)/(self.dx**2)
    #     u_new = u -(C)*((u_E-u_W)/2-(u_E-3*u+3*u_W-u_WW)/6) +r*(u_W -2*u +u_E)
    #     #u_new = (6/8)*u_W +3/8*u-1/8*u_WW
    #     return u_new


    def Iteration_Formula(self, u_WW, u_W, u_E, u):
        self.Pe = 50
        Pe = self.Pe
        C = self.dt/self.dx
        r = ((1/Pe)*self.dt)/(self.dx**2)
        u_new = u -(C)*((3/8*u_E +6/8*u-1/8*u_W)-(3/8*u +6/8*u_W-1/8*u_WW)) +r*(u_W -2*u +u_E)

        # u_new = ((3/8*u_E +6/8*u-1/8*u_W)+(3/8*u +6/8*u_W-1/8*u_WW))/2 


        return u_new
    
    # def Edge_Formula(self, u_W, u_E, u):
    #     self.Pe = 50
    #     Pe = self.Pe
    #     C = self.dt/self.dx
    #     r = ((1/Pe)*self.dt)/(self.dx**2)

    #     u_new = (r-C)*u_E + (1-2*r)*u + (r+C)*u_W
    #     return u_new
    

    def Iterative(self):
        u_Next = copy.deepcopy(self.u)

        for n in range(0, self.n_max):
            u_Next[0] = 0
            # u_Next[1] = self.Edge_Formula(self.u[0],self.u[2],self.u[1])
            u_Next[1] = self.Iteration_Formula(self.u[1]-2*(self.u[0]-self.u[1]), self.u[0], self.u[2], self.u[1])
            for i in range(2,self.i_max):
                u_Next[i] = self.Iteration_Formula(self.u[i-2], self.u[i-1], self.u[i+1], self.u[i])
            u_Next[self.i_max] = 1
            self.u = u_Next

        return self.u
    
 

###############################################################################################


    


def Runner(dx, dt, x_max, t_max):

    u_Final = 0
    Run = QUICK_Solver(dx, dt, x_max, t_max)
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
    plt.title('QUICK scheme vs. Exact Final Solution at t=%f'%t_max)
    plt.show()



def u_exact(x):
    Pe = 50
    return (np.exp(x*Pe)-1)/(np.exp(Pe)-1)



    


class Post_op(FTCS_Solver): #UNFINISH!!!!! DO NOT FORGRT
    pass

        





def main():

    x_max = 1

    t_max = 1

    dx = 1/20
    dt = 0.001

    dx_20 = 1/20
    u_20 = Runner(dx_20, dt, x_max, t_max)

    dx_50 = 1/50
    u_50 = Runner(dx_50, dt, x_max, t_max)

    dx_100 = 1/100
    u_100 = Runner(dx_100, dt, x_max, t_max)

    plot_result(x_max, t_max, dx_20, dx_50, dx_100, u_20, u_50, u_100)





    







if __name__ == '__main__':
    main()













