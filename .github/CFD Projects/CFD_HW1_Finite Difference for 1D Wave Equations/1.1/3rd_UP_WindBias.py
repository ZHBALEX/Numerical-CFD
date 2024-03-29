

import math

import numpy as np


import matplotlib.pyplot as plt
import matplotlib as mpl

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
        u_new = u + (math.pi/(2*10**4))*(u_W-u)
        return u_new

    def Iterative(self):

        for n in range(1, self.n_max):
            for i in range(0,self.i_max):
                self.u[i] = self.Iteration_Formula(self.u[i-1], self.u[i])

        return self.u
    
    def plot_result(self):
        x = np.arange(0,self.i_max)
        plt.plot(x,self.u)
        plt.show()
    
    

class UPwind3rdBias_Solver(UPwind1st_Solver):
    def Iteration_Formula(self, u, u_W, u_WW,u_E):
        u_new = u - ((3)*u -6*u_W +1*u_WW +2*u_E)*self.dt/(6*self.dx)
      
        return u_new
    
    def Iterative(self):
        u_next = copy.deepcopy(self.u)

        for n in range(1, self.n_max+1):
            for i in range(0,self.i_max):
                u_next[i] = self.Iteration_Formula(self.u[i],self.u[i-1],self.u[i-2],self.u[(i+1)%(self.i_max)])
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

    
    # def exactsoln(self):
    #     u_e = np.zeros(self.i_max +1)

    #     for i in range(0,self.i_max+1):
    #         u_e[i] = -math.sin(self.m*(self.t_max-i*self.dx))
    #     return u_e
    
    def exactsoln(self):
        u_e = np.zeros(self.i_max+1)
        for i in range(self.i_max+1):
            # 这里x是根据空间位置i计算的，确保与波的移动正确对应
            x = i * self.dx
            u_e[i] = math.sin(self.m * (x - self.t_max))
        return u_e

    def CompareTOexact(self):
        self.Error = np.zeros((self.i_max+1))
        u_e = self.exactsoln()
        for i in range(0,self.i_max+1):
                self.Error[i]+= abs(self.u_Full[i]-u_e[i])
    
    def plot_error(self):
        x = np.arange(0,self.i_max+1)
        plt.plot(x, self.Error)
        plt.show()


    def plot_result_compare(self):
        x = np.arange(0,self.i_max+1)

        u_e =  self.exactsoln()

        
        plt.plot(x*(2*math.pi/20),self.u_Full)
        plt.plot(x*(2*math.pi/20),u_e)

        print(u_e - self.u_Full)

        plt.show()






    


    

class Post_op(): #UNFINISH!!!!! DO NOT FORGRT
    
    pass

        





def main():

    x_max = 2*math.pi


    t_max = 0.05

    dx = 2*math.pi/20
    dt = 0.001

    m = 2



    


    # 3rd bias


    upwind3b = UPwind3rdBias_Solver(dx, dt, x_max, t_max,m)

    upwind3b.grid_generate()

    upwind3b.IC(m)

    u = upwind3b.Iterative()
    print(u)
    upwind3b.plot_result_compare()
    
    
    #upwind3b.CompareTOexact()
    #upwind3b.plot_error()










if __name__ == '__main__':
    main()













