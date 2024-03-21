

import math

import numpy as np


import matplotlib.pyplot as plt
import matplotlib as mpl

import copy

class UPwind1st_Solver:
    def __init__(self, dx, dt, x_max, t_max):
        self.dx = dx
        self.dt = dt
        self.x_max = x_max
        self.t_max = t_max

    def grid_generate(self):
        self.i_max = int(self.x_max/self.dx)
        self.n_max = int(self.t_max/self.dt)

        self.u = np.zeros(( self.i_max  ))

    
    def IC(self,m):
        self.m = m
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
    

    

class UPwind2nd_Solver(UPwind1st_Solver):
    def Iteration_Formula(self, u, u_W, u_WW):
        u_new = u + (3*u -4*u_W +u_WW)*self.dt/(2*self.dx)
      
        return u_new
    
    def Iterative(self):
        u_next = np.copy(self.u)

        for n in range(1, self.n_max):
            for i in range(0,self.i_max):
                u_next[i] = self.Iteration_Formula(self.u[i],self.u[i-1],self.u[i-2])
            self.u[:] = u_next[:]
        self.u_Full = self.getTHElastBACK(copy.deepcopy(self.u))
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





class UPwind3rd_Solver(UPwind1st_Solver):
    def Iteration_Formula(self, u, u_W, u_WW,u_WWW):
        u_new = u - ((11/6)*u -3*u_W +(3/2)*u_WW -(1/3)*u_WWW)*self.dt/(self.dx)
      
        return u_new
    
    def Iterative(self):
        u_next = copy.deepcopy(self.u)

        for n in range(1, self.n_max):
            for i in range(0,self.i_max):
                u_next[i] = self.Iteration_Formula(self.u[i],self.u[i-1],self.u[i-2],self.u[i-3])
            self.u[:] = np.copy(u_next)
        
        self.u_Full = self.getTHElastBACK(np.copy(self.u))
        return self.u
    

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
        u_e = copy.deepcopy(self.u_Full)
        # for i in range(0,self.i_max+1):
        #     u_e[i] = 0
        for i in range(0,self.i_max+1):
            u_e[i] = -math.sin(self.m*(self.t_max-i*self.dx))
        return u_e
    
    # def CompareTOexact(self):
    #     self.Error = np.zeros((self.i_max+1))
    #     u_e = self.exactsoln()
    #     for i in range(0,self.i_max+1):
    #             self.Error[i]+= abs(self.u_Full[i]-u_e[i])
    
    # def plot_error(self):
    #     x = np.arange(0,self.i_max+1)
    #     plt.plot(x, self.Error)
    #     plt.show()


    def plot_result_compare(self):
        x = np.arange(0,self.i_max+1)

        u_e =  self.exactsoln()

        print(self.u_Full) #TESTING ##############DELETE AFT TEST

        
        plt.plot(x*(2*math.pi/20),self.u_Full)
        plt.plot(x*(2*math.pi/20),u_e)

        print(u_e - self.u_Full)

        plt.show()





    


    

class Post_op(): #UNFINISH!!!!! DO NOT FORGRT
    pass

        





def main():

    x_max = 2*math.pi


    t_max = .5

    dx = 2*math.pi/20
    dt = 0.001

    m = 4


    upwind3 = UPwind3rd_Solver(dx, dt, x_max, t_max)

    upwind3.grid_generate()

    upwind3.IC(m)

    upwind3.Iterative()
    #upwind3.plot_result()
    upwind3.plot_result_compare()









if __name__ == '__main__':
    main()













