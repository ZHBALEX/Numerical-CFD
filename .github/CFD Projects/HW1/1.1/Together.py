

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



class UPwind2nd_Solver(UPwind1st_Solver):
    def Iteration_Formula(self, u, u_W, u_WW):
        u_new = u - (3*u -4*u_W +u_WW)*self.dt/(2*self.dx)
      
        return u_new
    
    def Iterative(self):
        u_next = np.copy(self.u)

        for n in range(1, self.n_max+1):
            for i in range(0,self.i_max):
                u_next[i] = self.Iteration_Formula(self.u[i],self.u[i-1],self.u[i-2])
            self.u[:] = u_next[:]
        self.u_Full = self.getTHElastBACK(self.u)
        return self.u_Full




class UPwind3rd_Solver(UPwind1st_Solver):
    def Iteration_Formula(self, u, u_W, u_WW,u_WWW):
        u_new = u - ((11/6)*u -3*u_W +(3/2)*u_WW -(1/3)*u_WWW)*self.dt/(self.dx)
      
        return u_new
    
    def Iterative(self):
        u_next = np.copy(self.u)

        for n in range(1, self.n_max+1):
            for i in range(0,self.i_max):
                u_next[i] = self.Iteration_Formula(self.u[i],self.u[i-1],self.u[i-2],self.u[i-3])
            self.u[:] = np.copy(u_next)
        
        self.u_Full = self.getTHElastBACK(np.copy(self.u))
        return self.u_Full
    



    
    

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
    


class Total_Compare(UPwind1st_Solver):
    def __init__(self, dx, dt, x_max, t_max, m, UPwind1, UPwind2, UPwind3, UPwind3b):
        super().__init__(dx, dt, x_max, t_max, m)
        self.UPwind1 = UPwind1
        self.UPwind2 = UPwind2
        self.UPwind3 = UPwind3
        self.UPwind3b = UPwind3b



    def plot_result_compare(self):
        x = np.arange(0,self.i_max+1)

        u_e =  self.exactsoln()

        plt.plot(x*(2*math.pi/20),self.UPwind1, label = "UPwind1")
        plt.plot(x*(2*math.pi/20),self.UPwind2, label = "UPwind2")
        plt.plot(x*(2*math.pi/20),self.UPwind3, label = "UPwind3")
        plt.plot(x*(2*math.pi/20),self.UPwind3b, label = "UPwind3b")
        plt.plot(x*(2*math.pi/20),u_e, label = "Exact")

        plt.title('Solutions at t_max=%.1f, m=%d' % (self.t_max, self.m))
        plt.xlabel('x')
        plt.ylabel('u')

        plt.legend()
        #plt.figure(figsize = (5,1))

        plt.show()

    def Error_Compute(self, u, u_e):
        Error = np.zeros(self.i_max+1)
        for i in range(0,self.i_max+1):
                Error[i]+= abs(u[i]-u_e[i])
        return Error

    def plot_Error_compare(self):
        x = np.arange(0,self.i_max+1)

        u_e =  self.exactsoln()

        EUPwind1 = self.Error_Compute(self.UPwind1, u_e)
        EUPwind2 = self.Error_Compute(self.UPwind2, u_e)
        EUPwind3 = self.Error_Compute(self.UPwind3, u_e)
        EUPwind3b = self.Error_Compute(self.UPwind3b, u_e)

        plt.plot(x*(2*math.pi/20),EUPwind1, label = "UPwind1")
        plt.plot(x*(2*math.pi/20),EUPwind2, label = "UPwind2")
        plt.plot(x*(2*math.pi/20),EUPwind3, label = "UPwind3")
        plt.plot(x*(2*math.pi/20),EUPwind3b, label = "UPwind3b")

        plt.title('Error at t_max=%.1f, m=%d' % (self.t_max, self.m))
        plt.xlabel('x')
        plt.ylabel('u')

        plt.legend()
        plt.show()


    




        





def main():

    x_max = 2*math.pi


    t_max = 1

    dx = 2*math.pi/20
    dt = 0.001

    m = 2


    upwind1st = UPwind1st_Solver(dx, dt, x_max, t_max,m)

    upwind1st.grid_generate()

    upwind1st.IC(m)

    U1 = upwind1st.Iterative()



    upwind2nd = UPwind2nd_Solver(dx, dt, x_max, t_max,m)

    upwind2nd.grid_generate()

    upwind2nd.IC(m)

    U2 = upwind2nd.Iterative()



    upwind3 = UPwind3rd_Solver(dx, dt, x_max, t_max, m)

    upwind3.grid_generate()

    upwind3.IC(m)

    U3 = upwind3.Iterative()
    




    


    # 3rd bias


    upwind3b = UPwind3rdBias_Solver(dx, dt, x_max, t_max,m)

    upwind3b.grid_generate()

    upwind3b.IC(m)

    U3b = upwind3b.Iterative()





    Total = Total_Compare(dx, dt, x_max, t_max, m, U1, U2, U3, U3b)
    Total.grid_generate()

    Total.IC(m)

    Total.plot_result_compare()
    #Total.plot_Error_compare()
    










if __name__ == '__main__':
    main()













