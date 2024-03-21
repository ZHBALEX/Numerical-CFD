
import numpy as np #use it in all class u,v

import math # use it in Init to get initial value
import copy # use it in TDMA and Update to avoid influence origional data

import matplotlib.pyplot as plt  # use it to draw vorticity contour



def TDMA(a, b, c, d):  
          # TDMA solver 
          # b: main diagional  a: down diaginoal c: up diagional 
          # lenth of a , b, c, d is as same as D, with no use of a[0], a[max], c[0], c[max]
           
                # | b0  c0   0  ...   0   |
                # | a1  b1  c1  ...   0   |
                # |  0  a2  b2  ...   0   |
                # | ... ... ... ...  ...  |
                # |  0  ... an-2 bn-2 cn-2|
                # |  0  ...  0   an-1 bn-1|

    do = d.copy()
    ao = a.copy()
    bo = b.copy()
    co = c.copy()
    N = len(d)
    xo = np.zeros(N)
    
    for rowi in range(1,N):
        k = ao[rowi]/bo[rowi-1]
        bo[rowi] -= co[rowi-1]*k
        do[rowi] -= do[rowi-1]*k
    
    xo[N-1] = do[N-1]/bo[N-1]
    for rowi in range(N-2,-1,-1):
        xo[rowi] = (do[rowi]-co[rowi]*xo[rowi+1])/bo[rowi]

    return xo



class ADI_Solver:
    def __init__(self, x_max, y_max, t_max, dx, dy, dt, mu):
        self.x_max = x_max
        self.y_max = y_max
        self.t_max = t_max
        self.dx = dx
        self.dy = dy
        self.dt = dt
        self.mu = mu

    def  Grid_Generate(self):
        self.i_max = int(self.x_max/self.dx)
        self.j_max = int(self.y_max/self.dy)
        self.n_max = int(self.t_max/self.dt)

        self.u = np.zeros([ self.i_max+1 , self.j_max+1 ], dtype = float)
        self.v = np.zeros([ self.i_max+1 , self.j_max+1 ], dtype = float)


    def Initialize(self):
        for i in range(0, self.i_max+1):
            for j in range(0, self.j_max+1):
                self.u[i][j], self.v[i][j] = self.initial_formula(i, j)

        # print(self.u) ## TESTING 

    def initial_formula(self, i, j):
        V_t = 0.25
        x_0 = 0.5
        y_0 = 0.5
        r_0 = 0.1

        r2 = (i*self.dx-x_0)**2+(j*self.dy-y_0)**2

        u_init = 1 - V_t * (j*self.dy-y_0) * math.exp( (1-( r2 )/(r_0**2)) / (2) )
        v_init = V_t * (i*self.dx-x_0) * math.exp( (1-( r2 )/(r_0**2)) / (2) )

        return u_init, v_init
    

    def BC(self, u , v):

        for i in range(0, self.i_max+1):
            u[i][0] = 1
            v[i][0] = 0
            u[i][-1] = 1
            v[i][-1] = 0
        
        for j in range(0, self.j_max+1):
            u[0][j] = 1
            v[0][j] = 0
            u[-1][j] = 1
            v[-1][j] = 0
            # u[self.i_max-1][j], v[self.i_max-1][j] = u[self.i_max][j], v[self.i_max][j]
        return u, v
    

    def Iterative(self):
        self.u, self.v = self.BC(self.u, self.v)
        
        for n in range(0,self.n_max):
            self.u, self.v = self.BC(self.u, self.v)
        
            u_new, v_new = self.k_iteration()

            self.u, self.v = u_new, v_new

            print(n)
        
        return self.u, self.v
    
    
    def k_iteration(self):

        uk = self.u.copy()
        vk = self.v.copy()
        ry = self.dt * self.mu / ( self.dy**2) 
        rx = self.dt * self.mu / ( self.dx**2)

        error = 100
        error_detail = self.v.copy()

        k = 0
        
        while error > 1e-6:
            uk_new, vk_new = self.BC(uk.copy(),vk.copy())
            for i in range(1,self.i_max):
                for j in range(1,self.j_max):
                    cox = self.dt * self.u[i,j] / ( self.dx)
                    coy = self.dt * self.v[i,j] / ( self.dy)
                     
                    uk_new[i,j] = self.u[i,j] \
                        -cox*(self.u[i+1,j]-self.u[i-1,j])/2 \
                            -coy*(self.u[i,j+1]-self.u[i,j-1])/2 \
                            + rx * (uk[i+1,j]+uk[i-1,j]-2*uk[i,j]) \
                            + ry * (uk[i,j+1]+uk[i,j-1]-2*uk[i,j])
                    
                    error_detail[i,j] = -uk_new[i,j] + self.u[i,j] \
                        -cox*(self.u[i+1,j]-self.u[i-1,j])/2 \
                            -coy*(self.u[i,j+1]-self.u[i,j-1])/2 \
                            + rx * (uk_new[i+1,j]+uk_new[i-1,j]-2*uk_new[i,j]) \
                            + ry * (uk_new[i,j+1]+uk_new[i,j-1]-2*uk_new[i,j])
            error = np.sum(np.abs(error_detail))
            uk = uk_new
            # k+=1 
            # print(k)
        self.u = uk_new

        error = 100

        k = 0

        while error > 1e-6:
            uk_new = uk.copy()
            vk_new = vk.copy()
            for i in range(1,self.i_max):
                for j in range(1,self.j_max):
                    # cox = self.dt * self.u[i,j] / ( self.dx)
                    # coy = self.dt * self.v[i,j] / ( self.dy)
                     
                    vk_new[i,j] = self.v[i,j] \
                        -cox*(self.v[i+1,j]-self.v[i-1,j])/2 \
                            -coy*(self.v[i,j+1]-self.v[i,j-1])/2 \
                            + rx * (vk[i+1,j]+vk[i-1,j]-2*vk[i,j]) \
                            + ry * (vk[i,j+1]+vk[i,j-1]-2*vk[i,j])
                    
                    error_detail[i,j] = -vk_new[i,j] + self.v[i,j] \
                        -cox*(self.v[i+1,j]-self.v[i-1,j])/2 \
                            -coy*(self.v[i,j+1]-self.v[i,j-1])/2 \
                            + rx * (vk_new[i+1,j]+vk_new[i-1,j]-2*vk_new[i,j]) \
                            + ry * (vk_new[i,j+1]+vk_new[i,j-1]-2*vk_new[i,j])
            error = np.sum(np.abs(error_detail))
            vk = vk_new
            # k+= 1
            # print(k)
        self.v = vk_new



        self.u, self.v = self.BC(self.u,self.v)

        print(error)



        return self.u, self.v
            



def Do_Solver(x_max, y_max, t_max, dx, dy, dt , mu):
     
    Try_1  = ADI_Solver(x_max, y_max, t_max, dx, dy, dt , mu)
    Try_1.Grid_Generate()
    Try_1.Initialize()
    u, v =Try_1.Iterative()

    return u,v


def Vorticity(u, v, t):
    # u,v = np.transpose(u), np.transpose(v)

    Ni, Nj = u.shape

    x = np.linspace(0, 2, Ni)
    y = np.linspace(0, 1, Nj)

    X, Y = np.meshgrid(x, y)

    omega = np.zeros((Ni, Nj))
    # omega[1:-1, 1:-1] = (v[1:-1, 2:] - v[1:-1, :-2]) / (x[2] - x[1]) - (u[2:, 1:-1] - u[:-2, 1:-1]) / (y[2] - y[1])

    for i in range(1,Ni-1):
         for j in range(1,Nj-1):
              omega[i,j] = (v[i+1,j] - v[i-1,j])/(x[2]-x[0])-(u[i,j+1] - u[i,j-1])/(y[2] - y[0])

    # omega = v

    omega = np.transpose(omega)





    d = x[2]-x[1]

    print(d)

    plt.figure(figsize=(12, 6))
    plt.contourf(X, Y, omega, levels=50, cmap='jet')
    plt.colorbar()
    plt.title('d={:.3f},t={:.2f}'.format(d,t))
    plt.xlabel('x')
    plt.ylabel('y')


    plt.tight_layout()
    # plt.savefig('3d{:.3f}t{:.2f}.png'.format(d,t))
    plt.show()

    


def main():
    x_max =2
    y_max =1
    t_max = 0.02

    # dx = 0.02
    # dy = 0.01
    # dt = 0.005


    dx = 0.02
    dy = 0.02
    dt = 0.001


    mu = 0.01
    # mu = 1


    # t_max = 0
    # u,v = Do_Solver(x_max, y_max, t_max, dx, dy, dt , mu)
    # # print(u)
    # Vorticity(u, v, t_max)

    # t_max = 0.5
    # u,v = Do_Solver(x_max, y_max, t_max, dx, dy, dt , mu)
    # # print(u)
    # Vorticity(u, v, t_max)

    # t_max = 1
    # u,v = Do_Solver(x_max, y_max, t_max, dx, dy, dt , mu)
    # # print(u)
    # Vorticity(u, v, t_max)

    t_max = 1
    u,v = Do_Solver(x_max, y_max, t_max, dx, dy, dt , mu)
    # print(u)
    Vorticity(u, v, t_max)

    # t_max = 2
    # u,v = Do_Solver(x_max, y_max, t_max, dx, dy, dt , mu)
    # # print(u)
    # Vorticity(u, v, t_max)



if __name__ == '__main__':
     main()
    # arg1 = [2,2,2,2,2]
    # arg2 = [1,1,1,1,1]
    # arg3 = arg1
    # arg4 = [1,2,3,4,5]
    # print(TDMA(arg1,arg2,arg3,arg4))





        





