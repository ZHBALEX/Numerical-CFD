
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
        
            u_new, v_new = self.Time_Step(self.u, self.v)

            self.u, self.v = u_new, v_new
            print(n)
        
        return self.u, self.v
                 






    def Time_Step(self, u, v): # update u, v from n to n+1



        u_half, v_half = u.copy() , v.copy()
        u_new, v_new = u.copy() , v.copy()

        # Update n to n + 1/2
        for j in range(1,self.j_max):

            r = (self.mu*self.dt)/(2*self.dx**2)

            a = np.full(self.i_max+1, -r , dtype = 'float' )
            b = np.full(self.i_max+1, 1+2*r , dtype = 'float' )

            a[0] , a[self.i_max] = 0,0
            b[0] , b[self.i_max] = 1,1
            c = a

            d = np.full(self.i_max+1, 1 , dtype = 'float' )
            

            for i in range(1,self.i_max):
                cox = self.dt * u[i,j] / (2 * self.dx)
                coy = self.dt * v[i,j] / (2 * self.dy)
                ry = self.dt * self.mu / (2 * self.dy**2)
                d[i] = u[i,j] \
                    - cox * (u[i+1,j] - u[i-1,j])/2 \
                        - coy * (u[i,j+1] - u[i,j-1])/2 \
                            + ry * (u[i,j+1] + u[i,j-1] - 2 * u[i,j])
                d[0], d[-1]  = u[0,j], u[-1,j]



            u_half[:,j] = TDMA(a,b,c,d)

            for i in range(1,self.i_max):
                cox = self.dt * u[i,j] / (2 * self.dx)
                coy = self.dt * v[i,j] / (2 * self.dy)
                ry = self.dt * self.mu / (2 * self.dy**2)
                d[i] = u[i,j] \
                    - cox * (v[i+1,j] - v[i-1,j])/2 \
                        - coy * (v[i,j+1] - v[i,j-1])/2 \
                            + ry * (v[i,j+1] + v[i,j-1] - 2 * v[i,j])
                d[0], d[-1]  = v[0,j], v[-1,j]


            v_half[:,j] = TDMA(a,b,c,d)

        # Update n + 1/2 to n+1
        for i in range(1,self.i_max):
            r = (self.mu*self.dt)/(2*self.dy**2)

            a = np.full(self.j_max+1, -r , dtype = 'float' )
            b = np.full(self.j_max+1, 1+2*r , dtype = 'float' )

            a[0] , a[self.j_max] = 0,0
            b[0] , b[self.j_max] = 1,1
            c = a

            d = np.full(self.j_max+1, 1 , dtype = 'float' )

            
            for j in range(1,self.j_max):
                cox = self.dt * u[i,j] / (2 * self.dx)
                coy = self.dt * v[i,j] / (2 * self.dy)
                rx = self.dt * self.mu / (2 * self.dx**2)
                d[j] = u[i,j] \
                    - cox * (u[i+1,j] - u[i-1,j])/2 \
                        - coy * (u[i,j+1] - u[i,j-1])/2 \
                            + rx * (u[i+1,j] + u[i-1,j] - 2 * u[i,j])
                d[0], d[-1]  = u_half[i,0], u_half[i,-1]


            u_new[i,:] = TDMA(a,b,c,d)


            for j in range(1,self.j_max):
                cox = self.dt * u_half[i,j] / (2 * self.dx)
                coy = self.dt * v_half[i,j] / (2 * self.dy)
                rx = self.dt * self.mu / (2 * self.dx**2)
                d[j] = v[i,j] \
                    - cox * (v_half[i+1,j] - v_half[i-1,j])/2 \
                        - coy * (v_half[i,j+1] -v_half[i,j-1])/2 \
                            + rx * (v_half[i+1,j] + v_half[i-1,j] - 2 * v_half[i,j])
                d[0], d[-1]  = v_half[i,0], v_half[i,-1]
            v_new[j,:] = TDMA(a,b,c,d)

        return u_new , v_new





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
    plt.savefig('3d{:.3f}t{:.2f}.png'.format(d,t))
    # plt.show()

    


def main():
    x_max =2
    y_max =1
    t_max = 1.5

    # dx = 0.02
    # dy = 0.01
    # dt = 0.005


    dx = 0.005
    dy = 0.005
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

    t_max = 1.5
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





        





