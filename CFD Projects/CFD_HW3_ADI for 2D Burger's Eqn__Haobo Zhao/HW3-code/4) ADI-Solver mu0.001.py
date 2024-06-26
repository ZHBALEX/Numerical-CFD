
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

    import copy
    do = copy.deepcopy(d)
    ao = copy.deepcopy(a)
    bo = copy.deepcopy(b)
    co = copy.deepcopy(c)
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
            
            u_new, v_new = copy.deepcopy(self.u), copy.deepcopy(self.v)
            
            u_new, v_new = self.K_iteration(self.u, self.v)

            self.u, self.v = u_new, v_new
            print(n)
        
        return self.u, self.v
                 






    def K_iteration(self, u, v):

        u_half, v_half = copy.deepcopy(u), copy.deepcopy(v)
        u_new, v_new =  copy.deepcopy(u), copy.deepcopy(v)


            # update u:
                # u{n} to u{n+1/2} line update:
        for j in range(1, self.j_max):  # flag=flag(0/1,0/1), 1:donig 0:other
            flag = (0,1)
            u_half[:,j] = self.UPdate( u , v, u, j, self.i_max, flag)
                # Update function: Update(u, v, doing(u/v), doing(row=j/column=i), doingAnother, flag)
            


        u_new, v_new = self.BC(u_new, v_new)
        u_half, v_half = self.BC(u_half, v_half )
            
        
            # update v:
        for j in range(1, self.j_max):  # flag=flag(0/1,0/1), 1:donig 0:other
            flag = (0,1)
            v_half[:,j] = self.UPdate( u_half , v, v, j, self.i_max, flag)
        

        u_new, v_new = self.BC(u_new, v_new)
        u_half, v_half = self.BC(u_half, v_half )

                
                # u{n+1/2} to u{n+1} column update:
        for i in range(1,self.i_max-1):
            flag = (1,0)
            u_new[i,:] = self.UPdate( u_half , v_half, u_half, i, self.j_max, flag)
        

        u_new, v_new = self.BC(u_new, v_new)
        u_half, v_half = self.BC(u_half, v_half )
            


        for i in range(1,self.i_max):
            flag = (1,0)
            v_new[i,:] = self.UPdate( u_new , v_half , v_half, i, self.j_max, flag)
        

        u_new, v_new = self.BC(u_new, v_new)
        u_half, v_half = self.BC(u_half, v_half )

        return u_new , v_new
    



    # def UPdate(self, u, v, U,  Idoing, Ianother_max, flag):

    #     delta = self.dx * flag[1] + self.dy*flag[0]


    #     r = (self.mu*self.dt)/(2*delta**2)


    #     # Creating a, b, c, D
    #     a = np.full(Ianother_max-1, -r, dtype = float)
    #     b = np.full(Ianother_max-1, 1 + 4*r, dtype = float)
    #     c = np.full(Ianother_max-1, -r, dtype = float)
    #     d = np.full(Ianother_max+1, 0, dtype = float)


    #     # updating D form [1,max-1]:
    #     if flag == (0, 1):
    #         j = Idoing
    #         for i in range(1,self.i_max):
    #             cox = (u[i,j]*self.dt/self.dx /2)
    #             coy = (v[i,j]*self.dt/self.dy /2)
               
    #             ry = (self.mu*self.dt/(2*self.dy**2))
    #             d[i] = U[i,j]-cox*((U[i+1,j]-U[i-1,j])/(2))-coy*((U[i,j+1]-U[i,j-1])/(2))+ry*(U[i,j+1]+U[i,j-1])
    #         d[0], d[-1] = U[0,j], U[-1,j]
                
    #     elif flag == (1,0):
    #          i = Idoing
    #          for j in range(1,self.j_max):
    #             cox = (u[i,j]*self.dt/self.dx /2)
    #             coy = (v[i,j]*self.dt/self.dy /2)
               
    #             rx = (self.mu*self.dt/(2*self.dx**2))
    #             d[j] = U[i,j]-cox*((U[i+1,j]-U[i-1,j])/(2))-coy*((U[i,j+1]-U[i,j-1])/(2))+rx*(U[i+1,j]+U[i-1,j])
    #          d[0], d[-1] = U[i,0], U[i,-1]
        
    #     # added boundary get full matrix:
    #     a = np.concatenate((np.array([0]), a, np.array([0])), axis=0)
    #     b = np.concatenate((np.array([1]), b, np.array([1])), axis=0)
    #     c = np.concatenate((np.array([0]), c, np.array([0])), axis=0)


           
    #     u_UPdate = TDMA(a, b, c, d)

    #         # | 1   0    0   ...   0      0    |
    #         # | a1  b1  c1   ...   0      0    |
    #         # |  0  a2  b2   ...   0      0    |
    #         # | ... ... ...  ...  ...    ...   |
    #         # |  0  0    0   ...  b(n-2) c(n-2)|
    #         # |  0  0    0   ...   0      1    |
        
    #     # print("TDMA")
    #     # print(u_UPdate)

    #     return u_UPdate

    def UPdate(self, u, v, U,  Idoing, Ianother_max, flag):

        delta = self.dx * flag[1] + self.dy*flag[0]


        r = (self.mu*self.dt)/(2*delta**2)


        # Creating a, b, c, D
        a = np.full(Ianother_max-1, -r, dtype = float)
        b = np.full(Ianother_max-1, 1 + 4*r, dtype = float)
        c = np.full(Ianother_max-1, -r, dtype = float)
        d = np.full(Ianother_max+1, 0, dtype = float)

        for ij in range(1,Ianother_max):
            i,j = Idoing*flag[0] + ij* flag[1] ,Idoing*flag[1] + ij* flag[0]
            cox = (u[i,j]*self.dt/self.dx /2)
            coy = (v[i,j]*self.dt/self.dy /2)
            rx = (self.mu*self.dt/(2*self.dx**2))
            ry = (self.mu*self.dt/(2*self.dy**2))
            rxy = rx * flag[1] + ry*flag[0]
            d[ij] = U[i,j]-cox*((U[i+1,j]-U[i-1,j])/(2))-coy*((U[i,j+1]-U[i,j-1])/(2))+rxy*(U[i+flag[0],j+flag[1]]+U[i-flag[0],j-flag[1]])
        d[0], d[-1] = U[i*flag[0],j*flag[1]], U[i*flag[0]-1,j*flag[1]-1]
        
        
        # added boundary get full matrix:
        a = np.concatenate((np.array([0]), a, np.array([0])), axis=0)
        b = np.concatenate((np.array([1]), b, np.array([1])), axis=0)
        c = np.concatenate((np.array([0]), c, np.array([0])), axis=0)


           
        u_UPdate = TDMA(a, b, c, d)

            # | 1   0    0   ...   0      0    |
            # | a1  b1  c1   ...   0      0    |
            # |  0  a2  b2   ...   0      0    |
            # | ... ... ...  ...  ...    ...   |
            # |  0  0    0   ...  b(n-2) c(n-2)|
            # |  0  0    0   ...   0      1    |
        
        # print("TDMA")
        # print(u_UPdate)

        return u_UPdate






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
    plt.title('mu = 0.001 d={:.3f},t={:.2f}'.format(d,t))
    plt.xlabel('x')
    plt.ylabel('y')


    plt.tight_layout()
    plt.savefig('mu3d{:.3f}t{:.2f}.png'.format(d,t))
    # plt.show()

    


def main():
    x_max =2
    y_max =1
    t_max = 1

    # dx = 0.02
    # dy = 0.01
    # dt = 0.005


    dx = 0.01
    dy = 0.01
    dt = 0.001


    mu = 0.001
    # mu = 1

    t_max = 0
    u,v = Do_Solver(x_max, y_max, t_max, dx, dy, dt , mu)
    # print(u)
    Vorticity(u, v, t_max)

    t_max = 0.5
    u,v = Do_Solver(x_max, y_max, t_max, dx, dy, dt , mu)
    # print(u)
    Vorticity(u, v, t_max)

    t_max = 1
    u,v = Do_Solver(x_max, y_max, t_max, dx, dy, dt , mu)
    # print(u)
    Vorticity(u, v, t_max)

    t_max = 1.5
    u,v = Do_Solver(x_max, y_max, t_max, dx, dy, dt , mu)
    # print(u)
    Vorticity(u, v, t_max)

    t_max = 2
    u,v = Do_Solver(x_max, y_max, t_max, dx, dy, dt , mu)
    # print(u)
    Vorticity(u, v, t_max)



if __name__ == '__main__':
     main()
    # arg1 = [2,2,2,2,2]
    # arg2 = [1,1,1,1,1]
    # arg3 = arg1
    # arg4 = [1,2,3,4,5]
    # print(TDMA(arg1,arg2,arg3,arg4))





        





