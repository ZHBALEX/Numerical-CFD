

import numpy as np
import matplotlib.pyplot as plt
plt.ioff()  # 关闭交互模式
import copy
import re


class mesher:
    def __init__(self, data):
        self.Lx = data["Lx"]
        self.Ly = data["Ly"]
        # self.Re = data["Re"]
        self.dxy = data["dx"]

    def parameter_grid(self): # generate large grid 
        Field1 = grid_generator( self.Ly,  self.Lx, self.dxy)
        N_rows, N_cols,a = np.shape(Field1)

        # save location y,x in field[j,i,0], 
        for j in range(0, N_rows):
            for i in range(0,N_cols):
                Field1[j,i,0] , Field1[j,i,1] = (j-1/2)*self.dxy, (i-1/2)*self.dxy
        # print(Field1[:,:,1])
        # save fluid(1)/Solid(0) in field[j,i,2]
        Field1[:,:,2] = i_fluid_detector(Field1[:,:,0], Field1[:,:,1], Field1[:,:,2])
        Field1[:,:,3] = 1 - Field1[:,:,2]

        # Field1.flags.writeable = False
        return Field1
    
def grid_generator( y_size,  x_size, grid_size):
        print(1.1)
        Nx , Ny = int(x_size/grid_size), int(y_size/grid_size)
        print(1.2)
        print(Ny,Nx)
        Field = np.full((Ny+2, Nx+2, 4), 0 ,dtype = 'float')
        print(1.3)
        return Field


def i_fluid_detector(field_y, field_x, field):
    N_rows, N_cols = np.shape(field)
    for j in range(0, N_rows):
        for i in range(0, N_cols ):
            result = 0
            result = circle_function(field_y[j,i], field_x[j,i])
            if result < 0:
                field[j,i] = 0
            else:
                field[j,i] = 1
    return field

def circle_function(y,x):
    R= 0.25
    x_c = data['Lx']*1/4
    y_c = data['Ly']*1/2
    result = (x-x_c)**2 + (y-y_c)**2 - R**2
    return result
        
        
def IBM_uv_coeff_change(c_W, c_E, c_S, c_N, c_P, iF, iS): # iF: ifluid. iS: isolid
        c_P[1:-1,1:-1] = c_P[1:-1,1:-1] -(
            iS[1:-1,:-2]*c_W[1:-1,1:-1] +
            iS[1:-1,2:]*c_E[1:-1,1:-1] +
            iS[:-2,1:-1]*c_S[1:-1,1:-1] +
            iS[2:,1:-1]*c_N[1:-1,1:-1] ) 
        
        c_W[1:-1,1:-1] = iF[1:-1,:-2] * c_W[1:-1,1:-1]
        c_E[1:-1,1:-1] = iF[1:-1,2:] * c_E[1:-1,1:-1]
        c_S[1:-1,1:-1] = iF[:-2,1:-1] * c_S[1:-1,1:-1]
        c_N[1:-1,1:-1] = iF[2:,1:-1] * c_N[1:-1,1:-1]
        
        
        c_P[1:-1,1:-1] = iF[1:-1,1:-1] * c_P[1:-1,1:-1] + iS[1:-1,1:-1]
        c_E[1:-1,1:-1] = Solid__coeff_correct(c_E, iF)
        c_W[1:-1,1:-1] = Solid__coeff_correct(c_W, iF)
        c_N[1:-1,1:-1] = Solid__coeff_correct(c_N, iF)
        c_S[1:-1,1:-1] = Solid__coeff_correct(c_S, iF)
        return c_W, c_E, c_S, c_N, c_P
 

def IBM_p_coeff_change(c_W, c_E, c_S, c_N, c_P, iF, iS): # iF: ifluid. iS: isolid
        c_P[1:-1,1:-1] = c_P[1:-1,1:-1] +(
            iS[1:-1,:-2]*c_W[1:-1,1:-1] +
            iS[1:-1,2:]*c_E[1:-1,1:-1] +
            iS[:-2,1:-1]*c_S[1:-1,1:-1] +
            iS[2:,1:-1]*c_N[1:-1,1:-1] ) 
        c_W[1:-1,1:-1] = iF[1:-1,:-2] * c_W[1:-1,1:-1]
        c_E[1:-1,1:-1] = iF[1:-1,2:] * c_E[1:-1,1:-1]
        c_S[1:-1,1:-1] = iF[:-2,1:-1] * c_S[1:-1,1:-1]
        c_N[1:-1,1:-1] = iF[2:,1:-1] * c_N[1:-1,1:-1]
        
        c_P[1:-1,1:-1] = iF[1:-1,1:-1] * c_P[1:-1,1:-1] + iS[1:-1,1:-1]
        c_E[1:-1,1:-1] = Solid__coeff_correct(c_E, iF)
        c_W[1:-1,1:-1] = Solid__coeff_correct(c_W, iF)
        c_N[1:-1,1:-1] = Solid__coeff_correct(c_N, iF)
        c_S[1:-1,1:-1] = Solid__coeff_correct(c_S, iF)
        return c_W, c_E, c_S, c_N, c_P

def Solid__coeff_correct(C, iF):
        return np.copy(C[1:-1,1:-1] * iF[1:-1,1:-1])

def TDMA(a, b, c, d):  # TDMA Solver, input abcd in same length
    # a[0], c[-1] not been used
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


def general_grid_generator(len_x, Nx, Ny):
    # Creating same size of u, v, U, V
    # where col 0, row 0, row -1, of U not used
    # where row 0, col 0, col -1, of V not used
    initial_Value = 0.01

    dx = len_x/Nx
    grid = np.full((Ny+2, Nx+2), initial_Value)
    grid[0,0] = grid[-1,-1] = grid[-1,0] = grid[0,-1] = np.inf
    u,v,U,V = data_copy(grid, grid, grid, grid)
    p = np.copy(grid)
    return u, v, U, V, p, dx

def data_copy(u, v, U, V):
    return np.copy(u), np.copy(v), np.copy(U), np.copy(V)



def Ghost_BC_p(p_input): # updating ghost point for p based on BC condition
    p = np.copy(p_input)
    p[-1,1:-1] =   p[-2,1:-1] # upper BC
    p[0, 1:-1] =   p[1,1:-1] # bottom BC
    p[1:-1, 0] =   p[1:-1,1] # left BC
    p[1:-1, -1] =  p[1:-1,-2] # right BC
    return p

def Ghost_BC_u(u): # updating ghost point for p based on BC condition
    u[-1,1:-1] =   u[-2,1:-1] # upper BC
    u[0, 1:-1] =   u[1,1:-1] # bottom BC
    u[1:-1, 0] =   2-u[1:-1,1] # left BC
    u[1:-1, -1] =   u[1:-1,-2] # right BC
    u = u[:,:]*iF[:,:]
    return u

def Ghost_BC_v(v): # updating ghost point for p based on BC condition
    v[-1,1:-1] =   -v[-2,1:-1] # upper BC
    v[0, 1:-1] =   -v[1,1:-1] # bottom BC
    v[1:-1, 0] =   -v[1:-1,1] # left BC
    v[1:-1, -1] =   v[1:-1,-2] # right BC
    v = v[:,:]*iF[:,:]
    return v

def UV_BC(U, V, u, iF):
    U[1:-1,1] = 1 # left
    # U[1:-1,-1] = np.copy(u[1:-1, -2]) # right
    V[1,:] = 0 # bottom
    V[-1,:] = 0 # upper
    # V[:,-1] = V[:,-2]

    U[1:-1,1:-1] = U[1:-1,1:-1] * iF[1:-1,1:-1] 
    U[1:-1,2:] = U[1:-1,2:] * iF[1:-1,1:-1] 

    V[1:-1,1:-1] = V[1:-1,1:-1] * iF[1:-1,1:-1] 
    V[2:,1:-1] = V[2:,1:-1] * iF[1:-1,1:-1]
    return U,V



def coeff_assemble(rows_cols, a_e, b_e, c_e, e_e, f_e):
    rows, cols = rows_cols[0], rows_cols[1]
    # input elements of abcdef, get line matrix of abcdef
    a = np.full((rows, cols), a_e ) # abc in same length of Ncols of u
    b = np.full((rows, cols), b_e )
    c = np.full((rows, cols), c_e )

    e = np.full((rows, cols), e_e ) # ef in same length of Nrows of u
    f = np.full((rows, cols), f_e )
    a[:,0]= a[:,-1] = c[:,0] = c[:,-1] = 0
    b[:,0] = b[:,-1] = 1
    e[:,-1] =e[:,0] = f[:,-1]= f[:,0]=0
    e[0,:] = e[-1,:] = f[0,:] = f[-1,:] = 0
    a[0,:] = a[-1,:] = c[0,:] = c[-1,:] = 0
    return a,b,c, e,f

def RHS_uv(f, f_old, U, V, U_old, V_old, c0, r): 
    # Update Right Hand Side
    f_P = f[1:-1,1:-1]
    f_W = f[1:-1,:-2]
    f_E = f[1:-1,2:]
    f_S = f[:-2,1:-1]
    f_N = f[2:,1:-1]

    f_P_old = f_old[1:-1,1:-1]
    f_W_old = f_old[1:-1,:-2]
    f_E_old = f_old[1:-1,2:]
    f_S_old = f_old[:-2,1:-1]
    f_N_old = f_old[2:,1:-1]

    convection_new = c0*( U[1:-1,2:] * (f_E + f_P)/2 -
                          U[1:-1,1:-1] * (f_W + f_P)/2 +
                          V[2:,1:-1]  * (f_N + f_P)/2 -
                          V[1:-1,1:-1] * (f_S + f_P)/2 )
    
    convection_old = c0*( U_old[1:-1,2:] * (f_E_old + f_P_old)/2 -
                          U_old[1:-1,1:-1] * (f_W_old + f_P_old)/2 +
                          V_old[2:,1:-1]  * (f_N_old + f_P_old)/2 -
                          V_old[1:-1,1:-1] * (f_S_old + f_P_old)/2 )
    
    diffusion = f_P + r*(f_E + f_W + f_N + f_S -4*f_P)
    return -3/2*convection_new + 1/2*convection_old + diffusion


def Residual(u, a, b, c, d, e, f):
    return (a[1:-1,1:-1] * u[1:-1,:-2] +
                      b[1:-1,1:-1] * u[1:-1,1:-1] +
                      c[1:-1,1:-1] * u[1:-1,2:] +
                      e[1:-1,1:-1] * u[:-2,1:-1] +
                      f[1:-1,1:-1] * u[2:,1:-1] -
                      d[:,:])

def LineSOR_u(a,b,c,d,e,f, u): 
    # input abcdef, and u, get :D
    # then use TDMA get new u
    # abc is in same size of N_rows of u, ef in same size of N_cols of u
    # ! Input d is LIKE inner field of u, need to cut each line to use !
    u_k = np.copy(u)
    u_k_new = np.copy(u)
    n_rows = np.size(u,0)
    w = 1.3
    D = np.zeros_like(a[1,:])
    Res = 100
    k = 0
    while Res>1e-6:
        u_k = np.copy(u_k_new) 
        for j in range(1,n_rows-1): # the number is based on d, which is (N-1)x(N-1)
            # D = np.zeros_like(a[1,:])
            D[1:-1] = d[j-1,:] - e[j,1:-1]*u_k_new[j-1,1:-1] - f[j,1:-1]*u_k_new[j+1,1:-1]
            D[0], D[-1] = 2-u_k_new[j,1], u_k_new[j,-2] # ghost point = ghost point
            u_k_new[j,:] = TDMA(a[j,:],b[j,:],c[j,:],D) # TDMA solve uk each line
            u_k_new[j,1:-1] = u_k[j,1:-1]*(1-w) +u_k_new[j,1:-1]*w # SOR term
        u_k_new = Ghost_BC_u(u_k_new)
        Res = np.max(np.abs(Residual(u_k_new, a, b, c, d, e, f)))
        # print("Res_u=")
        # print(Res) 
        u = np.copy(u_k)
        # print(Residual(u_k_new, a, b, c, d, e, f))
        k+=1
    # print("finish iteration")
    return u_k_new, k, Res

def LineSOR_v(a,b,c,d,e,f, u): 
    # input abcdef, and u, get D
    # then use TDMA get new u
    # abc is in same size of N_rows of u, ef in same size of N_cols of u
    # ! Input d is LIKE inner field of u, need to cut each line to use !
    u_k = np.copy(u)
    u_k_new = np.copy(u)
    n_rows = np.size(u,0)
    w = 1.3
    D = np.zeros_like(a[1,:])
    Res = 100
    k = 0
    while Res>1e-6:
        u_k = np.copy(u_k_new) 
        for j in range(1,n_rows-1): # the number is based on d, which is (N-1)x(N-1)
            # D = np.zeros_like(a[1,:])
            D[1:-1] = d[j-1,:] - e[j,1:-1]*u_k_new[j-1,1:-1] - f[j,1:-1]*u_k_new[j+1,1:-1]
            D[0], D[-1] = -u_k_new[j,1], u_k_new[j,-2] # ghost point = ghost point
            u_k_new[j,:] = TDMA(a[j,:],b[j,:],c[j,:],D) # TDMA solve uk each line
            u_k_new[j,1:-1] = u [j,1:-1]*(1-w) +u_k_new[j,1:-1]*w # SOR term
        u_k_new = Ghost_BC_v(u_k_new)
        Res = np.max(np.abs(Residual(u_k_new, a, b, c, d, e, f)))
        # print("Res_v=v")
        # print(Res) 
        u = np.copy(u_k)
        k+=1
        # print(Residual(u_k_new, a, b, c, d, e, f))
    # print("finish iteration")
    return u_k_new, k, Res

def LineSOR_p(a,b,c,d,e,f, u): 
    # input abcdef, and u, get D
    # then use TDMA get new u
    # abc is in same size of N_rows of u, ef in same size of N_cols of u
    # ! Input d is LIKE inner field of u, need to cut each line to use !
    u = Ghost_BC_p(u)
    u_k = np.copy(u)
    u_k_new = np.copy(u)
    n_rows = np.size(u,0)
    w = 1.95
    D = np.zeros_like(a[1,:])
    Res = 100
    k = 0
    while Res>1e-6:
        u_k = np.copy(u_k_new) 
        for j in range(1,n_rows-1): # the number is based on d, which is (N-1)x(N-1)
            # D = np.zeros_like(a[1,:])
            D[1:-1] = d[j-1,:] - e[j,1:-1]*u_k_new[j-1,1:-1] - f[j,1:-1]*u_k_new[j+1,1:-1]
            D[0], D[-1] = u_k_new[j,1], u_k_new[j,-2] # ghost point = ghost point
            u_k_new[j,:] = TDMA(a[j,:],b[j,:],c[j,:],D) # TDMA solve uk each line
            u_k_new[j,1:-1] = u_k[j,1:-1]*(1-w) +u_k_new[j,1:-1]*w # SOR term
        u_k_new = Ghost_BC_p(u_k_new)
        Res = np.max(np.abs(Residual(u_k_new, a, b, c, d, e, f)))/n_rows**2
        # print("Res_p=")
        # print(Res) 
        k+=1
    # print("finish iteration")
    return u_k_new, k, Res




def Solver(u, v, U, V , u_old, v_old, U_old, V_old,  p, data, parameter_field):

    Re = data["Re"]
    dx = data["dx"]
    dt = data["dt"]
    t = data["t"]
    t_0 = data['t_0']

    
    r = dt/(2 * Re * dx**2 )
    c0 = dt/(dx)
    

    rows_cols = np.shape(u)
    c_W_uv , c_P_uv , c_E_uv , c_S_uv , c_N_uv  =\
          coeff_assemble(rows_cols, -r, 1+4*r, -r, -r, -r)
    c_W_uv , c_E_uv , c_S_uv , c_N_uv , c_P_uv = \
          IBM_uv_coeff_change(c_W_uv, c_E_uv, c_S_uv, c_N_uv, c_P_uv, iF, iS)
    
    c_W_p , c_P_p , c_E_p , c_S_p , c_N_p  = \
          coeff_assemble(rows_cols, 1/(dx*dx), -4/(dx*dx), 1/(dx*dx), 1/(dx*dx), 1/(dx*dx))

    c_W_p , c_E_p , c_S_p , c_N_p , c_P_p = \
          IBM_p_coeff_change(c_W_p, c_E_p, c_S_p, c_N_p, c_P_p, iF, iS)
    
    RHS = np.zeros_like(u)
    



    U, V = UV_BC(U, V, u, iF)

    u = Ghost_BC_u(u)
    v = Ghost_BC_v(v)
    # u, v = Ghost_BC(u,v, iF)
    

    u_star, v_star, U_star, V_star = data_copy(u, v, U, V) 
    u_new, v_new, U_new, V_new = data_copy(u, v, U, V) 
    # u_old, v_old, U_old, V_old = data_copy(u, v, U, V)
    


    n_max = int(t/dt)

    Difference = 100

    n= int( t_0/dt)

    while n<n_max and Difference>1e-8 :
    # while n<n_max:
        U, V = UV_BC(U, V, u, iF) # Setup BC
        
        # Step 1: get u_star, v_star, U_star, V_star
        # u = outflow_correct(u, rows_cols[0])
        # u, v = Ghost_BC(u,v, iF)
        u = Ghost_BC_u(u)
        v = Ghost_BC_v(v)
        # u_star, v_star = Ghost_BC(u_star,v_star, iF)
        u_star = Ghost_BC_u(u_star)
        v_star = Ghost_BC_v(v_star)
        

        # get u_star:
        
        d = RHS_uv(u, u_old, U, V, U_old, V_old, c0, r)
        rhs= d[:,:]*iF[1:-1,1:-1]
        u_star ,a, Resu = LineSOR_u(c_W_uv , c_P_uv , c_E_uv, rhs , c_S_uv , c_N_uv, u)
        
        # u_star = outflow_correct(u_star, rows_cols[0])
        u_star = Ghost_BC_u(u_star)
        v_star = Ghost_BC_v(v_star)

        # get v_star:
        d = RHS_uv(v, v_old, U, V, U_old, V_old, c0, r) 
        rhs= d[:,:]*iF[1:-1,1:-1]
        v_star ,a, Resv= LineSOR_v(c_W_uv , c_P_uv , c_E_uv, rhs , c_S_uv , c_N_uv, v)

        u_star = Ghost_BC_u(u_star)
        v_star = Ghost_BC_v(v_star)
        # u_star, v_star = Ghost_BC(u_star,v_star, iF)
        

        # get U_star, V_star
        U_star[1:-1,2:] =  (u_star[1:-1,1:-1] + u_star[1:-1,2:])/2
        V_star[2:,1:-1] =  (v_star[2:,1:-1] + v_star[1:-1,1:-1])/2
        U_star, V_star = UV_BC(U_star, V_star, u_star, iF) # Setup BC

        Diff = np.sum(U_star[1:-1, -1]-U_star[1:-1,1])
        U_star[1:-1, -1] -= 1/(rows_cols[0]-2) * Diff

        
        # Step 2: get New Pressure P
        d = (
            U_star[1:-1, 2:] - U_star[1:-1, 1:-1] +
            V_star[2:, 1:-1] - V_star[1:-1, 1:-1]
        )/(dt*dx)

        # print(d[4,5]-(U_star[5,7]-U_star[5,6]+V_star[6,6] - V_star[5,6])/(dt*dx))

        rhs= d[:,:]*iF[1:-1,1:-1]
        p_new,k, Resp= LineSOR_p(c_W_p , c_P_p , c_E_p, rhs , c_S_p , c_N_p, p)
        # print(p_new-Ghost_BC_p(p_new))
        p_new = Ghost_BC_p(np.copy(p_new))

        # Step 3: get u_new, v_new, U_new , V_new
        # update u_new, 
        u_new[1:-1,1:-1] = u_star[1:-1, 1:-1] - dt/dx*(iF[1:-1,2:]*p_new[1:-1,2:] + iS[1:-1,2:]*p_new[1:-1,1:-1] -
                                                       iF[1:-1,:-2]*p_new[1:-1,:-2] - iS[1:-1,:-2]*p_new[1:-1,1:-1])/2 
        v_new[1:-1,1:-1] = v_star[1:-1, 1:-1] - dt/dx*(iF[2:,1:-1]*p_new[2:,1:-1] + iS[2:,1:-1]*p_new[1:-1,1:-1] -
                                                       iF[:-2,1:-1]*p_new[:-2,1:-1] - iS[:-2,1:-1]*p_new[1:-1,1:-1])/2

        u_new = Ghost_BC_u(u_new); v_new = Ghost_BC_v(v_new)
        # update U_new , V_new
        U_new[1:-1,2:] = U_star[1:-1,2:] - dt/dx*(p_new[1:-1,2:]-p_new[1:-1,1:-1])
        V_new[2:,1:-1] = V_star[2:,1:-1] - dt/dx*(p_new[2:,1:-1]-p_new[1:-1,1:-1])
        U_new, V_new = UV_BC(U_new, V_new, u_new, iF) # Setup BC

        u_old, v_old, U_old, V_old = data_copy(u, v, U, V)
        u, v, U, V = data_copy(u_new, v_new, U_new, V_new)
        # N = rows_cols[1]-2
        Difference = np.sum(np.abs(u[1:-1,1:-1] - u_old[1:-1,1:-1]) + 
                            np.abs(v[1:-1,1:-1] - v_old[1:-1,1:-1]) + 
                            np.abs(p_new[1:-1,1:-1] - p[1:-1,1:-1]) )
      

        p = np.copy(p_new)
        
        div = (U_new[1:-1, 2:] - U_new[1:-1, 1:-1]+V_new[2:, 1:-1] - V_new[1:-1, 1:-1])/dx
        div = np.linalg.norm(div)
        p = p[:,:]*iF[:,:]
        u = Ghost_BC_u(u)
        v = Ghost_BC_v(v)
        p = Ghost_BC_p(p)

        n +=1

        N = 1/data['dx']
        if n%2 == 0 or n==n_max:
            output_save(u,v,p, U, V, u_old,v_old, U_old, V_old, data['Lx'], data['Ly'] ,label=[N,Re,n*dt])

        if n%1000 == 0:
            draw_control(u,v,p, data['Lx'], data['Ly'], label=[N,Re,t, n*dt])
        print("time=%.3f, Difference=%.9f, p iteration times=%d, Res_p=%f, div=%.2f "%((n+1)*dt,Difference,k,Resp, div))
    N = 1/data['dx']
    output_save(u,v,p, U, V, u_old,v_old, U_old, V_old, data['Lx'], data['Ly'] ,label=[N,Re,n*dt])
    draw_control(u,v,p, data['Lx'], data['Ly'], label=[N,Re, n*dt])
    return u, v, p




def output_save(u,v,p, U, V, u_old,v_old, U_old, V_old, len_x, len_y ,label):
    N = label[0]
    Re = label[1]
    t = label[2]
    with open("result/save/output_N%d_Re=%d_t=%.2f_%dx%d.txt" % (int(N), int(Re), t, int(len_x), int(len_y)), 'w') as f:
        f.write(f'For grid size N= {label}:\n')
        f.write('Array u:\n')  # 添加注释
        f.write(','.join(map(str, u)) + '\n\n')  # 写入u数组,并在数组后添加换行以分隔
        f.write('Array v:\n')  # 添加注释
        f.write(','.join(map(str, v)) + '\n\n')  # 写入v数组,并在数组后添加换行以分隔
        f.write('Array p:\n')  # 添加注释
        f.write(','.join(map(str, p)) + '\n')  # 写入p数组

        f.write('Array U:\n')  # 添加注释
        f.write(','.join(map(str, U)) + '\n')  # 写入U数组
        f.write('Array V:\n')  # 添加注释
        f.write(','.join(map(str, V)) + '\n')  # 写入V数组

        f.write('Array u_old:\n')  # 添加注释
        f.write(','.join(map(str, u_old)) + '\n\n')  # 写入u数组,并在数组后添加换行以分隔
        f.write('Array v_old:\n')  # 添加注释
        f.write(','.join(map(str, v_old)) + '\n\n')  # 写入v数组,并在数组后添加换行以分隔

        f.write('Array U_old:\n')  # 添加注释
        f.write(','.join(map(str, U_old)) + '\n\n')  # 写入u数组,并在数组后添加换行以分隔
        f.write('Array V_old:\n')  # 添加注释
        f.write(','.join(map(str, V_old)) + '\n\n')  # 写入v数组,并在数组后添加换行以分隔
        


def draw_control(u,v,p, Lx, Ly, label):
    draw_streamline(u[1:-1,1:-1],v[1:-1,1:-1], Lx, Ly , label)
    draw_velocity_contour(u[1:-1,1:-1],v[1:-1,1:-1], Lx, Ly , label)
    draw_u_contour(u[1:-1,1:-1],v[1:-1,1:-1], Lx, Ly , label)
    draw_p_contour(p[1:-1,1:-1],v[1:-1,1:-1], Lx, Ly , label)





def draw_streamline(u, v, len_x, len_y ,label):
    rows, cols = np.shape(u)  
    X = np.linspace(0, len_x, cols)
    Y = np.linspace(0, len_y, rows)
    L = np.sqrt(len_x**2+len_y**2)
    plt.figure(figsize=(15*len_x/L, 15*len_y/L))
    plt.axis([0,len_x,0,len_y])
    plt.streamplot(X, Y, u, v,  density=1.5, minlength=1, arrowsize=0.5)
    plt.savefig(f'result/stline_N{label[0]}_Re={label[1]}_t={label[2]}_{len_x}x{len_y}.png')
    plt.close()
    # plt.show()


def draw_velocity_contour(u, v, len_x, len_y,label):
    levels=100 
    rows, cols = np.shape(u)
    X, Y = np.linspace(0, len_x, cols), np.linspace(0, len_y, rows)
    velocity_magnitude = np.sqrt(u**2 + v**2)
    L = np.sqrt(len_x**2+len_y**2)
    plt.figure(figsize=(15*len_x/L, 15*len_y/L))
    plt.contourf(X, Y, velocity_magnitude, levels=levels, cmap='jet')
    plt.colorbar(label='Velocity magnitude')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Velocity Magnitude Contour')
    plt.savefig(f'result/uv_N{label[0]}_Re={label[1]}_t={label[2]}_{len_x}x{len_y}.png')
    plt.close()
    # plt.show()


def draw_u_contour(u, v, len_x, len_y ,label):
    levels=100
    rows, cols = np.shape(u)
    X, Y = np.linspace(0, len_x, cols), np.linspace(0, len_y, rows)
    L = np.sqrt(len_x**2+len_y**2)
    plt.figure(figsize=(15*len_x/L, 15*len_y/L))
    plt.contourf(X, Y, u, levels=levels, cmap='jet')
    plt.colorbar(label='u magnitude')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('u Contour')
    plt.savefig(f'result/u_N{label[0]}_Re={label[1]}_t={label[2]}_{len_x}x{len_y}.png')
    plt.close()
    # plt.show()


def draw_p_contour(p, v, len_x, len_y ,label):
    levels=100
    rows, cols = np.shape(p)
    X, Y = np.linspace(0, len_x, cols), np.linspace(0, len_y, rows)
    L = np.sqrt(len_x**2+len_y**2)
    plt.figure(figsize=(15*len_x/L, 15*len_y/L))
    plt.contourf(X, Y, p, levels=levels, cmap='jet')
    plt.colorbar(label='u magnitude')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('u Contour')
    plt.savefig(f'result/p_N{label[0]}_Re={label[1]}_t={label[2]}_{len_x}x{len_y}.png')
    plt.close()
    # plt.show()



def load_2d_arrays_with_numpy(filename):
    # # filename = f'cyl_output_N=32.0, Re=150.0, t=10.0.txt'
    # filename = f'cyl_output_N=32, Re=150, t=50.txt'
    try:
        with open(filename, 'r') as file:
            content = file.read()
    except FileNotFoundError:
        print("文件未找到,请检查文件路径是否正确")
        return None, None, None

    # 辅助函数来提取并构建二维数组
    def extract_2d_array(array_content):
        array_data = []
        # 匹配方括号内的所有内容,包括可能的多行
        matches = re.findall(r'\[(.*?)]', array_content.replace('\n', ''), re.DOTALL)
        for match in matches:
            # 替换掉可能存在的换行符,然后解析数值
            formatted_match = ' '.join(match.split())
            numbers = np.fromstring(formatted_match, sep=' ')
            array_data.append(numbers)
        return np.array(array_data)

    u_data = re.search(r'Array u:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    v_data = re.search(r'Array v:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    U_data = re.search(r'Array U:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    V_data = re.search(r'Array V:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    u_old_data = re.search(r'Array u_old:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    v_old_data = re.search(r'Array v_old:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    U_old_data = re.search(r'Array U_old:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    V_old_data = re.search(r'Array V_old:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    p_data = re.search(r'Array p:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)

    if not (u_data and v_data and p_data):
        print("未能正确提取u, v, 或 p的数据。请检查数组标记和文件格式。")
        return np.array([]), np.array([]), np.array([])

    u = extract_2d_array(u_data.group(1))
    v = extract_2d_array(v_data.group(1))
    p = extract_2d_array(p_data.group(1))
    U = extract_2d_array(U_data.group(1))
    V = extract_2d_array(V_data.group(1))
    u_old = extract_2d_array(u_old_data.group(1))
    v_old = extract_2d_array(v_old_data.group(1))
    U_old = extract_2d_array(U_old_data.group(1))
    V_old = extract_2d_array(V_old_data.group(1))
    return u, v, U, V , u_old, v_old, U_old, V_old,  p


def Recover_result(filename):
    u, v, U, V , u_old, v_old, U_old, V_old,  p = load_2d_arrays_with_numpy(filename)
    pattern = r'output_N(\d+\.\d+)_Re=(\d+\.\d+)_t=(\d+\.\d+)_'
    match = re.search(pattern, filename)
    if match:
        N = float(match.group(1))
        Re = float(match.group(2))
        t_0 = float(match.group(3))
        # print("N =", N)
        # print("Re =", Re)
        # print("t =", t)
    else:
        print("No match found")
    Ly = 1/N* (np.size(u,0)-2)
    Lx = 1/N* (np.size(u,1)-2)
    dy = dx = 1/N
    dt = 0.005
    data_value = np.array([Re, Ly, Lx, dy, dx ,t_0 ,dt])
    keys = ["Re", "Ly", "Lx", "dy", "dx" ,"t_0" ,"dt"]
    data_re = dict(zip(keys, data_value))
    return data_re, u, v, U, V , u_old, v_old, U_old, V_old,  p





def main():
    print(1)

    
    
    

    
    N = 32
    Ly = 4
    Lx = 8

    # N = 30
    # Ly = 1
    # Lx = 2

    dy = dx = 1/N
    Re = 1000

    t = 111
    dt = 0.005
    t_0 = 0

    u_0 = 0

    data_re, u_0, v_0, U_0, V_0 , \
        u_old_0, v_old_0, U_old_0, V_old_0,  p_0 \
            = Recover_result(filename = f'result/save/output_N32.0_Re=1000.0_t=100.0_8.0x4.0222.txt')
    t_0 = data_re['t_0']
    Re = data_re['Re']
    Ly = data_re["Ly"]
    Lx = data_re["Lx"]
    dx = data_re['dx']
    dt = data_re['dt']
    print(Ly)
    print('data_read:',data_re)






    
    global data
    data_value = np.array(\
        [ Re , Ly , Lx , dy , dx , t , dt , t_0 ]\
            )
    keys = \
        ["Re","Ly","Lx","dy","dx","t","dt","t_0"]
    data = dict(zip(keys, data_value))
    print('data_init=',data)

    
    mesh1 =  mesher(data)
    parameter_field = mesh1.parameter_grid()



    global iF, iS
    iF = parameter_field[:,:,2]
    iS = parameter_field[:,:,3]
    

    u, v, U, V ,p = [np.zeros_like(parameter_field[:,:,0]) for _ in range(5)]
    u_old, v_old, U_old, V_old =\
          np.zeros_like(u), np.zeros_like(v), np.zeros_like(U), np.zeros_like(V)
    
    if np.sum(u_0) !=0:
        u, v, U, V , u_old, v_old, U_old, V_old,  p =u_0, v_0, U_0, V_0 , u_old_0, v_old_0, U_old_0, V_old_0,  p_0

    print(u)

    # u,v, p =Solver(u, v, U, V, p, data, parameter_field)
    u,v, p =Solver(u, v, U, V , u_old, v_old, U_old, V_old,  p, data, parameter_field)

    pass



if __name__ == '__main__':
    main()
        
        
















