# This solver is a new attempt using same size of u, U... , 
# keep updating ghost point and dont need to change formula

import numpy as np
import math
import matplotlib.pyplot as plt


def TDMA(a, b, c, d):
    n = len(d)
    c_prime = np.zeros(n-1)
    d_prime = np.zeros(n)
    x = np.zeros(n)

    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]
    
    for i in range(1, n-1):
        temp = b[i] - a[i] * c_prime[i-1]
        c_prime[i] = c[i] / temp
        d_prime[i] = (d[i] - a[i] * d_prime[i-1]) / temp
    
    d_prime[n-1] = (d[n-1] - a[n-1] * d_prime[n-2]) / (b[n-1] - a[n-1] * c_prime[n-2])
    
    x[n-1] = d_prime[n-1]
    for i in range(n-2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i+1]
    
    return x



def general_grid_generator(len_x, Nx, Ny):
    # Creating same size of u, v, U, V
    # where col 0, row 0, row -1, of U not used
    # where row 0, col 0, col -1, of V not used
    initial_Value = 0.01

    dx = len_x/Nx
    grid = np.full((Nx+2, Ny+2), initial_Value)
    u,v,U,V = data_copy(grid, grid, grid, grid)
    p = np.copy(grid)
    return u, v, U, V, p, dx

def data_copy(u, v, U, V):
    return np.copy(u), np.copy(v), np.copy(U), np.copy(V)

def Ghost_BC(u,v): # updating ghost point based on BC condition
    # (u_indomain + u_ghost)/2 = Boundary_value
    # u_ghost = 2 * Boundary_value - u_indomain
    u[-1,:] = 2*1 - u[-2,:] # upper BC
    v[-1,:] =  - v[-2,:]

    u[0, :] =  - u[1,:] # bottom BC
    v[0, :] =  - v[1,:]

    u[1:-1, 0] =  - u[1:-1,1] # left BC
    v[1:-1, 0] =  - v[1:-1,1]

    u[1:-1, -1] =  - u[1:-1,-2] # right BC
    v[1:-1, -1] =  - v[1:-1,-2]
    return u, v


def Ghost_BC_p(p): # updating ghost point for p based on BC condition
    p[-1,:] =   p[-2,:] # upper BC
    p[-1,:] =   p[-2,:]

    p[0, :] =   p[1,:] # bottom BC
    p[0, :] =   p[1,:]

    p[1:-1, 0] =   p[1:-1,1] # left BC
    p[1:-1, 0] =   p[1:-1,1]

    p[1:-1, -1] =   p[1:-1,-2] # right BC
    p[1:-1, -1] =   p[1:-1,-2]
    return p



def UV_BC(U, V):
    U[1:-1,1] = 0
    U[1:-1,-1] = 0
    V[1,1:-1] = 0
    V[-1,1:-1] = 0
    return U,V



def coeff_assemble(rows_cols, a_e, b_e, c_e, e_e, f_e):
    rows, cols = rows_cols[0], rows_cols[1]
    # input elements of abcdef, get line matrix of abcdef
    a = np.full(cols, a_e ) # abc in same length of Ncols of u
    b = np.full(cols, b_e )
    c = np.full(cols, c_e )

    e = np.full(rows, e_e ) # ef in same length of Nrows of u
    f = np.full(rows, f_e )
    a[0]= a[-1] = c[0] = c[-1] = 0
    b[0] = b[-1] = 1
    e[-1] =e[0] = f[-1]= f[0]=0
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

def Not_Use_UV_BC(U,V):
    U[0,:] = U[-1,:] = U[:,0] = V[:,0] = V[:,-1] = V[0,:] = -10000
    return U,V


def LineSOR(a, b, c, d, e, f, u):
    n_rows, n_cols = u.shape
    u_new = np.copy(u)
    Res = np.inf
    tolerance = 1e-1
    D = np.zeros(n_cols)
    w = 1.4
    while Res > tolerance:
        u_old = np.copy(u_new)
        for j in range(1, n_rows-1):
            D[1:-1] = d[j-1, :] - e[j-2]*u_old[j, :-2] - f[j-1]*u_old[j+1, 2:]
            D[0] = u[j, 0]
            D[-1] = u[j, -1]
            u_new[j, :] = u[j,:]*(1-w)+w*TDMA(a, b, c, D)

        Res = np.linalg.norm(u_new - u_old)
        if Res > 1e20:
            raise Exception("Not converging.")

    return u_new






def Solver(u, v, U, V, p, t, dt, dx, Re):
    r = dt/(2*Re*dx**2)
    c0 = dt/(dx)

    a_e, c_e, e_e, f_e = -r, -r, -r, -r
    b_e = 1+4*r
    rows_cols = np.shape(u)

    U, V = UV_BC(U, V)
    u, v = Ghost_BC(u,v)

    U, V = Not_Use_UV_BC(U, V) 


    u_star, v_star, U_star, V_star = data_copy(u, v, U, V) 
    u_new, v_new, U_new, V_new = data_copy(u, v, U, V) 
    u_old, v_old, U_old, V_old = data_copy(u, v, U, V)


    n_max = int(t/dt)

    Big_Res = 100

    n=0

    while Big_Res>1e-8:
        U, V = UV_BC(U, V) # Setup BC

        # Step 1: get u_star, v_star, U_star, V_star
        u, v = Ghost_BC(u,v)
        u_star, v_star = Ghost_BC(u_star,v_star)

        # get u_star:
        a,b,c,e,f  = coeff_assemble(rows_cols, a_e, b_e, c_e, e_e, f_e)
        d = RHS_uv(u, u_old, U, V, U_old, V_old, c0, r)
        u_star = LineSOR(a,b,c,d,e,f, u)

        # print("u* done")
        # print("DDDDDDDDDDDDDDDDDDDDDDDDDDDDD")

        u_star, v_star = Ghost_BC(u_star,v_star)

        # get v_star:
        a,b,c,e,f  = coeff_assemble(rows_cols, a_e, b_e, c_e, e_e, f_e)
        d = RHS_uv(v, v_old, U, V, U_old, V_old, c0, r) 
        v_star = LineSOR(a,b,c,d,e,f, v)

        # print("v* done")
        # print("DDDDDDDDDDDDDDDDDDDDDDDDDDDDD")


        u_star, v_star = Ghost_BC(u_star,v_star)

        # get U_star, V_star
        U_star[1:-1,2:] =  (u_star[1:-1,1:-1] + u_star[1:-1,2:])/2
        V_star[2:,1:-1] =  (v_star[2:,1:-1] + v_star[1:-1,1:-1])/2

        U_star, V_star = UV_BC(U_star, V_star) # Setup BC


        # Step 2: get New Pressure P
        a,b,c,e,f = coeff_assemble(rows_cols, 1, -4, 1, 1, 1)
        d = dx/dt * (
            U_star[1:-1, 2:] - U_star[1:-1, 1:-1] +
            V_star[2:, 1:-1] - V_star[1:-1, 1:-1]
        )
        p_new = LineSOR(a,b,c,d,e,f, p)
        p_new = Ghost_BC_p(p_new)

        # print("p done")
        # print("DDDDDDDDDDDDDDDDDDDDDDDDDDDDD")



        # Step 3: get u_new, v_new, U_new , V_new
        # update u_new, v_new
        u_new[1:-1,1:-1] = u_star[1:-1, 1:-1] - dt/dx*(p_new[1:-1,2:] - p_new[1:-1,:-2])/2 
        v_new[1:-1,1:-1] = v_star[1:-1, 1:-1] - dt/dx*(p_new[2:,1:-1] - p_new[:-2,1:-1])/2 
        u_new, v_new = Ghost_BC(u_new,v_new)
        # update U_new , V_new
        U_new[1:-1,2:] = U_star[1:-1,2:] - dt/dx*(p_new[1:-1,2:]-p_new[1:-1,1:-1])
        V_new[2:,1:-1] = V_star[2:,1:-1] - dt/dx*(p_new[2:,1:-1]-p_new[1:-1,1:-1])
        U_new, V_new = UV_BC(U_new, V_new) # Setup BC

        # print('new done')

        u_old, v_old, U_old, V_old = data_copy(u, v, U, V)
        u, v, U, V = data_copy(u_new, v_new, U_new, V_new)

        Big_Res = np.sum(np.abs(u - u_old)+np.abs(v - v_old)+np.abs(p_new - p))

        p = np.copy(p_new)

        # print(Big_Res)
        n +=1
        print("time=%.3f, Res=%f"%((n+1)*dt,Big_Res))
    
    return u, v, p

def draw(u, v, len_x, len_y):
    rows, cols = np.shape(u)  
    X = np.linspace(0, len_x, cols)
    Y = np.linspace(0, len_y, rows)
    plt.streamplot(X, Y, u, v)
    plt.show()

    
    









def main():
    N = 4
    t = 0.1


    dt = 0.001

    len_x = 1
    len_y = 1
    Nx = Ny = N
    Re = 100

    u, v, U, V, p, dx = general_grid_generator(len_x, Nx, Ny)

    u,v, p =Solver(u, v, U, V, p, t, dt, dx, Re)

    print(u)


    draw(u[1:-1,1:-1],v[1:-1,1:-1], len_x, len_y)




















if __name__ == '__main__':
    main()






