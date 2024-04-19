# This solver is a new attempt using same size of u, U... , 
# keep updating ghost point and dont need to change formula

import numpy as np
import math

import copy


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
    dx = len_x/Nx
    grid = np.full((Nx+2, Ny+2), 1)
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

def UV_BC(U, V):
    U[:,1] = 0
    U[:,-1] = 0
    V[1,:] = 0
    V[-1,:] = 0
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
    diffusion = r*(f_E + f_W + f_N + f_S -4*f_P)
    d = -3/2*convection_new + 1/2*convection_old + diffusion

    return d



def LineSOR(a,b,c,d,e,f, u): 
    # input abcdef, and u, get D
    # then use TDMA get new u
    # abc is in same size of N_rows of u, ef in same size of N_cols of u
    # ! Input d is LIKE inner field of u, need to cut each line to use !
    u_k = np.copy(u)
    u_k_new = np.copy(u)
    n_rows = np.size(u,0)
    w = 1
    D = np.zeros_like(a)
    Res = 100
    while Res>1e-6:
        u_k = np.copy(u_k_new)
        for j in range(1,n_rows-2):
            D[1:-1] = d[j,:] - e[j-1]*u_k[j-1,1:-1] - f[j+1]*u_k[j+1,1:-1]
            D[0], D[-1] = u_k[j,0], u_k[j,-1] 
            u_k_new[j,:] = TDMA(a,b,c,D) # TDMA solve uk each line
            # u_k_new[j,:] = u_k[j,:]*(1-w) + TDMA(a,b,c,D)*w # SOR term

        # Calculate Residual
        # for j in range(1,n_rows-2):
        #     Res += np.sum(np.abs(
        #         a[1:-1]*u_k[j,2:]+b[1:-1]*u_k[j,1:-1]+c[1:-1]*u_k[j,:-2]+ 
        #         e[j-1]*u[j,1:-1] + f[j+1]*u[j,1:-1] -d[j,:]
        #     ))
        Res = np.sum(np.abs(u_k_new - u_k))/100

        print("Res=")
        print(Res) 
        if Res > 10:
            print("aaaaaaaaaa") 
        # u = np.copy(u_k)
    return u_k_new



def Solver(u, v, U, V, p, t, dt, dx, Re):
    r = dt/(2*Re*dx**2)
    c0 = dt/(dx)

    a_e, c_e, e_e, f_e = -r, -r, -r, -r
    b_e = 1+4*r
    rows_cols = np.size(u,0),np.size(u,1)

    u_star, v_star, U_star, V_star = data_copy(u, v, U, V) 
    u_new, v_new, U_new, V_new = data_copy(u, v, U, V) 
    u_old, v_old, U_old, V_old = data_copy(u, v, U, V)

    n_max = int(t/dt)

    for n in range(0, n_max):
        U, V = UV_BC(U, V) # Setup BC

        # Step 1: get u_star, v_star, U_star, V_star
        u, v = Ghost_BC(u,v)
        u_star, v_star = Ghost_BC(u_star,v_star)

        # get u_star:
        a,b,c,e,f  = coeff_assemble(rows_cols, a_e, b_e, c_e, e_e, f_e)
        d = RHS_uv(u, u_old, U, V, U_old, V_old, c0, r)
        u_star = LineSOR(a,b,c,d,e,f, u)

        u_star, v_star = Ghost_BC(u_star,v_star)

        # get v_star:
        a,b,c,e,f  = coeff_assemble(rows_cols, a_e, b_e, c_e, e_e, f_e)
        d = RHS_uv(v, v_old, U, V, U_old, V_old, c0, r) 
        v_star = LineSOR(a,b,c,d,e,f, v)

        u_star, v_star = Ghost_BC(u_star,v_star)

        # get U_star, V_star
        U_star[1:-1,2:-2] =  (u_star[1:-1,2:-2] + u_star[1:-1,3:-1])/2
        V_star[2:-2,1:-1] =  (v_star[2:-2,1:-1] + v_star[3:-1,1:-1])/2

        U_star, V_star = UV_BC(U_star, V_star) # Setup BC


        # Step 2: get New Pressure P
        a,b,c,e,f = coeff_assemble(rows_cols, 1, -4, 1, 1, 1)
        d = dx/dt * (
            U_star[1:-1, 2:] - U_star[1:-1, 1:-1] +
            V_star[2:, 1:-1] - V_star[1:-1, 1:-1]
        )
        p_new = LineSOR(a,b,c,d,e,f, p)


        # Step 3: get u_new, v_new, U_new , V_new
        # update u_new, v_new
        u_new[1:-1,1:-1] = u_star[1:-1, 1:-1] - dt/dx*(p_new[1:-1,2:] - p_new[1:-1,:-2])/2 
        v_new[1:-1,1:-1] = v_star[1:-1, 1:-1] - dt/dx*(p_new[2:,1:-1] - p_new[:-2,1:-1])/2 
        u_new, v_new = Ghost_BC(u_new,v_new)
        # update U_new , V_new
        U_new[1:-1,2:-1] = U_star[1:-1,2:-1] - dt/dx*(p_new[1:-1,2:-1]-p_new[1:-1,1:-2])
        V_new[2:-1,1:-1] = V_star[2:-1,1:-1] - dt/dx*(p_new[2:-1,1:-1]-p_new[1:-2,1:-1])
        U_new, V_new = UV_BC(U_new, V_new) # Setup BC

        u_old = np.copy(u)
        v_old = np.copy(v)
        U_old = np.copy(U)
        V_old = np.copy(V)
        u = np.copy(u_new)
        v = np.copy(v_new)
        U = np.copy(U_new)
        V = np.copy(V_new)
        p = np.copy(p_new)

        print(n)
    
    return u, v, p









def main():
    N = 16
    t = 0.006


    dt = 0.001

    len_x = 1
    len_y = 1
    Nx = Ny = N
    Re = 100

    u, v, U, V, p, dx = general_grid_generator(len_x, Nx, Ny)

    u,v, p =Solver(u, v, U, V, p, t, dt, dx, Re)


    print(u)







    













    pass







if __name__ == '__main__':
    main()






