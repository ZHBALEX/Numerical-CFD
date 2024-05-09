# This solver is a new attempt using same size of u, U... , 
# keep updating ghost point and dont need to change formula

import numpy as np
import math


def general_grid_generator(len_x, Nx, Ny):
    # Creating same size of u, v, U, V
    # where col 0, row 0, row -1, of U not used
    # where row 0, col 0, col -1, of V not used
    d = len_x/Nx
    grid = np.full((Nx+2, Ny+2), 1)
    u,v,U,V = data_copy(grid, grid, grid, grid)
    p = np.copy(grid)
    return u, v, U, V, p

def data_copy(u, v, U, V):
    return np.copy(u), np.copy(v), np.copy(U), np.copy(V)

def Ghost_BC(u,v): # updating ghost point based on BC condition
    # (u_indomain + u_ghost)/2 = Boundary_value
    # u_ghost = 2 * Boundary_value - u_indomain
    u[-1,:] = 2*1 - u[-2,:] # upper BC
    v[-1,:] =  - v[-2,:]

    u[0, :] =  - u[1,:] # bottom BC
    v[0, :] =  - v[1,:]

    u[:, 0] =  - u[:,1] # left BC
    v[:, 0] =  - v[:,1]

    u[:, -1] =  - u[:,-2] # right BC
    v[:, -1] =  - v[:,-2]

    return u, v

def UV_BC(U, V):
    U[:,1] = 0
    U[:,-1] = 0
    V[1,:] = 0
    V[-1,:] = 0

    return U,V



def u_coeff_assemble(input, a_e, b_e, c_e, d_e, e_e, f_e):
    # input elements of abcdef, get line matrix of abcdef
    size = np.size(input)
    a = np.full(size, a_e )
    b = np.full(size, b_e )
    c = np.full(size, c_e )
    d = np.full(size, d_e )
    e = np.full(size, e_e )
    f = np.full(size, f_e )
    return a,b,c,d,e,f


def LineSOR(a,b,c,d,e,f, u):
    # input abcdef, and u, get D
    # then use TDMA get new u

    u_line = u.reshape(-1)
    number = len(u_line)
    D = np.copy(d)

    Res = 100
    while Res>1e-6:
        for I in range(0,number+1):
            D[I] = -e[I-1]*u_line[I-1] -f[I+1]*u_line[I+1] + d[I]
        u_line_new = TDMA(a,b,c,D)
        

        










def Solver(u, v, U, V, p, t, dt):
    a = np.full(len())
    
    # Step 1: get u_star, v_star, U_star, V_star
    u_star, v_star, U_star, V_star = data_copy(u, v, U, V) 
    u_new, v_new, U_new, V_new = data_copy(u, v, U, V) 

    n_max = int(t/dt)
    for n in range(0, n_max):
        # Step 1: get u_star, v_star, U_star, V_star
        # get u_star:



    
    


    
    





def main():
    N = 4
    t = 1

    dt = 0.001

    len_x = 1
    len_y = 1
    Nx = Ny = N
    u, v, U, V, p = general_grid_generator(len_x, Nx, Ny)







    













    pass







if '__name__' == '__main__':
    main()






