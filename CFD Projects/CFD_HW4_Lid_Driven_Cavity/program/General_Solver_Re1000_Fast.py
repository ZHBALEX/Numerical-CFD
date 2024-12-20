# This solver is a new attempt using same size of u, U... , 
# keep updating ghost point and dont need to change formula

import numpy as np
import math
import matplotlib.pyplot as plt


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
    a = np.full((rows, cols), a_e ) # abc in same length of Ncols of u
    b = np.full((rows, cols), b_e )
    c = np.full((rows, cols), c_e )

    e = np.full((rows, cols), e_e ) # ef in same length of Nrows of u
    f = np.full((rows, cols), f_e )
    a[:,0]= a[:,-1] = c[:,0] = c[:,-1] = 0
    b[:,0] = b[:,-1] = 1
    e[:,-1] =e[:,0] = f[:,-1]= f[:,0]=0
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




def Residual(u, a, b, c, d, e, f):
    Res = np.zeros_like(u)
    Res[1:-1,1:-1] = (a[1:-1,1:-1] * u[1:-1,:-2] +
                      b[1:-1,1:-1] * u[1:-1,1:-1] +
                      c[1:-1,1:-1] * u[1:-1,2:] +
                      e[1:-1,1:-1] * u[:-2,1:-1] +
                      f[1:-1,1:-1] * u[2:,1:-1] -
                      d[:,:])
    # print(Res)
    return Res



from numba import jit

@jit(nopython=True)
def Point_GS(a, b, c, d, e, f, u, w=1.9, tol=1e-6):
    rows, cols = u.shape
    u_new = np.empty_like(u)
    Res = np.inf
    while Res > tol:
        u_new[:] = u  # 创建一个工作数组，避免直接修改 u
        for i in range(1, rows-1):
            for j in range(1, cols-1):
                # 这里严格按照 SOR 迭代公式执行更新
                temp = (d[i, j] - a[i, j] * u[i, j-1] - c[i, j] * u[i, j+1] -
                        e[i, j] * u[i-1, j] - f[i, j] * u[i+1, j]) / b[i, j]
                u_new[i, j] = u[i, j] * (1-w) + w * temp

        Res = np.linalg.norm(u_new - u)
        if Res > 1e+20:
            raise Exception("NOT CONVERGE")
        u[:] = u_new  # 一次迭代结束后更新 u
    return u




def LineSOR(a,b,c,d,e,f, u): 
    # input abcdef, and u, get D
    # then use TDMA get new u
    # abc is in same size of N_rows of u, ef in same size of N_cols of u
    # ! Input d is LIKE inner field of u, need to cut each line to use !
    u_k = np.copy(u)
    u_k_new = np.copy(u)
    n_rows = np.size(u,0)
    w = 1
    D = np.zeros_like(a[1,:])
    Res = 100
    while Res>1e-6:
        u_k = np.copy(u_k_new) 
        for j in range(1,n_rows-1): # the number is based on d, which is (N-1)x(N-1)
            # D = np.zeros_like(a[1,:])
            D[1:-1] = d[j-1,:] - e[j,1:-1]*u_k_new[j-1,1:-1] - f[j,1:-1]*u_k[j+1,1:-1]
            D[0], D[-1] = u[j,0], u[j,-1] # ghost point = ghost point
            u_k_new[j,:] = TDMA(a[j,:],b[j,:],c[j,:],D) # TDMA solve uk each line
            # u_k_new[j,:] = u[j,:]*(1-w) + TDMA(a[j,:],b[j,:],c[j,:],D)*w # SOR term
        
        Res = np.linalg.norm(Residual(u_k_new, a, b, c, d, e, f))
        # print("Res=")
        # print(Res) 
        if Res > 1e+20: 
            print("NOT CONVERGE NOT CONVERGE NOT CONVERGE NOT CONVERGE") 
            quit()
        u = np.copy(u_k)
    # print("finish iteration")
    return u_k_new






def Solver(u, v, U, V, p, t, dt, dx, Re):
    r = dt/(2*Re*dx**2)
    c0 = dt/(dx)
    

    rows_cols = np.shape(u)
    a_uv,b_uv,c_uv,e_uv,f_uv  = coeff_assemble(rows_cols, -r, 1+4*r, -r, -r, -r)
    a_p,b_p,c_p,e_p,f_p  = coeff_assemble(rows_cols, 1, -4, 1, 1, 1)

    U, V = UV_BC(U, V)
    u, v = Ghost_BC(u,v)

    U, V = Not_Use_UV_BC(U, V) 


    u_star, v_star, U_star, V_star = data_copy(u, v, U, V) 
    u_new, v_new, U_new, V_new = data_copy(u, v, U, V) 
    u_old, v_old, U_old, V_old = data_copy(u, v, U, V)


    # n_max = int(t/dt)

    Difference = 100

    n=0

    while Difference>1e-6:
        U, V = UV_BC(U, V) # Setup BC

        # Step 1: get u_star, v_star, U_star, V_star
        u, v = Ghost_BC(u,v)
        u_star, v_star = Ghost_BC(u_star,v_star)

        # get u_star:
        
        d = RHS_uv(u, u_old, U, V, U_old, V_old, c0, r)
        u_star = LineSOR(a_uv,b_uv,c_uv, d , e_uv,f_uv, u)

        # print("u* done")
        # print("DDDDDDDDDDDDDDDDDDDDDDDDDDDDD")

        u_star, v_star = Ghost_BC(u_star,v_star)

        # get v_star:
        a,b,c,e,f  = coeff_assemble(rows_cols,  -r, 1+4*r, -r, -r, -r)
        d = RHS_uv(v, v_old, U, V, U_old, V_old, c0, r) 
        v_star = LineSOR(a_uv,b_uv,c_uv, d , e_uv,f_uv, v)

        # print("v* done")
        # print("DDDDDDDDDDDDDDDDDDDDDDDDDDDDD")


        u_star, v_star = Ghost_BC(u_star,v_star)

        # get U_star, V_star
        U_star[1:-1,2:] =  (u_star[1:-1,1:-1] + u_star[1:-1,2:])/2
        V_star[2:,1:-1] =  (v_star[2:,1:-1] + v_star[1:-1,1:-1])/2

        U_star, V_star = UV_BC(U_star, V_star) # Setup BC


        # Step 2: get New Pressure P
        d = dx/dt * (
            U_star[1:-1, 2:] - U_star[1:-1, 1:-1] +
            V_star[2:, 1:-1] - V_star[1:-1, 1:-1]
        )
        p_new = Point_GS(a_p,b_p,c_p,d, e_p,f_p, p)
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

        Difference = np.sum(np.abs(u - u_old)+np.abs(v - v_old)+np.abs(p_new - p))

        p = np.copy(p_new)

        # print(Big_Res)
        n +=1
        print("time=%.3f, Difference=%.9f"%((n+1)*dt,Difference))
    
    return u, v, p





def draw_streamline(u, v, len_x, len_y ,label):
    rows, cols = np.shape(u)  
    X = np.linspace(0, len_x, cols)
    Y = np.linspace(0, len_y, rows)
    plt.figure(figsize=(12, 12))
    plt.streamplot(X, Y, u, v)
    plt.savefig(f'Re1000_streamline_{label}.png')
    # plt.show()


def draw_velocity_contour(u, v, len_x, len_y,label):
    levels=100 
    rows, cols = np.shape(u)
    X, Y = np.linspace(0, len_x, cols), np.linspace(0, len_y, rows)
    velocity_magnitude = np.sqrt(u**2 + v**2)
    plt.figure(figsize=(12, 12))
    plt.contourf(X, Y, velocity_magnitude, levels=levels, cmap='jet')
    plt.colorbar(label='Velocity magnitude')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Velocity Magnitude Contour')
    plt.savefig(f'Re1000_velocity_contour_{label}.png')
    # plt.show()


def draw_u_contour(u, v, len_x, len_y ,label):
    levels=100
    rows, cols = np.shape(u)
    X, Y = np.linspace(0, len_x, cols), np.linspace(0, len_y, rows)
    plt.figure(figsize=(12, 12))
    plt.contourf(X, Y, u, levels=levels, cmap='jet')
    plt.colorbar(label='u magnitude')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('u Contour')
    plt.savefig(f'Re1000_u_contour_{label}.png')
    # plt.show()



def save(u,v,p,label):
    with open(f'Re1000_output_N={label}.txt', 'w') as f:
        f.write(f'For grid size N= {label}:\n')
        f.write('Array u:\n')  # 添加注释
        f.write(','.join(map(str, u)) + '\n\n')  # 写入u数组，并在数组后添加换行以分隔
        f.write('Array v:\n')  # 添加注释
        f.write(','.join(map(str, v)) + '\n\n')  # 写入v数组，并在数组后添加换行以分隔
        f.write('Array p:\n')  # 添加注释
        f.write(','.join(map(str, p)) + '\n')  # 写入p数组

    
    









def main():
    N = 64
    t = 20



    dt = 0.001

    len_x = 1
    len_y = 1
    Nx = Ny = N
    Re = 1000

    u, v, U, V, p, dx = general_grid_generator(len_x, Nx, Ny)

    u,v, p =Solver(u, v, U, V, p, t, dt, dx, Re)

    print(u)
    print(v)
    print(p)

    save(u,v,p,N)
    

    
    draw_streamline(u[1:-1,1:-1],v[1:-1,1:-1], len_x, len_y , N)
    draw_velocity_contour(u[1:-1,1:-1],v[1:-1,1:-1], len_x, len_y , N)
    draw_u_contour(u[1:-1,1:-1],v[1:-1,1:-1], len_x, len_y , N)




















if __name__ == '__main__':
    main()





