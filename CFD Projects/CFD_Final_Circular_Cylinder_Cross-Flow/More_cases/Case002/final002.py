

import numpy as np
import matplotlib.pyplot as plt
plt.ioff()  # 关闭交互模式
import copy


class mesher:
    def __init__(self, data):
        self.Lx = data["Lx"]
        self.Ly = data["Ly"]
        # self.Re = data["Re"]
        self.dxy = data["dy"]

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

        Field1.flags.writeable = False
        return Field1
    
def grid_generator( y_size,  x_size, grid_size):
        Nx , Ny = int(x_size/grid_size), int(y_size/grid_size)
        Field = np.full((Ny+2, Nx+2, 10), 0 ,dtype = 'float')
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

def IBM_u_outer_BC(c_W, c_E, c_S, c_N, c_P, RHS):
    # West Boundary: u_boundary = 1
    W = 1 
    print(RHS[1:-1,1])
    c_P[1:-1,1] -= c_W[1:-1,1] ; RHS[1:-1,1] -= 2*W*c_W[1:-1,1] ; 
    print(c_W[1:-1,1])
    print(RHS[1:-1,1])
    c_W[1:-1,1] = 0
    # North Boundary: du/dy = 0
    c_P[-2,1:-1] += c_N[-2,1:-1] ;     c_N[-2,1:-1] = 0
    # Sorth Boundary: du/dy = 0
    c_P[1,1:-1] += c_S[1,1:-1];        c_S[1,1:-1] = 0
    # East Boundary: du/dx = 0
    c_P[1:-1,-2] += c_E[1:-1,-2];      c_E[1:-1,-2] = 0
    
    # #North_West:
    # c_P[-2,1] -= (c_W[-2,1] - c_N[-2,1]) ; RHS[-2,1] -= 2*W*c_W[-2,1]; c_W[-2,1], c_N[-2,1] = 0 , 0
    # #South_West:
    # c_P[1,1] -= (c_W[1,1] - c_S[1,1])    ; RHS[1,1] -= 2*W*c_W[1,1]  ; c_W[1,1], c_S[1,1] = 0 , 0
    # #North_East:
    # c_P[-2,-2] += (c_E[-2,-2] + c_N[-2,-2]);                           c_E[-2,-2], c_N[-2,-2] = 0 , 0
    # #South_East:
    # c_P[1,-2] += (c_E[1,-2] + c_S[1,-2]);                              c_E[-2,-2], c_S[1,-2] = 0 , 0
    return c_W, c_E, c_S, c_N, c_P, RHS


def IBM_v_outer_BC(c_W, c_E, c_S, c_N, c_P, RHS):
    # West Boundary:
    W = 0 # West boundary
    c_P[2:-2,1] -= c_W[2:-2,1];   RHS[2:-2,1] -=2*W* c_W[2:-2,1];  c_W[2:-2,1] = 0
    # North Boundary:
    North = 0
    c_P[-2,2:-2] -= c_N[-2,2:-2]; RHS[-2,2:-2] -= 2*North*c_N[-2,2:-2]; c_N[-2,2:-2] = 0
    # Sorth Boundary:
    S = 0
    c_P[1,2:-2] += c_S[1,2:-2];   RHS[1,2:-2] -= 2*S*c_S[1,2:-2];  c_S[1,2:-2] = 0
    # East Boundary:
    c_P[2:-2,-2] += c_E[2:-2,-2];                                   c_E[2:-2,-2] = 0

    # #North_West:
    # c_P[-2,1] -= (c_W[-2,1] + c_N[-2,1]) ;   c_W[-2,1], c_N[-2,1] = 0 , 0
    # #South_West:
    # c_P[1,1] -= (c_W[1,1] + c_S[1,1])    ;   c_W[1,1], c_S[1,1] = 0 , 0
    # #North_East:
    # c_P[-2,-2] += (c_E[-2,-2] - c_N[-2,-2]); c_E[-2,-2], c_N[-2,-2] = 0 , 0
    # #South_East:
    # c_P[1,-2] += (c_E[1,-2] - c_S[1,-2]);    c_E[1,-2], c_S[1,-2] = 0 , 0
    return c_W, c_E, c_S, c_N, c_P, RHS


def IBM_p_outer_BC(c_W, c_E, c_S, c_N, c_P, RHS):
    # West Boundary:
    c_P[2:-2,1] += c_W[2:-2,1];     c_W[2:-2,1] = 0
    # North Boundary:
    c_P[-2,2:-2] += c_N[-2,2:-2];   c_N[-2,2:-2] = 0
    # Sorth Boundary:
    c_P[1,2:-2] += c_S[1,2:-2];     c_S[1,2:-2] = 0
    # East Boundary:
    c_P[2:-2,-2] += c_E[2:-2,-2];   c_E[2:-2,-2] = 0

    # #North_West:
    # c_P[-2,1] += (c_W[-2,1] + c_N[-2,1]) ;   c_W[-2,1], c_N[-2,1] = 0 , 0
    # #South_West:
    # c_P[1,1] += (c_W[1,1] + c_S[1,1])    ;   c_W[1,1], c_S[1,1] = 0 , 0
    # #North_East:
    # c_P[-2,-2] += (c_E[-2,-2] + c_N[-2,-2]); c_E[-2,-2], c_N[-2,-2] = 0 , 0
    # #South_East:
    # c_P[1,-2] += (c_E[1,-2] + c_S[1,-2]);    c_E[1,-2], c_S[1,-2] = 0 , 0
    return c_W, c_E, c_S, c_N, c_P, RHS







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

# def outflow_correct(u, Ny):
#     Diff = np.sum(u[1:-1,-2]-u[1:-1,1])
#     u[:, -2] -= 1/(Ny-2) * Diff
#     return np.copy(u)

# def Ghost_BC(u,v, iF): # updating ghost point based on BC condition
#     # (u_indomain + u_ghost)/2 = Boundary_value
#     # u_ghost = 2 * Boundary_value - u_indomain
#     u[-1,1:-1] =  u[-2,1:-1] # upper BC
#     v[-1,1:-1] =  -v[-2,1:-1]

#     u[0,1:-1] =  u[1,1:-1] # bottom BC
#     v[0, 1:-1] =  -v[1,1:-1]

#     u[1:-1, 0] =  2 - u[1:-1,1] # left BC
#     v[1:-1, 0] =   -v[1:-1,1]

#     u[1:-1, -1] =   u[1:-1,-2] # right BC
#     v[1:-1, -1] =   v[1:-1,-2]
    
#     u = u[:,:]*iF[:,:]
#     v = v[:,:]*iF[:,:]
#     return u, v


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

# def Not_Use_UV_BC(U,V):
#     U[0,:] = U[-1,:] = U[:,0] = V[:,0] = V[:,-1] = V[0,:] = -10000
#     return U,V




def Residual(u, a, b, c, d, e, f):
    return (a[1:-1,1:-1] * u[1:-1,:-2] +
                      b[1:-1,1:-1] * u[1:-1,1:-1] +
                      c[1:-1,1:-1] * u[1:-1,2:] +
                      e[1:-1,1:-1] * u[:-2,1:-1] +
                      f[1:-1,1:-1] * u[2:,1:-1] -
                      d[:,:])

def Residual_a(u, a, b, c, d, e, f):
    return (a[1:-1,1:-1] * u[1:-1,:-2] +
                      b[1:-1,1:-1] * u[1:-1,1:-1] +
                      c[1:-1,1:-1] * u[1:-1,2:] +
                      e[1:-1,1:-1] * u[:-2,1:-1] +
                      f[1:-1,1:-1] * u[2:,1:-1] -
                      d[1:-1,1:-1])



def Math_Matirx_Vector(f, c_W, c_P, c_E, c_S, c_N):
    return (c_W[1:-1,1:-1] * f[1:-1,:-2] +
                      c_P[1:-1,1:-1] * f[1:-1,1:-1] +
                      c_E[1:-1,1:-1] * f[1:-1,2:] +
                      c_S[1:-1,1:-1] * f[:-2,1:-1] +
                      c_N[1:-1,1:-1] * f[2:,1:-1] )

def Math_dot_product(a,b):
    return a[:]*b[:]

def Math_Linag_Solver_BPiCGSTAB(f, c_W, c_E, c_P, c_S, c_N, RHS):
    x_0 = np.copy(f)
    
    r_0 = RHS - Math_Matirx_Vector(x_0, c_W, c_P, c_E, c_S, c_N) # 1. r0 = b-Ax0
    r_0_a = np.copy(r_0) #2.
    rho_0 = Math_dot_product(r_0_a,r_0) #3.
    p_0 = np.copy(r_0) #4.

    p_old = np.copy(r_0)
    rho_old = np.copy(rho_0)
    x_old = np.copy(x_0)
    r_old = np.copy(r_0)
    Residual = 100
    tolerance = 1e-6
    while Residual > tolerance: # 5.
        v = Math_Matirx_Vector(p_old, c_W, c_P, c_E, c_S, c_N) #1.
        
        alpha = rho_old/Math_dot_product(r_0_a,v) #2.
        h = x_old + alpha*p_old
        s = r_old -alpha*v
        if  np.linalg.norm(s)<tolerance: # step 5.
            return h
        
        t = Math_Matirx_Vector(s, c_W, c_P, c_E, c_S, c_N)
        w = Math_dot_product(t,s)/Math_dot_product(t,t)
        x = h + w*s #8.
        r = s-w*t
        if np.linalg.norm(r)<tolerance: # step 5.
            return x
        
        rho = Math_dot_product(r_0, r)
        beta = (rho/rho_old)*(alpha/w)
        p_old = r + beta(p_old - w*v) # 13.
        
        rho_old = rho
        x_old = x
        r_old = r
        Residual = np.linalg.norm(r)

    return x

# def Point_GS(a,b,c,d,e,f,u):
#     u_k = np.copy(u)
#     u_k_new = np.copy(u)
#     n_rows = np.size(u,0)
#     w = 1
#     Res = 100
#     k = 0
#     while Res>1e-8:
#         u_k = np.copy(u_k_new) 
#         u_k_new[1:-1,1:-1] = (-a[1:-1,1:-1] * u_k_new[1:-1,:-2] +
#                               -c[1:-1,1:-1] * u_k_new[1:-1,2:] +
#                               -e[1:-1,1:-1] * u_k_new[:-2,1:-1] +
#                               -f[1:-1,1:-1] * u_k_new[2:,1:-1] +
#                               d[:,:])/b[1:-1,1:-1]
#         u_k_new[1:-1] = u[1:-1]*(1-w) + u_k_new[1:-1]*w # SOR term
        
#         # Res = np.sqrt(np.mean(np.square((Residual(u_k_new, a, b, c, d, e, f)))))
#         Res = np.max(np.abs(Residual(u_k_new, a, b, c, d, e, f)))
#         # print("Res=")
#         # print(Res) 
#         u = np.copy(u_k)
#         k+=1
#     # print("finish iteration")
        
#     return u_k_new, k, Res



# def Point_GS_uv(a,b,c,d,e,f,u):
#     u_k = np.copy(u)
#     u_k_new = np.copy(u)
#     n_rows = np.size(u,0)
#     w = 1
#     Res = 100
#     k = 0
#     while Res>1e-8:
#         u_k = np.copy(u_k_new) 
#         u_k_new[1:-1,1:-1] = (-a[1:-1,1:-1] * u_k_new[1:-1,:-2] +
#                               -c[1:-1,1:-1] * u_k_new[1:-1,2:] +
#                               -e[1:-1,1:-1] * u_k_new[:-2,1:-1] +
#                               -f[1:-1,1:-1] * u_k_new[2:,1:-1] +
#                               d[:,:])/b[1:-1,1:-1]
#         u_k_new[1:-1] = u[1:-1]*(1-w) + u_k_new[1:-1]*w # SOR term
        
#         # Res = np.sqrt(np.mean(np.square((Residual(u_k_new, a, b, c, d, e, f)))))
#         Res = np.max(np.abs(Residual(u_k_new, a, b, c, d, e, f)))
#         # print("Res=")
#         # print(Res) 
#         u = np.copy(u_k)
#         k+=1
#     # print("finish iteration")
        
#     return u_k_new, k, Res




# def LineSOR(a,b,c,d,e,f, u): 
#     # input abcdef, and u, get D
#     # then use TDMA get new u
#     # abc is in same size of N_rows of u, ef in same size of N_cols of u
#     # ! Input d is LIKE inner field of u, need to cut each line to use !
#     u_k = np.copy(u)
#     u_k_new = np.copy(u)
#     n_rows = np.size(u,0)
#     w = 1
#     D = np.zeros_like(a[1,:])
#     Res = 100
#     k = 0
#     while Res>1e-6:
#         u_k = np.copy(u_k_new) 
#         for j in range(1,n_rows-1): # the number is based on d, which is (N-1)x(N-1)
#             # D = np.zeros_like(a[1,:])
#             D[1:-1] = d[j-1,:] - e[j,1:-1]*u_k_new[j-1,1:-1] - f[j,1:-1]*u_k[j+1,1:-1]
#             D[0], D[-1] = u[j,0], u[j,-1] # ghost point = ghost point
#             u_k_new[j,:] = TDMA(a[j,:],b[j,:],c[j,:],D) # TDMA solve uk each line
#             # u_k_new[j,:] = u[j,:]*(1-w) + TDMA(a[j,:],b[j,:],c[j,:],D)*w # SOR term
#             u
#         Res = np.linalg.rms(Residual(u_k_new, a, b, c, d, e, f))
#         # print("Res=")
#         # print(Res) 
#         u = np.copy(u_k)
#         k+=1
#     # print("finish iteration")
#     return u_k_new, k, Res


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





def Solver(u, v, U, V, p, data, parameter_field):

    Re = data["Re"]
    dx = data["dx"]
    dt = data["dt"]
    t = data["t"]

    # global iF, iS

    # iF = parameter_field[:,:,2]
    # iS = parameter_field[:,:,3]

    
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
    u_old, v_old, U_old, V_old = data_copy(u, v, U, V)
    


    n_max = int(t/dt)

    Difference = 100

    n=0

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


        # print((p_new[1,3]+p_new[1,1]+p_new[0,2]+p_new[2,2]-4*p_new[1,2]) 
        #       - dx/dt * (U_star[1,3]-U_star[1,2]+V_star[2,2] - V_star[1,2]))
        
        # print(c_P_p[5,6])
        # print((p_new[5,6]+p_new[5,5]+p_new[6,6]+p_new[4,6]-4*p_new[5,6]) 
        #       - dx/dt * (U_star[5,7]-U_star[5,6]+V_star[6,6] - V_star[5,6]))
        

        # print((p_new[-2,:-2]+p_new[-2,2:]+p_new[-1,1:-1]+p_new[-3,1:-1]-3*p_new[-2,1:-1]) 
        #       - dx/dt * (U_star[-2,2:]-U_star[-2,1:-1]+ V_star[-1,1:-1] - V_star[-3,1:-1]))

        # Step 3: get u_new, v_new, U_new , V_new
        # update u_new, 
        u_new[1:-1,1:-1] = u_star[1:-1, 1:-1] - dt/dx*(iF[1:-1,2:]*p_new[1:-1,2:] + iS[1:-1,2:]*p_new[1:-1,1:-1] -
                                                       iF[1:-1,:-2]*p_new[1:-1,:-2] - iS[1:-1,:-2]*p_new[1:-1,1:-1])/2 
        v_new[1:-1,1:-1] = v_star[1:-1, 1:-1] - dt/dx*(iF[2:,1:-1]*p_new[2:,1:-1] + iS[2:,1:-1]*p_new[1:-1,1:-1] -
                                                       iF[:-2,1:-1]*p_new[:-2,1:-1] - iS[:-2,1:-1]*p_new[1:-1,1:-1])/2
        
        # u_new = outflow_correct(u_new, rows_cols[0])
        # u_new, v_new = Ghost_BC(u_new,v_new, iF)
        u_new = Ghost_BC_u(u_new)
        v_new = Ghost_BC_v(v_new)

        
        # update U_new , V_new
        U_new[1:-1,2:] = U_star[1:-1,2:] - dt/dx*(p_new[1:-1,2:]-p_new[1:-1,1:-1])
        V_new[2:,1:-1] = V_star[2:,1:-1] - dt/dx*(p_new[2:,1:-1]-p_new[1:-1,1:-1])
        U_new, V_new = UV_BC(U_new, V_new, u_new, iF) # Setup BC

        # Diff = np.sum(U_new[1:-1, -1]-U_new[1:-1,1])
        # U_new[1:-1, -1] -= 1/(rows_cols[0]-2) * Diff

        # print((p_new[-2,:-2]+p_new[-2,2:]+p_new[-1,1:-1]+p_new[-3,1:-1]-3*p_new[-2,1:-1]) 
        #       - dx/dt * (U_new[-2,2:]-U_new[-2,1:-1]+ 0 - V_new[-3,1:-1]))


        

        # print('new done')

        u_old, v_old, U_old, V_old = data_copy(u, v, U, V)
        u, v, U, V = data_copy(u_new, v_new, U_new, V_new)
        # N = rows_cols[1]-2
        Difference = np.sum(np.abs(u[1:-1,1:-1] - u_old[1:-1,1:-1]) + 
                            np.abs(v[1:-1,1:-1] - v_old[1:-1,1:-1]) + 
                            np.abs(p_new[1:-1,1:-1] - p[1:-1,1:-1]) )
      

        p = np.copy(p_new)
        
        div = (U_new[1:-1, 2:] - U_new[1:-1, 1:-1]+V_new[2:, 1:-1] - V_new[1:-1, 1:-1])/dx
        # print(div[-1,:])

        # X, Y = np.linspace(0, 1, rows_cols[0]-2), np.linspace(0, 1, rows_cols[1]-2)
        
        # plt.contourf(Y, X, div, levels=10, cmap='jet')
        # plt.show()

        # X, Y = np.linspace(0, 1, rows_cols[0]-2), np.linspace(0, 1, rows_cols[1]-2)
        
        # plt.contourf(Y, X, div, levels=10, cmap='jet')
        # plt.show()
        
        div = np.linalg.norm(div)
        
        
        p = p[:,:]*iF[:,:]
        # u = outflow_correct(u, rows_cols[0])
        u = Ghost_BC_u(u)
        v = Ghost_BC_v(v)
        p = Ghost_BC_p(p)

        n +=1

        N = 1/data['dx']
        if n%1000 == 0:
            draw_control(u,v,p, data['Lx'], data['Ly'], label=[N,Re,t, n*dt])
        

        
        print("time=%.3f, Difference=%.9f, p iteration times=%d, Res_p=%f, div=%.2f "%((n+1)*dt,Difference,k,Resp, div))
        
        # div = (U_new[1:-1, 2:] - U_new[1:-1, 1:-1]+V_new[2:, 1:-1] - V_new[1:-1, 1:-1])/dx
        # X, Y = np.linspace(0, 1, rows_cols[0]-2), np.linspace(0, 1, rows_cols[1]-2)
        # levels = np.arange(-10,10,0.1)
        # plt.contourf(Y, X, div, levels=levels, cmap='jet')
        # plt.colorbar(label='div magnitude')
        # plt.show()
    return u, v, p


def draw_control(u,v,p, Lx, Ly, label):
    draw_streamline(u[1:-1,1:-1],v[1:-1,1:-1], Lx, Ly , label)
    draw_velocity_contour(u[1:-1,1:-1],v[1:-1,1:-1], Lx, Ly , label)
    draw_u_contour(u[1:-1,1:-1],v[1:-1,1:-1], Lx, Ly , label)
    draw_p_contour(p[1:-1,1:-1],v[1:-1,1:-1], Lx, Ly , label)
    save(u,v,p,Lx, Ly ,label)





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



def save(u,v,p,len_x, len_y ,label):
    with open(f'result/output_N{label[0]}_Re={label[1]}_t={label[2]}_{len_x}x{len_y}.txt', 'w') as f:
        f.write(f'For grid size N= {label}:\n')
        f.write('Array u:\n')  # 添加注释
        f.write(','.join(map(str, u)) + '\n\n')  # 写入u数组,并在数组后添加换行以分隔
        f.write('Array v:\n')  # 添加注释
        f.write(','.join(map(str, v)) + '\n\n')  # 写入v数组,并在数组后添加换行以分隔
        f.write('Array p:\n')  # 添加注释
        f.write(','.join(map(str, p)) + '\n')  # 写入p数组







def main():
    print(1)
    
    N = 32
    Ly = 4
    Lx = 8

    # N = 30
    # Ly = 1
    # Lx = 2

    dy = dx = 1/N
    Re = 150

    t = 100
    dt = 0.01

    global data
    data_value = np.array([Re, Ly, Lx, dy, dx ,t ,dt])
    keys = ["Re", "Ly", "Lx", "dy", "dx" ,"t" ,"dt"]
    data = dict(zip(keys, data_value))
    print(data)

    mesh1 =  mesher(data)
    parameter_field = mesh1.parameter_grid()

    global iF, iS

    iF = parameter_field[:,:,2]
    iS = parameter_field[:,:,3]

#   # test ifluid
#     plt.imshow(parameter_field[:,:,3], cmap = 'gray')
#     plt.show()  


    u, v, U, V ,p = [np.zeros_like(parameter_field[:,:,0]) for _ in range(5)]


    u,v, p =Solver(u, v, U, V, p, data, parameter_field)

    # print(u)
    # print(v)
    # print(p)

    save(u,v,p,Lx, Ly ,label=[N,Re,t])

    draw_streamline(u[1:-1,1:-1],v[1:-1,1:-1], Lx, Ly , label=[N,Re,t])
    draw_velocity_contour(u[1:-1,1:-1],v[1:-1,1:-1], Lx, Ly , label=[N,Re,t])
    draw_u_contour(u[1:-1,1:-1],v[1:-1,1:-1], Lx, Ly , label=[N,Re,t])
    draw_p_contour(p[1:-1,1:-1],v[1:-1,1:-1], Lx, Ly , label=[N,Re,t])


    
    pass



if __name__ == '__main__':
    main()
        
        
















