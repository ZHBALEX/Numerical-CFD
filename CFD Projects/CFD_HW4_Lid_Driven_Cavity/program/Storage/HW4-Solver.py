

import numpy as np
import math



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


def grid_generator(x_len, y_len, N):
    # create cell domin: cell_field
    # cell boundary domain: cell_bU_field, 
    cell_field = np.full((N,N), 1, dtype = 'int')
    cell_bU_field = np.full((N,N+1), 1, dtype = 'int')
    cell_bV_field = np.full((N+1,N), 1, dtype = 'int')
    return cell_field, cell_bU_field, cell_bV_field
    
def init(u, v, U, V):
    u[:,:],v[:,:],U[:,:],V[:,:] = 0, 0, 0, 0
    return u,v,U,V

def BC(U,V):
    U[0,:], V[0,:] = 0, 0 # bottom boundary
    U[-1,:], V[-1,:] = 1, 0 # top boundary
    U[:,0], V[:,0] = 0, 0 # left boundary
    U[:,-1], V[:,-1] = 0, 0 # right boundary
    return U, V


def LineSOR( u,  U, V, u_ex, U_ex, V_ex,r,C):
    Res = 100
    u_star = np.copy(u)
    rows = np.size(u,0)
    cols = np.size(u,1)
    w = 1.7

    a = np.full(cols, -r/2, dtype = 'float' )
    b = np.full(cols,1+r,dtype = 'float'  )
    c = a
    d = np.full(cols,1+r,dtype = 'float'  )
    
    while Res > 1e-6:
        for j in range(1, rows-1):
            for i in range(1,cols-1): # calculate d
                d2 = 3/2 * C * (U[j,i+1]*u[j,i+1]-U[j,i-1]*u[j,i-1])/2 \
                    + 3/2 * C * (V[j+1,i]*u[j+1,i]-V[j-1,i]*u[j-1,i])/2
                
                d3 = 1/2 * C * (\
                    U_ex[j,i+1]*u_ex[j,i+1]\
                        -U[j,i-1]*u_ex[j,i-1]\
                        )/2 \
                    + 1/2 * C * (\
                        V_ex[j,i]*u_ex[j+1,i]\
                            -V_ex[j-1,i]*u_ex[j-1,i]\
                                )/2
                
                d4 = 1/2 * r * (-4*u[j, i] + u[j+1, i] + u[j-1, i] + u[j, i+1] + u[j, i-1])
                d[i] = d2 + d3 + d4 - r/2 * (u_star[j, i+1] + u_star[j, i-1] -2*u_star[j, i])
            u_star[j,:] = TDMA(a,b,c,d)

        u_star[:,:] = (1-w)*u[:,:] + w*u_star[:,:]
        
        Res = 0

        # Calculate Residual
        for j in range(1,rows-1):
            for i in range(1,cols-1):
                d2 = 3/2 * C * (U[j,i+1]*u_star[j,i+1]-U[j,i-1]*u_star[j,i-1])/2 \
                    + 3/2 * C * (V[j+1,i]*u_star[j+1,i]-V[j-1,i]*u_star[j-1,i])/2
                
                d3 = 1/2 * C * (U_ex[j,i+1]*u_ex[j,i+1]-U_ex[j,i-1]*u_ex[j,i-1])/2 \
                    + 1/2 * C * (V_ex[j+1,i]*u_ex[j+1,i]-V_ex[j-1,i]*u_ex[j-1,i])/2
                d4 = 1/2 * r * (-4*u_star[j, i] + u_star[j+1, i] + u_star[j-1, i] + u_star[j, i+1] + u_star[j, i-1])
                Res += abs(u_star[j,i]-r/2*(-4*u_star[j,i] + u_star[j+1,i] + u_star[j-1,i] + u_star[j,i+1] + u_star[j,i-1])-d2-d3-d4)
        print(Res)
    return u_star   




def pLineSOR(p, U_star, V_star, d, dt ):
        rows = np.size(p,0)
        cols = np.size(p,1)
        k = d/dt
        w = 1.7

        Res = 100
        p_star = np.copy(p)

        a = np.full(cols, 1, dtype = 'float' )
        b = np.full(cols,-2,dtype = 'float'  )
        c = a
        d = np.full(cols,1,dtype = 'float'  )

        
        while Res > 1e-6:
            for j in range(1, rows-1):
                for i in range(1,cols-1): # calculate d
                    d1 = k*(U_star[j,i+1]-U_star[j,i]+V_star[j+1,i]-V_star[j,i])
                    d2 = -(p[j,i+1]+p[j,i-1]-2*p[j,i])
                    d[i] = d1+d2
                p_star[j,:] = TDMA(a,b,c,d)
            
            p_star = (1-w)*p + w*p_star
            
            Res = 0

            # Calculate Residual
            for j in range(1,rows-1):
                for i in range(1,cols-1):
                    d1 = k*(U_star[j,i+1]-U_star[j,i]+V_star[j+1,i]-V_star[j,i])
                    d2 = -(p_star[j,i+1]+p_star[j,i-1]-2*p_star[j,i])
                    Res += abs((-4*p_star[j,i] + p_star[j+1,i] + p_star[j-1,i] + p_star[j,i+1] + p_star[j,i-1])-d1)
            print(Res)
        return p_star 
        
    

def solver(Re, u, v, p,  U, V, t, d, dt):
    t_steps = int(t/dt)
    
    u_ex, v_ex = np.copy(u), np.copy(v)
    U_ex, V_ex = np.copy(U), np.copy(V)
    
    for n in range(0, t_steps):
        U, V = BC(U,V)
        

        u_star, v_star, U_star, V_star = \
            np.copy(u), np.copy(v), np.copy(U), np.copy(V)
        
        u_new, v_new, U_new, V_new = \
            np.copy(u), np.copy(v), np.copy(U), np.copy(V)
        
        # step 1
        r = dt/Re/d**2
        c = dt/d

        
    
        u_star = LineSOR( u,  U, V, u_ex, U_ex, V_ex,r, c)
        v_star = LineSOR( v,  U, V, v_ex, U_ex, V_ex,r, c)

        for j in range(0,np.size(U,0)):
            for i in range(1, np.size(U,1)-1):
                U_star[j,i] = ( u_star[j,i-1] + u_star[j,i] )/2
        
        for j in range(1,np.size(V,0)-1):
            for i in range(0, np.size(V,1)):
                V_star[j,i] = ( v_star[j-1,i] + v_star[j,i] )/2
        
        # step 2
        
        p_new = pLineSOR(p, U_star, V_star, d, dt )

        
        
        
        # step 3
        u_new = np.copy(u_star)
        for j in range(1,np.size(u,0)-1):
            for i in range(1, np.size(u,1)-1):
                u_new[j, i] = u_star[j,i] - dt/d*(p[j,i+1]-p[j,i-1])/2
                v_new[j, i] = v_star[j,i] - dt/d*(p[j+1,i]-p[j-1,i])/2
        
        for j in range(0,np.size(U,0)):
            for i in range(1, np.size(U,1)-1):
                U_new[j,i] = U_star[j,i] - dt/d * (p[j,i]-p[j,i-1])

        for j in range(1,np.size(V,0)-1):
            for i in range(0, np.size(V,1)):
                V_new[j,i] = V_star[j,i] - dt/d * (p[j,i]-p[j-1,i])

        u_ex, v_ex, U_ex, V_ex = u, v, U, V
        u, v, U, V = u_new, v_new, U_new, V_new

        print(n)

    return u,v,U,V





def main():
    x_len = 1
    y_len = 1
    N = 16 # grid cell number, Np = Nx * Ny
    Re = 1
    
    dt = 0.001
    d = x_len/N

    t = 1
    
    cell_field, U_field, V_field = grid_generator(x_len, y_len, N)
    
    u, v = np.copy(cell_field), np.copy(cell_field)
    U, V = np.copy(U_field), np.copy(V_field)
    u,v, U, V = init(u, v, U, V)
    U, V = BC(U, V)
    p = np.copy(u)

    u,v,U,V = solver(Re, u, v, p,  U, V, t, d, dt)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    print("YOU WING")
    
    
    
    
    


if __name__ == '__main__':
    main()