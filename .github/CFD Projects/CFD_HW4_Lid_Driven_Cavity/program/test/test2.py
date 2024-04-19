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

def LineGS(a,b,c,d,e,f, u): 
    # input abcdef, and u, get D
    # then use TDMA get new u
    # abc is in same size of N_rows of u, ef in same size of N_cols of u
    # ! Input d is LIKE inner field of u, need to cut each line to use !
    u_k = np.copy(u)
    u_k_new = np.copy(u)
    n_rows = np.size(u,0)
    w = 1.5
    D = np.zeros_like(a)
    Res = 100
    while Res>1e-6:
        print(u_k)
        u_k = np.copy(u_k_new) 
        for j in range(0,n_rows-2): # the number is based on d, which is (N-1)x(N-1)
            D[1:-1] = d[j,:] - e[j+1]*u_k_new[j,1:-1] - f[j+2]*u_k[j+2,1:-1]
            D[0], D[-1] = u[j+1,0], u[j+1,-1] # ghost point = ghost point
            u_k_new[j+1,:] = TDMA(a,b,c,D) # TDMA solve uk each line
            # u_k_new[j,:] = u_k[j,:]*(1-w) + TDMA(a,b,c,D)*w # SOR term
        
        Res = np.sum(np.abs(u_k_new - u_k))
        print("Res=")
        print(Res) 
        if Res > np.sum(np.abs(u))*1e+20: 
            print("NOT CONVERGE NOT CONVERGE NOT CONVERGE NOT CONVERGE") 
            quit()
        u_k = np.copy(u_k_new) 
        
    u = np.copy(u_k_new)
    print("finish iteration")
    return u_k_new

N = 8

a = np.full(N, 1)
b = np.full(N, -4)
c = np.full(N, 1)
d = np.full((N-2,N-2), 0)

e = np.full(N, 1)
f = np.full(N, 1)

a[-1] = c[0] = 0
b[0] = b[-1] = 1
f[0]  = e[-1] = 0

u = np.full((N,N),0, dtype='float')

BC = 1
u[0,:] = BC
u[-1,:] = BC
u[:,0] = u[:,-1] = BC

# u[5,6] = 1

result = LineGS(a,b,c,d,e,f, u)

print(result)
