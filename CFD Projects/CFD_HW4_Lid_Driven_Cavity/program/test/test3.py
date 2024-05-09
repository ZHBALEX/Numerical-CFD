import numpy as np

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

def LineGS(a, b, c, d, e, f, u):
    n_rows, n_cols = u.shape
    u_new = np.copy(u)
    Res = np.inf
    tolerance = 1e-6

    while Res > tolerance:
        u_old = np.copy(u_new)
        for j in range(1, n_rows-1):
            D = np.zeros(n_cols)
            D[1:-1] = d[j-1, :] - e[j-1]*u_old[j, 1:-1] - f[j-1]*u_old[j+1, 1:-1]
            D[0] = u[j, 0]
            D[-1] = u[j, -1]
            u_new[j, :] = TDMA(a, b, c, D)

        Res = np.linalg.norm(u_new - u_old)
        if Res > 1e20:
            raise Exception("Not converging.")
        print(Res)

    return u_new

# Example usage with your specified parameters
N = 8
a = np.full(N, 1)
b = np.full(N, -4)
c = np.full(N, 1)
d = np.zeros((N-2, N-2))
e = np.full(N-2, 1)
f = np.full(N-2, 1)

# Adjusting the ends for the boundary conditions
a[-1], c[0] = 0, 0
b[0], b[-1] = 1, 1

u = np.full((N, N), 0, dtype='float')
BC = 0
u[0, :], u[-1, :], u[:, 0], u[:, -1] = BC, BC, BC, BC

u[5,5] = -100000000
result = LineGS(a, b, c, d, e, f, u)
print(result)
