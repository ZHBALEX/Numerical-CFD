def TDMA(a, b, c, D):  # TDMA solver

		import copy
		Do = copy.deepcopy(D)
		ao = copy.deepcopy(a)
		bo = copy.deepcopy(b)
		co = copy.deepcopy(c)
		
		Lx = len(D) 
		co[1] = co[1] / bo[1]
		Do[1] = Do[1] / bo[1]
		for nx in range(2, Lx-1):
			bo[nx] -= co[nx - 1] * ao[nx]
			Do[nx] -= Do[nx - 1] * ao[nx]
			co[nx] = co[nx] / bo[nx]
			Do[nx] = Do[nx] / bo[nx]
		for nx in reversed(range(1, Lx-1)):
			Do[nx] -= co[nx] * Do[nx + 1]
		return Do


def solver_exact(Lx, Lt, m, dt, dx, Te):  # exact solution solver
    import math

    Te = [[0 for _ in range(0, Lx + 1)] for _ in range(0, Lt + 1)]
    for nt in range(0, Lt+1):
        for nx in range(0, Lx+1):
            Te[nt][nx] = math.exp(-(m**2) * nt * dt) * math.sin(m * nx * dx)
    return Te


def solver_FE(Lt, Lx, r, Tfe):
    for nt in range(0, Lt):
        for nx in range(1, Lx):
            Tfe[nt + 1][nx] = (
                r * (Tfe[nt][nx + 1] + Tfe[nt][nx - 1]) + (1 - 2 * r) * Tfe[nt][nx]
            )
    return Tfe


def solver_BE(Lx, Lt, r, Tbe):
    a_be = [-r 				for _ in range(0, Lx )]
    b_be = [1 + 2 * r for _ in range(0, Lx )]
    c_be = [-r 				for _ in range(0, Lx )]
    for nt in range(1, Lt+1):
        Dbe = Tbe[nt - 1].copy()
        Tbe[nt] = TDMA(a_be, b_be, c_be, Dbe)
    return Tbe


def solver_CN(Lx, Lt, r, Tcn):
    a_cn = [-r 						for _ in range(0, Lx )]
    b_cn = [(2 * (1 + r)) for _ in range(0, Lx )]
    c_cn = [-r 						for _ in range(0, Lx )]
    a1 = [r 						for _ in range(0, Lx )]
    b1 = [(2 * (1 - r)) for _ in range(0, Lx )]
    c1 = [r 						for _ in range(0, Lx )]

    for nt in range(1, Lt+1):
        Dcn = Tcn[nt -1].copy()
        
        Dcn[1] = b1[1] * Tcn[nt -1][1] + c1[1] * Tcn[nt -1][2]
        for nx in range(2, Lx - 1):
            Dcn[nx] = a1[nx] * Tcn[nt -1][nx - 1] + b1[nx] * Tcn[nt -1][nx] + c1[nx] * Tcn[nt -1][nx + 1]
        Dcn[Lx-1] = a1[Lx-1]*Tcn[nt -1][Lx-2] + b1[Lx-1]*Tcn[nt -1][Lx-1]
        Tcn[nt] = TDMA(a_cn, b_cn, c_cn, Dcn)
    return Tcn



def error(T1,Te,Lt,Lx):
    E = [[0 for _ in range(0, Lx + 1)] for _ in range(0, Lt + 1)]

    for nt in range(1, Lt + 1):
        for nx in range(1, Lx + 1):
            E[nt][nx] = T1[nt][nx] - Te[nt][nx]
    return E


def Serror(E,T1,Te,Lx,Lt      ):
    import math
    sE = [0 for _ in range(0, Lt + 1)]
    E = error(T1,Te,Lt,Lx)
    for nt in range(0,Lt+1):
        for nx in range(0,Lx+1):
             sE[nt] += E[nt][nx]**2

    for nt in range(1,Lt+1):     
        sE[nt] = math.sqrt(sE[nt]/Lx)
    return sE






'''
def Serror(E,T1,Te,Lx,Lt      ):
    import math
    sE = [0 for _ in range(0, Lt + 1)]
    E = error(T1,Te,Lt,Lx)
    for nt in range(0,Lt+1):
        for nx in range(0,Lx+1):
             sE[nt] += E[nt][nx]**2

    for nt in range(1,Lt+1):     
        sE[nt] = math.sqrt(sE[nt]/Lx)
    return sE
'''




def plot(Tfe, Tbe, Tcn, Te, Lx, Lt, dx, dt,t_max,m,r):
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib as mpl
    import math

		
    Efe = error(Tfe,Te,Lt,Lx)
    Ebe = error(Tbe,Te,Lt,Lx)
    Ecn = error(Tcn,Te,Lt,Lx)

    sEfe = Serror(Efe,Tfe,Te,Lx,Lt      )
    sEbe = Serror(Ebe,Tbe,Te,Lx,Lt      )
    sEcn = Serror(Ecn,Tcn,Te,Lx,Lt      )
    

    
    x = [0]
    t = [0]
    for i in range(1, Lx + 1):
        x.append(i * dx)
    for i in range(1, Lt + 1):
        t.append(i * dt)

    plt.figure(figsize=(120,6))  
    plt.subplot(1,3,1)

    time = Lt


    plt.plot(t, Tfe[time][4], color="green",label = "FE")
    plt.plot(t, Tbe[time][4], color="black",label = "BE")
    plt.plot(t, Te[time][4], color="red",label = "Exact")
    plt.plot(t, Tcn[time][4], color="blue",label = "CN")

    plt.legend()
    plt.xlabel('x')
    plt.ylabel('T')
    
    
    
    plt.title("The solutions at t=%f, m=%f, r=%f"%(Lt*dt,m,r))


    
    

    plt.show()
    return


######################################################


def main():
    import math
    import copy

    pi = math.pi

    t_max = 0.5
    x_max = 2 * pi
    r = 1/3
    dx = 2 * pi / 20
    dt = r * dx**2
    Lx = int(x_max / dx)
    Lt = int(t_max / dt)

    m = 2
    T = [[0 for _ in range(0, Lx + 1)] for _ in range(0, Lt + 1)]

    for nt in range(1, Lt + 1):
        T[nt][0] = 0
        T[nt][Lx] = 0
    for nx in range(0, Lx +1):
        T[0][nx] = math.sin(m * nx * dx)

    Tfe = copy.deepcopy(T)
    Tbe = copy.deepcopy(T)
    Tcn = copy.deepcopy(T)
    Te = copy.deepcopy(T)

    Tfe = solver_FE(Lt, Lx, r, Tfe)
    Tbe = solver_BE(Lx, Lt, r, Tbe)
    Tcn = solver_CN(Lx, Lt, r, Tcn)
    Te = solver_exact(Lx, Lt, m, dt, dx, Te)


    plot(Tfe, Tbe, Tcn, Te, Lx, Lt, dx, dt,t_max,m,r)
    
    print(Lt*dt)

    pass


if __name__ == "__main__":
    main()





    






