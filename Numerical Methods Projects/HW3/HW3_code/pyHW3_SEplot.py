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
    sE = 0
    E = error(T1,Te,Lt,Lx)
    for nx in range(0,Lx+1):
        sE += E[Lt][nx]**2
        
    
    sE = math.sqrt(sE/Lx)
    return sE




def plot(SEfe, SEbe, SEcn, Te, Lx, Lt, dx_list, dt,t_max,m,r):
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib as mpl
    import math

    SEfe = SEfe[1:]
    SEbe = SEbe[1:]
    SEcn = SEcn[1:]
    dx_list = dx_list[1:]
    ref_x = np.array(dx_list)
    ref_y = 16**(-1)*ref_x**2 
    ref_y2 = 52**(-1)*ref_x**2 
    

    
    x = [0]
    t = [0]
    for i in range(1, Lx + 1):
        x.append(i * dx_list)
    for i in range(1, Lt + 1):
        t.append(i * dt)

    plt.figure(figsize=(10,6))  

    
    
    


    plt.loglog(dx_list, SEfe, '.',color="green",label="Forward Euler")
    plt.loglog(dx_list, SEbe, '.', color="black",label="Backward Euler")
    plt.loglog(dx_list, SEcn, '.',color="blue",label="C-N Method")
    plt.loglog(ref_x, ref_y, color="red", linestyle="--", label="Reference Line (Slope=2)")
    plt.loglog(ref_x, ref_y2, color="red", linestyle="--", label="Reference Line (Slope=2)")
    plt.legend()
    plt.xlabel('dx')
    plt.ylabel('S Error')

    plt.title("Spatial errors")
    plt.grid()
    
    
    

    plt.show()
    return




######################################################


def main():
    import math
    import copy

    pi = math.pi

    t_max = 1
    x_max = 2 * pi
    r = 1/3

    m = 2


    N = []

    for i in range(12,20,1):
         N.append(i)

    dx = [ 0 for _ in range(0,len(N))]
    dt = [ 0 for _ in range(0,len(N))]

    SEfe = [0]
    SEbe = [0]
    SEcn = [0]

    

    for i in range(1,len(dx)):
         
         dx[i] = x_max/N[i]
         dt[i] = r * dx[i]**2
         dx1 = dx[i]
         dt1 = dt[i]
         Lx = int(x_max / dx1)
         Lt = int(t_max / dt1)
         T = [[0 for _ in range(0, Lx + 1)] for _ in range(0, Lt + 1)]
         for nt in range(1, Lt + 1):
            T[nt][0] = 0
            T[nt][Lx] = 0
         for nx in range(0, Lx +1):
            T[0][nx] = math.sin(m * nx * dx1)

         Tfe = copy.deepcopy(T)
         Tbe = copy.deepcopy(T)
         Tcn = copy.deepcopy(T)
         Te = copy.deepcopy(T)
         

         Tfe = solver_FE(Lt, Lx, r, Tfe)
         Tbe = solver_BE(Lx, Lt, r, Tbe)
         Tcn = solver_CN(Lx, Lt, r, Tcn)
         Te = solver_exact(Lx, Lt, m, dt1, dx1, Te)

         Efe = error(Tfe,Te,Lt,Lx)
         Ebe = error(Tbe,Te,Lt,Lx)
         Ecn = error(Tcn,Te,Lt,Lx)

         sEfe = Serror(Efe,Tfe,Te,Lx,Lt      )
         sEbe = Serror(Ebe,Tbe,Te,Lx,Lt      )
         sEcn = Serror(Ecn,Tcn,Te,Lx,Lt      )
         SEfe.append(sEfe)
         SEbe.append(sEbe)
         SEcn.append(sEcn)
    dx_list = dx.copy()
    plot(SEfe, SEbe, SEcn,Te, Lx, Lt, dx_list, dt1,t_max,m,r)



         


    pass


if __name__ == "__main__":
    main()
