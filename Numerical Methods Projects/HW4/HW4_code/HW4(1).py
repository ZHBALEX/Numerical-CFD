
import copy #for Python, need copy op to avoid 

def Jacobi(P_input, x_len, y_len): #Jacobi method Pin-->Pout
    P_old = copy.deepcopy(P_input)
    P_new = copy.deepcopy(P_input)
    for j in range(1,y_len-1):
        for i in range(1,x_len-1): 
            P_new[j][i] = \
                1/4 * (  P_old[j][i-1] + P_old[j][i+1] + P_old[j-1][i] + P_old[j+1][i] )
    return P_new

def GS(P_input, x_len, y_len): #Gauss Sadiel method Pin-->Pout
    Pg = copy.deepcopy(P_input)
    for j in range(1,y_len-1):
        for i in range(1,x_len-1): 
            Pg[j][i] = \
                1/4 * (  Pg[j][i-1] + P_input[j][i+1] + Pg[j-1][i] + P_input[j+1][i] )
    return Pg





def Res(r,Pin,d,x_len,y_len): #calcuate residual (Pin, source)-->Rout
    P = copy.deepcopy(Pin)
    rc = copy.deepcopy(r)
    for j in range(1,y_len-1):
        for i in range(1,x_len-1): 
            rc[j][i] = \
                ( P[j][i-1] + P[j][i+1] + P[j-1][i] + P[j+1][i] - 4*P[j][i] )/(d**2)
    return rc

def Error(r_in, x_len, y_len): #calculate error Rin-->e out
    r = copy.deepcopy(r_in)
    e = 0
    for j in range(1,y_len-1):
        for i in range(1,x_len-1): 
            e += abs(r[j][i])**2
    return e


import matplotlib.pyplot as plt
#import numpy as np
import matplotlib as mpl

def eploting(   ej, egs, ej1, egs1,ej2, egs2,   d,  x_len, y_len,   k     ):
    x = [0 for _ in range(0, k+1)]
    for c in range(0,k+1):
        x[c] = c*d
    plt.plot(x, ej, color = 'orangered', label = 'Jacobi IG=0')
    plt.plot(x, egs, color = 'maroon', label = 'Gauss-Seidel IG=0')
    plt.plot(x, ej1, color = 'lawngreen', label = 'Jacobi IG=xiyj')
    plt.plot(x, egs1, color = 'olivedrab', label = 'Gauss-Seidel IG=xiyj')
    plt.plot(x, ej2, color = 'deepskyblue', label = 'Jacobi IG=random')
    plt.plot(x, egs2, color = 'teal', label = 'Gauss-Seidel IG=random')
    plt.legend()

    #plt.yscale('log')
    plt.xlabel('k')
    plt.ylabel('Error')
    #plt.ylim(0,1)
    plt.ylim(0,10**(-6))

    plt.title("k=%d"%k)





    plt.show()




########################## main function below ########


def main():
    import math


    #import grid data
    pi = math.pi
    x_max = 2*pi
    y_max = 2*pi
    d = 2*pi/20
    dx = dy = d
    x_len = int(x_max/dx+1)
    y_len = int(y_max/dy+1)

    # x_len: the number of x points
    # the first x point is x[0], the last x point is x[x_len -1]

    k = 2000

    #import Boundary Condition (BC)
    P = [[0 for _ in range(0,x_len+1)] for _ in range(0, y_len+1)]

    for i in range(0,x_len):
        sin = math.sin
        P[0][i] = sin(2*i*dx) + sin(5*i*dx) + sin(7*i*dx)
        P[y_len-1][i] = 0
    for j in range(0,y_len):
        P[j][0] = 0
        P[j][x_len-1] = 0

    #import Initial Guess (IG = 0)
    IG = 0 #change IG for different guess
    for j in range(1,y_len-1):
        for i in range(1,x_len-1):
            P[j][i] = IG
    
    #P_input = copy.deepcopy(P)
    r = copy.deepcopy(P)
    e = [0 for _ in range(0,k+1)]
    ej = copy.deepcopy(e)
    egs = copy.deepcopy(e)

    P_Jin = copy.deepcopy(P)
    P_GSin = copy.deepcopy(P)
    

    for t in range(0,k+1):
        
        P_Jacobi = Jacobi(P_Jin, x_len, y_len)
        P_GS = GS(P_GSin, x_len, y_len)

        P_Jin = P_Jacobi
        P_GSin = P_GS

        r_J = Res(r,P_Jacobi,d,x_len,y_len)
        r_GS = Res(r,P_GS,d,x_len,y_len)
        ej[t] = Error(r_J, x_len, y_len)
        egs[t] = Error(r_GS, x_len, y_len)
    
    ######################## Ploting initial guess is xiyj
     #change IG for different guess
    for j in range(1,y_len-1):
        for i in range(1,x_len-1):
            P[j][i] = j*i*d*d
    
    #P_input = copy.deepcopy(P)
    r1 = copy.deepcopy(P)
    e1 = [0 for _ in range(0,k+1)]
    ej1 = copy.deepcopy(e1)
    egs1 = copy.deepcopy(e1)

    P_Jin1 = copy.deepcopy(P)
    P_GSin1 = copy.deepcopy(P)
    

    for t in range(0,k+1):
        
        P_Jacobi1 = Jacobi(P_Jin1, x_len, y_len)
        P_GS1 = GS(P_GSin1, x_len, y_len)

        P_Jin1 = P_Jacobi1
        P_GSin1 = P_GS1

        r_J1 = Res(r1,P_Jacobi1,d,x_len,y_len)
        r_GS1 = Res(r1,P_GS1,d,x_len,y_len)
        ej1[t] = Error(r_J1, x_len, y_len)
        egs1[t] = Error(r_GS1, x_len, y_len)

    ######################## Ploting initial guess is xiyj
     #change IG for different guess
    import random
    for j in range(1,y_len-1):
        for i in range(1,x_len-1):
            P[j][i] = random.randint(-1,1)
    
    #P_input = copy.deepcopy(P)
    r2 = copy.deepcopy(P)
    e2 = [0 for _ in range(0,k+1)]
    ej2 = copy.deepcopy(e2)
    egs2 = copy.deepcopy(e2)

    P_Jin2 = copy.deepcopy(P)
    P_GSin2 = copy.deepcopy(P)
    

    for t in range(0,k+1):
        
        P_Jacobi2 = Jacobi(P_Jin2, x_len, y_len)
        P_GS2 = GS(P_GSin2, x_len, y_len)

        P_Jin2 = P_Jacobi2
        P_GSin2 = P_GS2

        r_J2 = Res(r2,P_Jacobi2,d,x_len,y_len)
        r_GS2 = Res(r2,P_GS2,d,x_len,y_len)
        ej2[t] = Error(r_J2, x_len, y_len)
        egs2[t] = Error(r_GS2, x_len, y_len)


    
    
    
    eploting(  ej, egs, ej1, egs1,ej2, egs2, d,   x_len, y_len,     k   )
    #print(ej)






# main function operator
if __name__ == "__main__":
    main()





























