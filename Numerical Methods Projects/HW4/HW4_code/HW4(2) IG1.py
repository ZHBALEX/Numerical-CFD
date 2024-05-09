
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



def SORJacobi(P_inputSOR, x_len, y_len, w): #Jacobi SOR method Pin-->Pout
    P_oldSOR = copy.deepcopy(P_inputSOR)
    P_newSOR = copy.deepcopy(P_inputSOR)
    for j in range(1,y_len-1):
        for i in range(1,x_len-1): 
            P_newSOR[j][i] = \
                1/4 * (  P_oldSOR[j][i-1] + P_oldSOR[j][i+1] + P_oldSOR[j-1][i] + P_oldSOR[j+1][i] ) *w + (1-w)*P_oldSOR[j][i]
    return P_newSOR


def SORGS(P_inputSOR, x_len, y_len,w): #Gauss Sadiel SOR method Pin-->Pout
    PgSOR = copy.deepcopy(P_inputSOR)
    for j in range(1,y_len-1):
        for i in range(1,x_len-1): 
            PgSOR[j][i] = \
                1/4 * (  PgSOR[j][i-1] + P_inputSOR[j][i+1] + PgSOR[j-1][i] + P_inputSOR[j+1][i] )*w + (1-w)*P_inputSOR[j][i]
    return PgSOR





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
            e += abs(r[j][i])
    return e


import matplotlib.pyplot as plt
#import numpy as np
import matplotlib as mpl

def nploting(  nj, ngs ):
    w = [0 for _ in range(0, 20)]
    for c in range(1,20):
        w[c] = c*0.1
    
    plt.plot(w, nj, color = 'orangered', label = 'Jacobi SOR')
    plt.plot(w, ngs, color = 'maroon', label = 'Gauss-Seidel SOR')

    plt.legend()
    plt.xlabel('w')
    plt.ylabel('K number of Iteration to get 10^(-3)')
    #plt.ylim(0,5)
    #plt.ylim(0,10**(-10))





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

    k = 1000

    #import Boundary Condition (BC)
    P = [[0 for _ in range(0,x_len+1)] for _ in range(0, y_len+1)]

    for i in range(0,x_len):
        sin = math.sin
        P[0][i] = sin(2*i*dx) + sin(5*i*dx) + sin(7*i*dx)
        P[y_len-1][i] = 0
    for j in range(0,y_len):
        P[j][0] = 0
        P[j][x_len-1] = 0

    #import Initial Guess (IG)
    for j in range(1,y_len-1):
        for i in range(1,x_len-1):
            P[j][i] = j*i*d*d
    
    #P_input = copy.deepcopy(P)
    r = copy.deepcopy(P)
    e = [0 for _ in range(0,k+1)]
    ej = copy.deepcopy(e)
    egs = copy.deepcopy(e)

    P_Jin = copy.deepcopy(P)
    P_GSin = copy.deepcopy(P)


    en = 10**(-3)

    #calculate iteration number of Jacobi
    for t in range(0,k+1):
        
        P_Jacobi = Jacobi(P_Jin, x_len, y_len)

        P_Jin = P_Jacobi

        r_J = Res(r,P_Jacobi,d,x_len,y_len)
 
        ej = Error(r_J, x_len, y_len)

        if ej <= en:
            nj = t
            break
        if t == k:
            nj= 200
    
    #calculate interation number of Gauss-Seidel
    for t in range(0,k+1):
        
        P_GS = GS(P_GSin, x_len, y_len)

        P_GSin = P_GS


        r_GS = Res(r,P_GS,d,x_len,y_len)

        egs = Error(r_GS, x_len, y_len)
        if egs <= en:
            ngs = t
            break
        if t == k:
            ngs = 200

    #calculate iteration times of different w
    nj = [0 for _ in range(0,20)]
    ngs = [0 for _ in range(0,20)]

    for wc in range(1,20):
        w = wc*0.1
        r = copy.deepcopy(P)

        P_Jin = copy.deepcopy(P)
        P_GSin = copy.deepcopy(P)


        for t in range(0,k+1):
            
            
            SORP_Jacobi = SORJacobi(P_Jin, x_len, y_len, w)

            P_Jin = SORP_Jacobi

            SORr_J = Res(r,SORP_Jacobi,d,x_len,y_len)

            SORej = Error(SORr_J, x_len, y_len)
            if SORej <=en:
                nj[wc] = t
                break
            if t == k:
                nj[wc] = None
        
        for t in range(0,k+1):
            
            
            SORP_GS = SORGS(P_GSin, x_len, y_len, w)

            P_GSin = SORP_GS

            SORr_GS = Res(r,SORP_GS,d,x_len,y_len)

            SORegs = Error(SORr_GS, x_len, y_len)
            if SORegs <=en:
                ngs[wc] = t
                break
            if t == k:
                ngs[wc] = None


   
    nploting( nj , ngs)
    #print(SORegs)
    
    #print(ej)






# main function operator
if __name__ == "__main__":
    main()





























