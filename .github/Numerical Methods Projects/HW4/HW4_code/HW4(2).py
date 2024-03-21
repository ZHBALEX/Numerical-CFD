
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

def eploting(   ej, egs, SORej1, SORegs1, SORej2, SORegs2,SORej3, SORegs3,  d,  x_len, y_len,   k     ):
    x = [0 for _ in range(0, k+1)]
    for c in range(0,k+1):
        x[c] = c*d
    
    plt.plot(x, ej, color = 'orangered', label = 'Jacobi')
    plt.plot(x, egs, color = 'maroon', label = 'Gauss-Seidel')
    plt.plot(x, SORej1, color = 'darkolivegreen', label = 'Jacobi SOR w=0.5')
    plt.plot(x, SORegs1, color = 'gold', label = 'Gauss-Seidel SOR w=0.5')
    plt.plot(x, SORej2, color = 'lightblue', label = 'Jacobi SOR w=0.7')
    plt.plot(x, SORegs2, color = 'darkslategrey', label = 'Gauss-Seidel SOR w=0.7')   
    plt.plot(x, SORej3, color = 'mediumvioletred', label = 'Jacobi SOR w=1.5')
    plt.plot(x, SORegs3, color = 'thistle', label = 'Gauss-Seidel SOR w=1.5')
    plt.legend()
    plt.yscale('log')
    plt.xlabel('k')
    plt.ylabel('Error')
    
    #plt.ylim(0,5)
    plt.ylim(0,10**(-1))

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

    k = 500

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

    ##FOr w = 0.5 SOR #####################################
    w = 0.5
    r = copy.deepcopy(P)
    e = [0 for _ in range(0,k+1)]
    
    SORej1 = copy.deepcopy(e)
    SORegs1 = copy.deepcopy(e)

    P_Jin = copy.deepcopy(P)
    P_GSin = copy.deepcopy(P)


    for t in range(0,k+1):
        
        
        SORP_Jacobi = SORJacobi(P_Jin, x_len, y_len, w)
        SORP_GS = SORGS(P_GSin, x_len, y_len, w)

        P_Jin = SORP_Jacobi
        P_GSin = SORP_GS

        SORr_J = Res(r,SORP_Jacobi,d,x_len,y_len)
        SORr_GS = Res(r,SORP_GS,d,x_len,y_len)

        SORej1[t] = Error(SORr_J, x_len, y_len)
        SORegs1[t] = Error(SORr_GS, x_len, y_len)

    ##FOr w = 0.7 SOR #####################################
    w = 0.7
    r = copy.deepcopy(P)
    e = [0 for _ in range(0,k+1)]
    
    SORej2 = copy.deepcopy(e)
    SORegs2 = copy.deepcopy(e)

    P_Jin = copy.deepcopy(P)
    P_GSin = copy.deepcopy(P)


    for t in range(0,k+1):
        
        
        SORP_Jacobi = SORJacobi(P_Jin, x_len, y_len, w)
        SORP_GS = SORGS(P_GSin, x_len, y_len, w)

        P_Jin = SORP_Jacobi
        P_GSin = SORP_GS

        SORr_J = Res(r,SORP_Jacobi,d,x_len,y_len)
        SORr_GS = Res(r,SORP_GS,d,x_len,y_len)

        SORej2[t] = Error(SORr_J, x_len, y_len)
        SORegs2[t] = Error(SORr_GS, x_len, y_len)
    
    ## FOr w = 0.5 SOR #####################################
    w = 1.5
    r = copy.deepcopy(P)
    e = [0 for _ in range(0,k+1)]
    
    SORej3 = copy.deepcopy(e)
    SORegs3 = copy.deepcopy(e)

    P_Jin = copy.deepcopy(P)
    P_GSin = copy.deepcopy(P)


    for t in range(0,k+1):
        
        
        SORP_Jacobi = SORJacobi(P_Jin, x_len, y_len, w)
        SORP_GS = SORGS(P_GSin, x_len, y_len, w)

        P_Jin = SORP_Jacobi
        P_GSin = SORP_GS

        SORr_J = Res(r,SORP_Jacobi,d,x_len,y_len)
        SORr_GS = Res(r,SORP_GS,d,x_len,y_len)

        SORej3[t] = Error(SORr_J, x_len, y_len)
        SORegs3[t] = Error(SORr_GS, x_len, y_len)











        '''
        for wc in range(1,8):
            w = wc*0.25
            SORP_Jacobi = SORJacobi(P_Jin, x_len, y_len, w)
            SORP_GS = SORGS(P_GSin, x_len, y_len, w)

            SORP_Jin = SORP_Jacobi
            SORP_GSin = SORP_GS

            SORr_J = Res(r,SORP_Jacobi,d,x_len,y_len)
            SORr_GS = Res(r,SORP_GS,d,x_len,y_len)
            SORej[t] = Error(SORr_J, x_len, y_len)
            SORegs[t] = Error(SORr_GS, x_len, y_len)
        '''


    
    



    
    
    
    eploting(  ej, egs, SORej1, SORegs1, SORej2, SORegs2,SORej3, SORegs3,d,   x_len, y_len,     k   )
    #print(ej)






# main function operator
if __name__ == "__main__":
    main()





























