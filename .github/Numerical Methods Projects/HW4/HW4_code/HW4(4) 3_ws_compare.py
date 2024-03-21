
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

def eploting(   ej, eSORj, eSORj1,    d,    k     ):
    x = [0 for _ in range(0, k+1)]
    for c in range(0,k+1):
        x[c] = c*d
    
    f = [0 for _ in range(0,k+1)]

    

    #plt.plot(x, f, color = 'maroon', label = 'factor with Jacobi')
    plt.plot(x, ej, color = 'green', label = 'Un-relaxed Jacobi')
    #plt.plot(x, eSORj, color = 'red', label = 'SRJ-2w')
    plt.plot(x, eSORj1, color = 'blue', label = 'SRJ-3w')


    plt.legend()

    plt.yscale('log')
    plt.xlabel('k')
    plt.ylabel('Error')
    plt.ylim(0,1)
    #plt.ylim(0,10**(-6))

    plt.title("k=%d"%k)





    plt.show()



def wlist0(w1,w2,w3,k):
    wlist = [w3 for _ in range(0,k+1)]

    for s in range(0,2):
        wlist[s] = w1
    for s in range(2,3):
        wlist[s] = w2

    return wlist





def altw(w1,w2,k):
    # wlist = [w1,w2,w3]

    wlist = [w1,w2]
    wselect =wlist[k%len(wlist)]
    
    return wselect





########################## main function below ########
########################## main function below ########
########################## main function below ########
########################## main function below ########




def main():
    import math


    
    ########################## Grid Information v ######
    ########################## Grid Information v ######
    
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

    k = 400

    #import Boundary Condition (BC)
    P = [[0 for _ in range(0,x_len+1)] for _ in range(0, y_len+1)]

    for i in range(0,x_len):
        sin = math.sin
        P[0][i] = sin(2*i*dx) + sin(5*i*dx) + sin(7*i*dx)
        P[y_len-1][i] = 0
    for j in range(0,y_len):
        P[j][0] = 0
        P[j][x_len-1] = 0
    
    ########################## Grid Information ^ ######



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


    #en = 10**(-3)
    en = 1


    #calculate iteration times of different w

    w1o = 0
    w2o = 0
    ko = 0
    SORejo = 0

        #P_input = copy.deepcopy(P)
    r = copy.deepcopy(P)
    e = [0 for _ in range(0,k+1)]
    ej = copy.deepcopy(e)

    P_Jin = copy.deepcopy(P)  
    P_Jin0 = copy.deepcopy(P)     
    P_Jin1 = copy.deepcopy(P)    


    eSORj = copy.deepcopy(e)
    eSORj0 = copy.deepcopy(e)
    eSORj1 = copy.deepcopy(e)



    for t in range (0,k+1):

        P_Jacobi = Jacobi(P_Jin, x_len, y_len)

        P_Jin = P_Jacobi

        r_J = Res(r,P_Jacobi,d,x_len,y_len)

        ej[t] = Error(r_J, x_len, y_len)
    
        ################# Plot SRJ w1=0.6 w2=1.1 #################
        
        w = altw(0.7,2.1,t) # Alt fn. odd iteration doing w1, even iteration doing w2
                 
        SORP_Jacobi0 = SORJacobi(P_Jin0, x_len, y_len, w)

        P_Jin0 = SORP_Jacobi0

        SORr_J0 = Res(r,SORP_Jacobi0,d,x_len,y_len)

        eSORj[t] = Error(SORr_J0, x_len, y_len)

        ################# Plot SRJ w1=0.6 w2=1.1 #################
        
        w1 = wlist0(9.3,1.1,1.0,k)[t] # Alt fn. odd iteration doing w1, even iteration doing w2
                 
        SORP_Jacobi1 = SORJacobi(P_Jin1, x_len, y_len, w1)

        P_Jin1 = SORP_Jacobi1

        SORr_J1 = Res(r,SORP_Jacobi1,d,x_len,y_len)

        eSORj1[t] = Error(SORr_J1, x_len, y_len)

    
    eploting(  ej, eSORj, eSORj1,    d,    k     )





# main function operator
if __name__ == "__main__":
    main()




























