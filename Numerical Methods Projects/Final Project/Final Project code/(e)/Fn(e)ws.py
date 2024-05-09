import copy 
import matplotlib.pyplot as plt
#import numpy as np
import matplotlib as mpl


import tqdm
from tqdm import tqdm


##########################################################################################

def real_True(number): # turn complex number to negative
    if isinstance(number, int):
        return number
    elif isinstance(number, float):
        return number
    elif isinstance(number, complex):
        return -100

    
def int_True(number): # return int
    if isinstance(number, complex):
        return 0
    elif (number % 1==0):
        return -1
    else:
        return 0
    


def IBcboundary(i, j ,  x_len, y_len, d, r): # calculate boundary points for each i,j
    
    x = i * d
    y = j * d
    i_B1a_0= (x_len/2 - ( (r**2) - (y-y_len/2)**2 )**(1/2))/d
    i_B2b_0=(x_len/2 + ( (r**2) - (y-y_len/2)**2 )**(1/2))/d
    j_B1a_0=(y_len/2 - ( (r**2) - (x-x_len/2)**2 )**(1/2))/d
    j_B2b_0=(y_len/2 + ( (r**2) - (x-x_len/2)**2 )**(1/2))/d


    i_B1a= int(real_True(i_B1a_0)) 
    i_B2b=int(real_True(i_B2b_0)) +1 + int_True(i_B2b_0)
    j_B1a=int(real_True(j_B1a_0))
    j_B2b=int(real_True(j_B2b_0))+1 + int_True(j_B2b_0)



    return i_B1a, i_B2b, j_B1a, j_B2b
##########################################################################################




##########################################################################################

def IBcondiitionor(  i_B1a, i_B2b, j_B1a, j_B2b, j_max, i_max): # use each EP i,j to classify each EP
    GeneralP = [[0 for _ in range(0,j_max+1)] for _ in range(0,i_max+1)]
    P_xbl0 = copy.deepcopy(GeneralP)
    P_xbh0 = copy.deepcopy(GeneralP)
    P_ybl0 = copy.deepcopy(GeneralP)
    P_ybh0 = copy.deepcopy(GeneralP)


    C_xbl = copy.deepcopy(GeneralP)
    C_xbh = copy.deepcopy(GeneralP)
    C_ybl = copy.deepcopy(GeneralP)
    C_ybh = copy.deepcopy(GeneralP)

    C_edge = copy.deepcopy(GeneralP)

    C_xlyl = copy.deepcopy(GeneralP)
    C_xlyh = copy.deepcopy(GeneralP)
    C_xhyl = copy.deepcopy(GeneralP)
    C_xhyh = copy.deepcopy(GeneralP)

    

    a = 1
    b = 10
    c = 41
    d = 51


    for j in range(1,j_max):
        if i_B1a[j]>0:
            P_xbl0[j][i_B1a[j]] = a
        if i_B2b[j]>0:
            P_xbh0[j][i_B2b[j]] = b
    for i in range(1,i_max):
        if j_B1a[i]>0:
            P_ybl0[j_B1a[i]][i] = c
        if j_B2b[i]>0:
            P_ybh0[j_B2b[i]][i] = d

    for j in range(1,j_max):
        for i in range(1,i_max):
            if ((P_xbl0[j][i] + P_xbh0[j][i]==a+b) or (P_ybl0[j][i] + P_ybh0[j][i]==c+d)): # find out grid point just on the edge -->C_edge
                C_edge[j][i] = 1  

            if (P_xbl0[j][i] + P_ybl0[j][i] == a+c):     # find out where the edge cut x and y line, where x low y low -->C_xlyl
                C_xlyl[j][i] = 1 - C_edge[j][i]
            if (P_xbl0[j][i] + P_ybh0[j][i] == a+d):   # find out where the edge cut x and y line, where x low y high -->C_xlyh
                C_xlyh[j][i] = 1 - C_edge[j][i]
            if (P_xbh0[j][i] + P_ybl0[j][i] == b+c):   # find out where the edge cut x and y line, where x high y low -->C_xhyl
                C_xhyl[j][i] = 1 - C_edge[j][i]
            if (P_xbh0[j][i] + P_ybh0[j][i] == b+d):   # find out where the edge cut x and y line, where x high y high -->C_xhyh
                C_xhyh[j][i] = 1 - C_edge[j][i]
            
            special = C_edge[j][i] +  C_xlyl[j][i] + C_xlyh[j][i] + C_xhyl[j][i]  + C_xhyh[j][i]
        
    for j in range(1,j_max):
        for i in range(1,int(i_max/2)):

            C_xbl[j][i] = P_xbl0[j][i]//a - C_edge[j][i] -  C_xlyl[j][i] - C_xlyh[j][i] # find out where edge ONLY cut x low lines 
    
    for j in range(1,j_max):
        for i in range(int(i_max/2),i_max):

            C_xbh[j][i] = P_xbh0[j][i]//b - C_edge[j][i]  - C_xhyl[j][i] - C_xhyh[j][i] # find out where edge ONLY cut x high lines 

    for j in range(1,int(j_max/2)):
        for i in range(1,i_max):
            C_ybl[j][i] = P_ybl0[j][i]//c - C_edge[j][i] - C_xlyl[j][i]  - C_xhyl[j][i] # find out where edge ONLY cut y low lines

    for j in range(int(j_max/2),j_max):
        for i in range(1,i_max):
            C_ybh[j][i] = P_ybh0[j][i]//d - C_edge[j][i]  - C_xlyh[j][i]  - C_xhyh[j][i] # find out where edge ONLY cut y high lines 

    

    return C_xbl, C_xbh, C_ybl, C_ybh, C_edge, C_xlyl, C_xlyh, C_xhyl, C_xhyh
##########################################################################################








##########################################################################################
def Cinner_and_IBi(i_max, j_max, i_Bmin, i_Bmax, j_Bmin, j_Bmax, C_xbl, C_xbh, C_ybl, C_ybh, C_edge, C_xlyl, C_xlyh, C_xhyl, C_xhyh):

    C_total = copy.deepcopy(C_xbl)
    C_inner = copy.deepcopy(C_xbl)

    for j in range(0,j_Bmax+1):
        for i in range(0,i_Bmax+1):
            C_total[j][i] = C_xbl[j][i] + C_xbh[j][i] + C_ybl[j][i] + C_ybh[j][i] + C_edge[j][i] + C_xlyl[j][i] + C_xlyh[j][i] + C_xhyl[j][i] + C_xhyh[j][i]
            C_inner[j][i] = 1 - C_total[j][i]

    i_l = [0 for _ in range(0, j_max+1)]
    i_h = [0 for _ in range(0, j_max+1)]



    for j in range(0,j_Bmin):
        for i in range(0,i_max):
            C_inner[j][i] = 0

    for j in range(j_Bmax+1,j_max):
        for i in range(0,i_max):
            C_inner[j][i] = 0

    
    for j in range(j_Bmin,j_Bmax+1):
        for i in range(0,i_max+1):
            if C_total[j][i] == 1:
                i_l[j] = i
                break
            else: C_inner[j][i] -= 1
        for i in range(i_max,i_Bmin,-1):
            if C_total[j][i] == 1:
                i_h[j] = i
                break
            else: C_inner[j][i] -= 1

    for j in range(0,j_max+1):
        for i in range(0,i_max+1):
            if C_inner[j][i] < 0: C_inner[j][i] = 0
    
    return C_inner, C_total, i_l, i_h # C_total means the total inner boundary, while C_inner is the points in the inner boundary
##########################################################################################


    
##########################################################################################
def calculator_GSSOR(P_inputSOR, j_limlow, j_limhigh, i_limlow, i_limhigh , w): #Gauss Sadiel SOR method Pin-->Pout
    PgSOR = copy.deepcopy(P_inputSOR)
    for j in range(j_limlow, j_limhigh+1):
        for i in range(i_limlow, i_limhigh+1): 
            PgSOR[j][i] = \
                1/4 * (  PgSOR[j][i-1] + P_inputSOR[j][i+1] + PgSOR[j-1][i] + P_inputSOR[j+1][i] )*w + (1-w)*P_inputSOR[j][i]
    return PgSOR




def Point_Calculator_IBGSSOR(P, PgSOR, j, i, w, C_inner, C_total, C_xbl, C_xbh, C_ybl, C_ybh, C_edge, C_xlyl, C_xlyh, C_xhyl, C_xhyh):
    

    PgSOR[j][i] = \
        (1/4 * (  PgSOR[j][i-1] + P[j][i+1] + PgSOR[j-1][i] + P[j+1][i] )*w + (1-w)*P[j][i])*(1-C_total[j][i])*(1-C_inner[j][i]) +(\
            C_edge[j][i] \
        + C_xbl[j][i] * IB1_cal(PgSOR[j][i-1],PgSOR[j-1][i],P[j+1][i])\
        + C_xbh[j][i] * IB1_cal(P[j][i+1],PgSOR[j-1][i],P[j+1][i])\
        + C_ybl[j][i] * IB1_cal(PgSOR[j-1][i],PgSOR[j][i-1],P[j][i+1])\
        + C_ybh[j][i] * IB1_cal(P[j+1][i],PgSOR[j][i-1],P[j][i+1])\
        + C_xlyl[j][i] * IB2_cal(PgSOR[j][i-1],PgSOR[j-1][i])\
        + C_xlyh[j][i] * IB2_cal(PgSOR[j][i-1],P[j+1][i])\
        + C_xhyl[j][i] * IB2_cal(P[j][i+1],PgSOR[j-1][i])\
        + C_xhyh[j][i] * IB2_cal(P[j][i+1],P[j+1][i]) ) *(1-C_inner[j][i])\
        + C_inner[j][i] * 1
    
    return PgSOR[j][i]
##########################################################################################

##########################################################################################
def IB1_cal(a,b,c):
    P = (a*4/3+8/3+b+c)/6
    return P
##########################################################################################

##########################################################################################
def IB2_cal(a,b):
    P = (a+b)/6 + 2/3
    return P
##########################################################################################
    
        
    








##########################################################################################
def GSSOR_Solver(P_input, i_max, j_max, i_Bmin, i_Bmax, j_Bmin, j_Bmax, i_l, i_h, w , \
                 C_inner, C_total, C_xbl, C_xbh, C_ybl, C_ybh, C_edge, C_xlyl, C_xlyh, C_xhyl, C_xhyh):

    # Solve for low GP
    P = copy.deepcopy(P_input)
    
    P = calculator_GSSOR(P, 1 , j_Bmin-1, 1 , i_max-1 , w)

    # Solve for EP + GP in middle
    PgSOR = copy.deepcopy(P)
    for j in range(j_Bmin, j_Bmax+1):


        # solve for left GP
        for i in range(1, i_l[j]): 
            PgSOR[j][i] = \
                1/4 * (  PgSOR[j][i-1] + P[j][i+1] + PgSOR[j-1][i] + P[j+1][i] )*w + (1-w)*P[j][i]
        P = PgSOR


        # solve for middle (including point judgement (mathematically))
        for i in range(i_l[j], i_h[j]+1):
            PgSOR[j][i] = Point_Calculator_IBGSSOR(P, PgSOR, j, i, w,\
                                               C_inner, C_total, C_xbl, C_xbh, C_ybl, C_ybh, C_edge, C_xlyl, C_xlyh, C_xhyl, C_xhyh)
        P = PgSOR
            
        # solve for right GP
        for i in range(i_h[j]+1, i_max):
            PgSOR[j][i] = \
                1/4 * (  PgSOR[j][i-1] + P[j][i+1] + PgSOR[j-1][i] + P[j+1][i] )*w + (1-w)*P[j][i]
        P = PgSOR
    
    # solve for high GP
    P = calculator_GSSOR(P, j_Bmax+1 , j_max-1 ,1, i_max-1 , w)

    return P
    
##########################################################################################



def Res(r,Pin,d,i_l,i_h,i_max,j_max): #calcuate residual (Pin, source)-->Rout
    P = copy.deepcopy(Pin)
    rc = copy.deepcopy(r)
    for j in range(1,j_max):
        for i in range(1,i_l[j]): 
            rc[j][i] = \
                ( P[j][i-1] + P[j][i+1] + P[j-1][i] + P[j+1][i] - 4*P[j][i] ) #*(d**2)
        for i in range(i_h[j]+1,i_max): 
            rc[j][i] = \
                ( P[j][i-1] + P[j][i+1] + P[j-1][i] + P[j+1][i] - 4*P[j][i] ) #*(d**2)
    return rc

def Error(r_in, i_max, j_max): #calculate error Rin-->e out
    r = copy.deepcopy(r_in)
    e = 0
    for j in range(1,j_max):
        for i in range(1,i_max): 
            e += abs(r[j][i])   #**2
    return e



############


def eploting(   Kw, w_min, w_max, w_step  ):
    x = copy.deepcopy(Kw)

    for c in range(0, len(Kw)):
        x[c] = w_min + c*w_step
    plt.plot(x, Kw, color = 'orangered', )
    plt.legend()

    #plt.yscale('log')
    plt.xlabel('w')
    plt.ylabel('Iteration Steps')
    #plt.ylim(0,1)
    #plt.ylim(0,10**(-6))

    plt.title('convergence number with overrelaxization parameter')





    plt.show()


 
        

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


def main():

    ########## input mathematical domain inforamtion
    x_len = y_len = len = 1
    r = 0.25
    x_c = x_len/2
    y_c = y_len/2


    ######### input numerical grid information
    size = 96

    d = len/size


    ######### create grid information
    i_max = int(x_len/d)  # i in [0,i_max]
    j_max = int(y_len/d)


    P_0 = [[0 for _ in range(0,j_max+1)] for _ in range(0,i_max+1)]

    import random
    for j in range(1,y_len-1):
        for i in range(1,x_len-1):
            P_0[j][i] = random.randint(-1,1)




    ######### define inner boundary limit
    i_Bmin = int((x_c-r)/d)
    i_Bmax = int((x_c+r)/d) +1 +  int_True((x_c+r)/d)
    j_Bmin = int((y_c-r)/d)
    j_Bmax = int((y_c+r)/d) +1 +  int_True((y_c+r)/d)


    ######### define inner boundary
    iline = [0 for _ in range(0,i_max+1)]
    jline = [0 for _ in range(0,j_max+1)]

    i_B1a = copy.deepcopy(iline)
    i_B2b = copy.deepcopy(iline)
    j_B1a = copy.deepcopy(jline)
    j_B2b = copy.deepcopy(jline)



    for j in range(j_Bmin,j_Bmax+1):
        for i in range(i_Bmin, i_Bmax+1):
            i_B1a[j], i_B2b[j], j_B1a[i], j_B2b[i] = IBcboundary(i, j ,  x_len, y_len, d, r)


    C_xbl, C_xbh, C_ybl, C_ybh, C_edge, C_xlyl, C_xlyh, C_xhyl, C_xhyh = IBcondiitionor(  i_B1a, i_B2b, j_B1a, j_B2b, j_max, i_max)

    C_inner, C_total, i_l, i_h = Cinner_and_IBi(i_max, j_max, i_Bmin, i_Bmax, j_Bmin, j_Bmax, \
                                                C_xbl, C_xbh, C_ybl, C_ybh, C_edge, C_xlyl, C_xlyh, C_xhyl, C_xhyh)     




    ##################################################### Iteration Part ###################################################################################

    #w =1
    Kmax = 4000

    refE = 10**(-3)


    r1_0 = copy.deepcopy(P_0)
    e1_0 = [0 for _ in range(0,Kmax+1)]



    


    w_min = 1
    w_max = 2
    w_step = 0.1

    Kw = [0 for _ in range(  0, int((w_max-w_min)/w_step+1)  )]

    
    for w_t in tqdm(range( 0, int((w_max-w_min)/w_step+1)  ), desc = 'General'):

        w = w_min + w_t*w_step

        P_input = copy.deepcopy(P_0)
        r1 = copy.deepcopy(r1_0)
        e1 = copy.deepcopy(e1_0)

        for k in tqdm(range(0,Kmax+1)):

            P = GSSOR_Solver(P_input, i_max, j_max, i_Bmin, i_Bmax, j_Bmin, j_Bmax, i_l, i_h, w , \
                    C_inner, C_total, C_xbl, C_xbh, C_ybl, C_ybh, C_edge, C_xlyl, C_xlyh, C_xhyl, C_xhyh)
            P_input = copy.deepcopy(P)

            r1 = Res(r1,P,d,i_l,i_h,i_max,j_max)
            e1[k] = Error(r1, i_max, j_max)

            if e1[k] < refE:
                Kw[w_t] = k
                break
            if k == Kmax:
                Kw[w_t] = None


    eploting(   Kw, w_min, w_max, w_step  )
    print(Kw)
    


    




















'''



'''



if __name__ == "__main__":
    main()





















