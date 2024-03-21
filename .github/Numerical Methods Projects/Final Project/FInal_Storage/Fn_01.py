import copy 



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
            if ((P_xbl0[j][i] + P_xbh0[j][i]==a+b) or (P_ybl0[j][i] + P_ybh0[j][i]==c+d)): # find out grid point on the edge
                C_edge[j][i] = 1 

            if (P_xbl0[j][i] + P_ybl0[j][i] == a+c):     # find out where the edge cut x and y line, where x low y low
                C_xlyl[j][i] = 1 - C_edge[j][i]
            if (P_xbl0[j][i] + P_ybh0[j][i] == a+d):   # find out where the edge cut x and y line, where x low y high
                C_xlyh[j][i] = 1 - C_edge[j][i]
            if (P_xbh0[j][i] + P_ybl0[j][i] == b+c):   # find out where the edge cut x and y line, where x high y low
                C_xhyl[j][i] = 1 - C_edge[j][i]
            if (P_xbh0[j][i] + P_ybh0[j][i] == b+d):   # find out where the edge cut x and y line, where x high y high
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




def IBi( C_xbl, C_xbh, C_ybl, C_ybh, C_edge, C_xlyl, C_xlyh, C_xhyl, C_xhyh):

    C_total = copy.deepcopy(C_xbl)
    for j in range(0,j_Bmax+1):
        for i in range(0,i_Bmax+1):
            C_total[j][i] = C_xbl[j][i] + C_xbh[j][i] + C_ybl[j][i] + C_ybh[j][i] + C_edge[j][i] + C_xlyl[j][i] + C_xlyh[j][i] + C_xhyl[j][i] + C_xhyh[j][i]

    i_l = [0 for _ in range(0, j_max+1)]
    i_h = [0 for _ in range(0, j_max+1)]

    
    for j in range(j_Bmin,j_Bmax+1):
        for i in range(i_Bmin,i_Bmax+1):
            if C_total[j][i] == 1:
                i_l[j] = i
                break
        for i in range(j_Bmax,j_Bmin+1,-1):
            if C_total[j][i] == 1:
                i_h[j] = i
                break
    return i_l, i_h






        




def calculator_GSSOR(P_inputSOR, x_len, y_len,w): #Gauss Sadiel SOR method Pin-->Pout
    PgSOR = copy.deepcopy(P_inputSOR)
    for j in range(1,y_len-1):
        for i in range(1,x_len-1): 
            PgSOR[j][i] = \
                1/4 * (  PgSOR[j][i-1] + P_inputSOR[j][i+1] + PgSOR[j-1][i] + P_inputSOR[j+1][i] )*w + (1-w)*P_inputSOR[j][i]
    return PgSOR


    











########## input mathematical domain inforamtion
x_len = y_len = len = 1
r = 0.25
x_c = x_len/2
y_c = y_len/2


######### input numerical grid information
size = 32

d = len/size


######### create grid information
i_max = int(x_len/d)  # i in [0,i_max]
j_max = int(y_len/d)


P = [[0 for _ in range(0,j_max+1)] for _ in range(0,i_max+1)]




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


#print(C_edge)



C_total = copy.deepcopy(C_xbl)
for j in range(0,j_Bmax+1):
        for i in range(0,i_Bmax+1):
            C_total[j][i] = C_xbl[j][i] + C_xbh[j][i] + C_ybl[j][i] + C_ybh[j][i] + C_edge[j][i] + C_xlyl[j][i] + C_xlyh[j][i] + C_xhyl[j][i] + C_xhyh[j][i]

'''
C_total = copy.deepcopy(C_xbl)
for j in range(0,j_Bmax+1):
        for i in range(0,i_Bmax+1):
            C_total[j][i] = C_edge[j][i] + C_xlyl[j][i] + C_xlyh[j][i] + C_xhyl[j][i] + C_xhyh[j][i]
#print(C_total)
'''

i_l, i_h =IBi( C_xbl, C_xbh, C_ybl, C_ybh, C_edge, C_xlyl, C_xlyh, C_xhyl, C_xhyh)



#print(i_h)








print(C_total)



'''

# Test 
print(i_B1a)
print(i_B2b)

print(j_B1a)
print(j_B2b)



print(i_Bmin)



print(i_Bmax)
print(i_max)

'''
































