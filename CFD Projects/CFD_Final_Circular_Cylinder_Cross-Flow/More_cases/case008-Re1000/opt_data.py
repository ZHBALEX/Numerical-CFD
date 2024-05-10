import numpy as np
import matplotlib.pyplot as plt

from N32_Re1000_8x4 import *

def Math_FFT(value):
    # compute DFT with optimized FFT
    w = np.abs(np.fft.fft(value))
    n = np.size(w)

    half = int(n/2)
    x = np.linspace(0,1, np.size(w))
    x = x/dt

    w = w[0:half]

    w_max = np.argmax(w[1:])/(n-1)/dt


    print("w_max = ", w_max)

    plt.plot(x[:half],w[:half])
    plt.show()

    


    return w_max






dt = 0.01

lift_coeff = Get_Lift_coeff()
drag_coeff = Get_Drag_coeff()

lift_coeff[:] = lift_coeff[:]
drag_coeff[:] = drag_coeff[:]

print('mean lift_coeff=', np.max(lift_coeff))
print('mean drag_coeff=', np.mean(drag_coeff))


w1 = Math_FFT(lift_coeff)
w2 = Math_FFT(drag_coeff)



# print('lift coefficient frequency f=', f_cl)
# print('drag coefficient frequency f=', f_cd)


dt = 0.01
St = 0.4 * 0.5/1 
print('St=', St)













