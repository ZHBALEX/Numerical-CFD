import math
import numpy as np
import matplotlib.pyplot as plt

def main():
    x_max = 2 * math.pi
    t_max = 60
    dx_values = [0.005, 0.01, 0.02]
    m = 2

    def sin(x):
        return np.sin(x)

    def cos(x):
        return np.cos(x)

    x = np.arange(0, math.pi + 0.005, 0.005)
    y_exact = x  # Exact k

    plt.plot(y_exact, x, color='blue', label='Exact k', linestyle='-')

    y1 = np.sin(x)
    y2 = np.sqrt(2*(1-np.cos(x)))
    
    plt.plot(x, y1, color='pink', label='1st Central Difference k*', linestyle='-')
    plt.plot( x, y2, color='red', label='2nd Central Difference k*', linestyle='-')
    
    
    plt.xlabel('kΔx')
    plt.ylabel('k*Δx')
    plt.title('Modified Wavenumber')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()