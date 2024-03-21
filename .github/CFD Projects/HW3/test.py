import numpy as np

import math

# i = 5
# j = 5
# a = np.full(((i+1),(j+1)), 1, dtype = float)
# a[3,2] = 4
# print(a)

# F = np.full(3, 1, dtype = 'float')
# F[2] = 0
# print(F[-1])

# print(math.abs(1))

# u = np.array([[1,1,1],[1,0,1]])
# print(u[:,1])
a = np.array([0,1])
b = np.array([1,2])
c = np.abs(a-b)
print(c)