import numpy as np
import matplotlib.pyplot as plt

def solve_burgers_FTCS(dx, dt, x_max, t_max, Pe):
    # 网格点数和时间步数
    Nx = int(x_max / dx) + 1
    Nt = int(t_max / dt)
    
    # 网格和时间步初始化
    x = np.linspace(0, x_max, Nx)
    u = np.zeros(Nx)
    u_new = np.zeros(Nx)
    
    # 初始条件
    u[:] = 0  # u(x, 0) = 0
    
    # FTCS迭代
    for n in range(1, Nt+1):
        for i in range(1, Nx-1):
            u_new[i] = u[i] + dt * (1/Pe) * (u[i+1] - 2*u[i] + u[i-1]) / (dx**2) - dt * u[i] * (u[i+1] - u[i-1]) / (2*dx)
        # 边界条件
        u_new[0] = 0  # u(0, t) = 0
        u_new[-1] = 1  # u(1, t) = 1
        u[:] = u_new
    
    return x, u

# 参数
Pe = 50
x_max = 1
t_max = 0.1 # 示例结束时间
dx = 1/20  # 最细的网格
dt = 0.0001  # 选择一个合适的时间步长

# 解决问题
x, u = solve_burgers_FTCS(dx, dt, x_max, t_max, Pe)

# 绘图
plt.plot(x, u, label=f'dx={dx}')
plt.xlabel('x')
plt.ylabel('u')
plt.title('Solution of Viscous Burger\'s Equation using FTCS')
plt.legend()
plt.show()
