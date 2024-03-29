import sympy as sp

KR, rho, VMin, n = sp.symbols('KR rho VMin n')
r_0 = 3
R = 1
k = sp.pi
x = 2
Delta_x = 0.1  
Q0 = 0

n_max = int(x / Delta_x)

sum_expr = sp.Sum((rho * Q0 + 2 * sp.pi * (r_0 * n * Delta_x - R / k * sp.cos(k * n * Delta_x)) * VMin) / 
                  ((r_0 + R * sp.sin(k * n * Delta_x))**4), (n, 1, n_max))

delta_P = -KR / sp.pi**2 * sum_expr

delta_P_evaluated = delta_P.doit()

print("ΔP =")
sp.pprint(delta_P_evaluated, use_unicode=True)

