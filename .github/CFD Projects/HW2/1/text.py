from sympy import *
init_printing()
x,a,b,k,d = symbols('x a b k d', real=True)

expr = (x+d)/(a+b*sin(k*x))**2
integral = integrate(expr,x)

pprint(integral)