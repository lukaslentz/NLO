#
#
import numpy as np
import sympy as sp
#
xi, alpha, gamma, f, Omega, t = sp.symbols('xi alpha gamma f Omega t')
x = sp.Function('x')(t)
z = sp.Function('z')(t)
#
dgl = sp.diff(z+x,t,2) + xi*sp.diff(z+x,t) + alpha*(z+x) + gamma*(z+x)**3-f*sp.cos(Omega*t)
dgl = dgl.expand()
dgl = dgl.subs(f*sp.cos(Omega*t),sp.diff(x,t,2) + xi*sp.diff(x,t) + alpha*x + gamma*x**3)
#
print(dgl)
