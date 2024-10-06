#
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
#
t = sp.Symbol('t')
x = sp.Function('x')(t)
z = sp.Function('z')(t)
#
xi, alpha, gamma, f, Omega = 0.1, 1, -0.1, 0.2, 0.7
#
#   HB-Solutions
#*********************************************************************************************************
#
a,b = sp. symbols('a b')
x_hb = sp.Function('x')(t)
x_hb = a*sp.cos(Omega*t)+b*sp.sin(Omega*t)
#
res = (sp.diff(x_hb,t,2) + xi*sp.diff(x_hb,t) + alpha*x_hb + gamma*x_hb**3-f*sp.cos(Omega*t)).rewrite(sp.cos,sp.exp).rewrite(sp.sin,sp.exp).expand().rewrite(sp.exp,sp.cos).rewrite(sp.exp,sp.sin).simplify()
eq1 = res.coeff(sp.cos(Omega*t))
eq2 = res.coeff(sp.sin(Omega*t))
hb_all_sols = sp.solve([eq1,eq2],[a,b])
hb_sols = []
eps = 1e-19
for sol in hb_all_sols:
    if abs(sp.im(sol[0])) < eps and abs(sp.im(sol[1])) < eps:
        hb_sols.append([sp.re(sol[0]),sp.re(sol[1])])
hb_xs = [ sp.lambdify(t,x_hb.subs({a:round(sol[0],10),b:round(sol[1],10)})) for sol in hb_sols]
#
#   Stability
#*********************************************************************************************************
#
def get_A(t, sol):
    return([[0,1],[-alpha-3*gamma*sol(t)**2,-xi]])
#
def f(t, state, sol):
    return np.dot(get_A(t, sol),state)
#
def get_M(sol):   
    T = np.pi/Omega
    sol_1 =  solve_ivp(f, [0, T], [1, 0], method = 'RK45', max_step = 0.01, args = [sol])
    sol_2 =  solve_ivp(f, [0, T], [0, 1], method = 'RK45', max_step = 0.01, args = [sol])
    return [[sol_1.y[0][-1],sol_1.y[1][-1]],[sol_2.y[0][-1],sol_2.y[1][-1]]]
#
def is_stable(sol):
    eigs,_ = np.linalg.eig(get_M(sol))
    return 1 > np.max(np.abs(eigs))

stability_bool = [is_stable(sol) for sol in hb_xs]

print(stability_bool)




# print((x_hb**2).rewrite(sp.cos,sp.exp).rewrite(sp.sin,sp.exp).expand().rewrite(sp.exp,sp.cos).rewrite(sp.exp,sp.sin).simplify())
# #
# T = np.pi/Omega
# t_vals = np.arange(0,T,0.1)
# #
# As = get_A(t_vals,hb_xs[0])
# #
# fig = plt.figure()
# #
# ax = fig.add_subplot(1, 1, 1)
# ax.plot(t_vals,As[1][0],color='blue')
# #
# plt.show()
