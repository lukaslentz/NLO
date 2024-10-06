#
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from tqdm import tqdm
#
t = sp.Symbol('t')
x = sp.Function('x')(t)
z = sp.Function('z')(t)
#
xi, alpha, gamma, f = 0.1, 1, -0.1, 0.2
#
#   HB-Solutions
#*********************************************************************************************************
#
def hb_solutions(Omega):
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
    return  [sp.lambdify(t,x_hb.subs({a:sol[0],b:sol[1]})) for sol in hb_sols],[sp.sqrt((a**2+b**2).subs({a:sol[0],b:sol[1]})) for sol in hb_sols]

#
#   Stability
#*********************************************************************************************************
#
def get_A(t, sol):
    return([[0,1],[-alpha-3*gamma*sol(t)**2,-xi]])
#
def F(t, state, sol):
    return np.dot(get_A(t, sol),state)
#
def get_M(sol,T):   
    sol_1 =  solve_ivp(F, [0, T], [1, 0], method = 'RK45', max_step = 0.01, args = [sol])
    sol_2 =  solve_ivp(F, [0, T], [0, 1], method = 'RK45', max_step = 0.01, args = [sol])
    return [[sol_1.y[0][-1],sol_1.y[1][-1]],[sol_2.y[0][-1],sol_2.y[1][-1]]]
#
def is_stable(sol,Omega):
    eigs,_ = np.linalg.eig(get_M(sol,np.pi/Omega))
    return 1 > np.max(np.abs(eigs))
#
Omega_vals = np.arange(0.3,1.3,0.01)
stable_Omega = []
instable_Omega = []
stable_amp = []
instable_amp = []
#
for Omega in  tqdm(Omega_vals):
    hb_sol, hb_amp = hb_solutions(Omega)
    for i,sol in enumerate(hb_sol):
        if is_stable(sol,Omega):
            stable_Omega.append(Omega)
            stable_amp.append(hb_amp[i])
        else:
            instable_Omega.append(Omega)
            instable_amp.append(hb_amp[i])
#
fig = plt.figure()
#
ax = fig.add_subplot(1, 1, 1)
ax.plot(stable_Omega,stable_amp,marker= 'o',markerfacecolor='blue',markeredgecolor='blue',linestyle = 'None')
ax.plot(instable_Omega,instable_amp,marker= 'o',markerfacecolor='gray',markeredgecolor='gray',linestyle = 'None')
ax.set_xlabel('Erregerkreisfrequenz $\Omega$')
ax.set_ylabel('Amplitude $C$')
ax.legend(['stabil','instabil'],loc='upper right')
#
plt.show()
