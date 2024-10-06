#
#*******************************************NLO_WS23/24*******************************************
#                     Störungsrechnung nach Lindstedt/Poincare für den Duffingschwinger
#*************************************************************************************************
#Ertellt von: L.Lentz@umwelt-campus.de
#*************************************************************************************************
#
#Benötigte Bibliotheken laden
#******************************************************************************
import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from scipy.integrate import solve_ivp
#
def cut_order(exp,coeff,n):
   return exp.subs(mu,0) + sum([coeff**i*exp.coeff(coeff**i) for i in range(1,n+1)])
#
def time_set(omega):
    T = 2*np.pi/omega
    t_end = 5*T
    delta = T/100
    t_span = (0.0, t_end)
    t = np.arange(0.0, t_end, delta)
    return [t_span, t, T]
#
def duffing(t, state, om0, mu):
    x, v = state
    dx = v
    dv = -om0**2*x-mu*x**3
    return [dx, dv]
#
#Symbole definieren
t = Symbol('t')#Zeit
mu = symbols('mu' , negative=True)#Buchhaltungsparameter
omega = symbols('omega', positive=True)#Parameter der Dgl
om0, om1, om2 = symbols('om0 om1 om2', positive=True)
#
#Ansatz für omega0
eqom = Eq(omega,om0+mu*om1+mu**2*om2)
sol_om0 = solve(eqom,om0)[0]
#
#Funktionen definieren
x = Function('x')(t)#gesuchte Funktion
x0 = Function('x0')(t)#nullte Näherung
x1 = Function('x1')(t)#erste Näherung
x2 = Function('x2')(t)#zweite Näherung
#
#Ansatz für x
x = x0 + mu*x1 + mu**2*x2
#
#Dgl
dgl = expand(diff(x,t,2) + sol_om0**2*x + mu*x**3)
#
#Terme nach Größenordnungen von eps sortieren
dgl_0 = dgl.subs({mu:0})
dgl_1 = diff(dgl,mu).subs(mu,0)
dgl_2 = 1/2*diff(dgl,mu,2).subs(mu,0)
#
#Nullte Näherung lösen
#******************************************************************************
sol_x0 = dsolve(dgl_0,x0).rhs
#
A = Symbol('A', positive=True)
#
eq0_1 = Eq(sol_x0.subs(t,0),A)
eq0_2 = Eq(diff(sol_x0,t).subs(t,0),0)
const_0 = solve([eq0_1,eq0_2],symbols('C1 C2'))
#
sol_x0 = sol_x0.subs(const_0)
#
#Erste Näherung lösen
#******************************************************************************
dgl_1 = (dgl_1.subs(x0,sol_x0)).rewrite(sin, exp).rewrite(cos, exp).expand().rewrite(exp, sin).simplify()
#
sak_0 = dgl_1.coeff(cos(omega*t))
sol_sak_0 = solve(sak_0,om1)[0]
#
om0_eq_1 = Eq(cut_order((sol_om0**2).subs({om2:0,om1:sol_sak_0}).expand(),mu,1),om0**2)
sol_omega_1 = solve(om0_eq_1,omega)[0]
#
sol_x1 = dsolve(dgl_1.subs(om1,sol_sak_0),x1).rhs
#
eq1_1 = Eq(sol_x1.subs(t,0),0)
eq1_2 = Eq(diff(sol_x1,t).subs(t,0),0)
const_1 = solve([eq1_1,eq1_2],symbols('C1 C2'))
#
sol_x1 = sol_x1.subs(const_1)
#
#Zweite Näherung lösen
#******************************************************************************
dgl_2 = (dgl_2.subs({x0:sol_x0,x1:sol_x1,om1:sol_sak_0})).rewrite(sin, exp).rewrite(cos, exp).expand().rewrite(exp, sin).simplify()
#
sak_1 = dgl_2.coeff(cos(omega*t))
sol_sak_1 = solve(sak_1,om2)[0]
#
om0_eq_2 = Eq(cut_order((sol_om0**2).subs({om2:sol_sak_1,om1:sol_sak_0}).expand(),mu,2),om0**2)
sol_omega_2 = solve(om0_eq_2,omega)[1]
#
#Graphische Auswertung
#******************************************************************************
#Parameter
p_om = 1
p_mu = -1/6*p_om**2
p_A = 2
#
#Referenz Lösung aus numerischer Integration
[t_span, tvals, T] = time_set(0.5)
results_ivp_duffing = solve_ivp(duffing, t_span, [p_A,0], method='RK45', t_eval=tvals, args=(p_om,p_mu))
x_ref = results_ivp_duffing.y[0] 
#
#Referenz Lösung aus Datei
from csv import DictReader
with open('C:/SeafileContainer/Seafile/Meine Bibliothek/Lehre/NLO/Programme/Aufgabe 1/freq.csv', 'r') as f:
    ref_amp = [float(row["amp"]) for row in DictReader(f,delimiter='\t')]
with open('C:/SeafileContainer/Seafile/Meine Bibliothek/Lehre/NLO/Programme/Aufgabe 1/freq.csv', 'r') as f:   
    ref_freq = [float(row["freq"]) for row in DictReader(f,delimiter='\t')]
#
#Amplituden
sol_1 = x.subs({x0:sol_x0,x1:sol_x1,x2:0,omega:sol_omega_1}).subs({om0:p_om ,mu:p_mu,A:p_A})
sol_2 = x.subs({x0:sol_x0,x1:sol_x1,x2:0,omega:sol_omega_2}).subs({om0:p_om ,mu:p_mu,A:p_A})
#
x_sol_1 = lambdify(t, sol_1) 
x_sol_2 = lambdify(t, sol_2) 
#
#tvals = np.arange(0,100,0.01)
x_1 = x_sol_1(tvals)
x_2 = x_sol_2(tvals)
#
fig = plt.figure()
#
ax = fig.add_subplot(1, 1, 1)
ax.plot(tvals,x_1,color='blue')
ax.plot(tvals,x_2,color='red')
ax.plot(tvals,x_ref,color='black')
ax.set_xlabel('Zeit t')
ax.set_ylabel('Auslenkung x')
ax.legend(['1.Ordnung','2.Ordnung','Referenz'],loc='upper right')
ax.set_title('Störungsrechnung nach Lindstedt und Poincare')
# 
#
#Kreisfrequenz
om_sol_1 = lambdify(A,sol_omega_1.subs({mu:p_mu,om0:p_om}))
om_sol_2 = lambdify(A,sol_omega_2.subs({mu:p_mu,om0:p_om}))
# 
As = np.arange(0,2.1,0.1)
om_1 = om_sol_1(As)
om_2 = om_sol_2(As)
#
fig = plt.figure()
#
ax = fig.add_subplot(1, 1, 1)
ax.plot(om_1,As,color='blue')
ax.plot(om_2,As,color='red')
ax.plot(ref_freq,ref_amp,marker = 'o',markerfacecolor='black',markeredgecolor='black',linestyle = 'None')
ax.set_xlabel('Kreisfrequenz $\omega$')
ax.set_ylabel('Auslenkung x')
ax.legend(['1.Ordnung','2.Ordnung','Referenz'],loc='upper right')
ax.set_title('Störungsrechnung nach Lindstedt und Poincare')
#
plt.show()