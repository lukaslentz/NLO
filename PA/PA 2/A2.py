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
#Symbole definieren
t = Symbol('t')#Zeit
mu = symbols('mu' , negative=True)#Buchhaltungsparameter
omega = symbols('omega', positive=True)#Parameter der Dgl
om0, om1, om2 = symbols('om0 om1 om2', positive=True)
#
#Ansatz für omega0
eqom = Eq(omega, om0+mu*om1+mu**2*om2) # Gleichung für omega
sol_om0 = solve(eqom, om0)[0] # aufgelöst nach omega0
#
# #Funktionen definieren
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
dgl_1 = diff(dgl,mu).subs({mu:0})
dgl_2 = 1/2*diff(dgl,mu,2).subs({mu:0})
#
#Nullte Näherung lösen
#******************************************************************************
sol_x0 = dsolve(dgl_0,x0).rhs
# 
A = symbols('A', positive=True)
#
# Bestimmung der Integrationskonstanten
eq0_1 = Eq(sol_x0.subs(t,0),A)# Gleichung für x0(0)=A
eq0_2 = Eq(diff(sol_x0,t).subs(t,0),0)# Gleichung für v0(0)=0
const_0 = solve([eq0_1,eq0_2],symbols('C1 C2'))# Bestimmung der Konstanten
# 
sol_x0 = sol_x0.subs(const_0)
# 
#Erste Näherung lösen
#******************************************************************************
dgl_1 = (dgl_1.subs(x0,sol_x0)).rewrite(cos,exp).expand().rewrite(exp,cos).simplify() # Lösung für x0 in dgl_1 einsetzen und trigonometrische Potenzen auswerten
#
sak_0 = dgl_1.coeff(cos(omega*t)) # Säkularterm
sol_sak_0 = solve(sak_0,om1)[0] # om1 aus Säkularterm bestimmen 
#
sol_x1 = dsolve(dgl_1.subs(om1,sol_sak_0),x1).rhs # Lösung für x1 unter Berücksichtigung der Säkularterme
# 
# Bestimmung der Integrationskonstanten
eq1_1 = Eq(sol_x1.subs(t,0),0)# Gleichung für x1(0)=0
eq1_2 = Eq(diff(sol_x1,t).subs(t,0),0)# Gleichung für v1(0)=0
const_1 = solve([eq1_1,eq1_2],symbols('C1 C2'))# Bestimmung der Konstanten
#
sol_x1 = sol_x1.subs(const_1)
# 
# 
#
om0_eq_1 = Eq(cut_order((sol_om0**2).subs({om1:sol_sak_0,om2:0}).expand(),mu,1),om0**2) # Gleichung für om0 aufstellen, dabei mu nur bis hin zu linearen Termen berücksichtigen
sol_omega_1 = solve(om0_eq_1,omega)[0] # Gleichung nach omega lösen, liefert omega(A) in 1.Näherung
# 
#Zweite Näherung lösen
#******************************************************************************
dgl_2 = dgl_2.subs({x0:sol_x0,x1:sol_x1,om1:sol_sak_0}).rewrite(cos,exp).expand().rewrite(exp,cos).simplify() # Lösung für x0 und x1 in dgl_2 einsetzen und trigonometrische Potenzen auswerten
#
sak_1 = dgl_2.coeff(cos(omega*t)) # Säkularterm
sol_sak_1 = solve(sak_1,om2)[0] # om2 aus Säkularterm bestimmen 
#
om0_eq_2 = Eq(cut_order((sol_om0**2).subs({om1:sol_sak_0,om2:sol_sak_1}).expand(),mu,2),om0**2) # Gleichung für om1 aufstellen, dabei mu nur bis hin zu quadratischen Termen berücksichtigen
sol_omega_2 = solve(om0_eq_2,omega)[1] # Gleichung nach omega lösen, liefert omega(A) in 2.Näherung
#
#Graphische Auswertung
#******************************************************************************
#Parameter
p_om = 1
p_mu = -1/6*p_om**2
p_A = 2
#
#Amplituden
sol_1 = x.subs({x0:sol_x0,x1:sol_x1,x2:0,omega:sol_omega_1}).subs({om0:p_om,mu:p_mu,A:p_A})
sol_2 = x.subs({x0:sol_x0,x1:sol_x1,x2:0,omega:sol_omega_2}).subs({om0:p_om,mu:p_mu,A:p_A})
#  
#
x_sol_1 = lambdify(t, sol_1) 
x_sol_2 = lambdify(t, sol_2) 
#
tvals = np.arange(0,100,0.01)
x_1 = x_sol_1(tvals)
x_2 = x_sol_2(tvals)
#
fig = plt.figure()
#
ax = fig.add_subplot(1,1,1)
ax.plot(tvals,x_1,color='blue')
ax.plot(tvals,x_2,color='red')
ax.set_xlabel('Zeit t')
ax.set_ylabel('Auslenkung x')
ax.set_title('Störungsrechnung nach L/P')
ax.legend(['1.Ordnung','2.Ordnung'],loc='upper right')
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
ax = fig.add_subplot(1,1,1)
ax.plot(om_1,As,color='blue')
ax.plot(om_2,As,color='red')
ax.set_xlabel('Kreisfrequenzen')
ax.set_ylabel('Amplitude')
ax.set_title('Störungsrechnung nach L/P')
ax.legend(['1.Ordnung','2.Ordnung'],loc='upper right')
# 
plt.show()