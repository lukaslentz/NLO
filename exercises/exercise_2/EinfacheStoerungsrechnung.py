#
#*******************************************NLO_WS23/24*******************************************
#                           Einfache Störungsrechnung für den Duffingschwinger
#*************************************************************************************************
#Ertellt von: L.Lentz@umwelt-campus.de
#*************************************************************************************************

#Benötigte Bibliotheken laden
#******************************************************************************
import numpy as np
import matplotlib.pyplot as plt
from sympy import *


#Symbole definieren
t = Symbol('t')#Zeit
mu = symbols('mu')#Buchhaltungsparameter
omega = symbols('omega', positive=True)#Parameter der Dgl

#Funktionen definieren
x = Function('x')(t)#gesuchte Funktion
x0 = Function('x0')(t)#nullte Näherung
x1 = Function('x1')(t)#erste Näherung

#Ansatz
x = x0 + mu*x1 

#Dgl
dgl = expand(diff(x,t,2)+omega**2*x+mu*x**3)

#Terme nach Größenordnungen von eps sortieren
dgl0 = dgl.subs({mu:0})
dgl1 = diff(dgl,mu).subs(mu,0)
dgl2 = 1/2*diff(dgl,mu,2).subs(mu,0)

#Nullte Näherung lösen
solx0 = dsolve(dgl0,x0).rhs

A = Symbol('A', positive=True)

eq0_1 = Eq(solx0.subs(t,0),A)
eq0_2 = Eq(diff(solx0,t).subs(t,0),0)
const0 = solve([eq0_1,eq0_2],symbols('C1 C2'))

solx0 = solx0.subs(const0)


#Erste Näherung lösen
solx1 = dsolve(dgl1.subs(x0,solx0),x1).rhs

eq1_1 = Eq(solx1.subs(t,0),0)
eq1_2 = Eq(diff(solx1,t).subs(t,0),0)
const1 = solve([eq1_1,eq1_2],symbols('C1 C2'))

solx1 = solx1.subs(const1)

sol = x.subs({x0:solx0,x1:solx1}).subs({omega:1,mu:-1/6,A:1})
psol = lambdify(t, sol) 

t = np.arange(0,100,0.01)
y = psol(t)

fig = plt.figure()
#
ax = fig.add_subplot(1, 1, 1)
ax.plot(t,y,color='blue')
ax.set_xlabel('Zeit t')
ax.set_ylabel('Auslenkung x')
ax.set_title('Einfache Störungsrechnung')
#
plt.show()


print(psol(1))