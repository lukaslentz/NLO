#
#*******************************************NLO_WS24/25*******************************************
#                           Einfache Störungsrechnung für den Duffingschwinger
#*************************************************************************************************
# Ertellt von: L.Lentz@umwelt-campus.de
#*************************************************************************************************
#
# Benötigte Bibliotheken laden
import numpy as np
import matplotlib.pyplot as plt
from sympy import *
#
#
# Symbole definieren
t = ...# Zeit
mu = ...# Buchhaltungsparameter
omega0 = ...# Parameter der Dgl
#
# Funktionen definieren
x = ...# gesuchte Funktion
x0 = ...# nullte Näherung
x1 = ...# erste Näherung
#
# Ansatz
x = ...
#
# Dgl
dgl = ...
#
# Terme nach Größenordnungen von eps sortieren
dgl_0 = ...
dgl_1 = ...
dgl_2 = ...
#
# Nullte Näherung lösen
sol_x0 = ... # Lösung der DGL
#
A = ...
#
# Bestimmung der Integrationskonstanten
eq0_1 = ...# Gleichung für x0(0)=A
eq0_2 = ...# Gleichung für v0(0)=0
const_0 = ...# Bestimmung der Konstanten
# 
sol_x0 = ...
# 
# 
#
# Erste Näherung lösen
sol_x1 = ...
#
# Bestimmung der Integrationskonstanten
eq1_1 = ...# Gleichung für x1(0)=0
eq1_2 = ...# Gleichung für v1(0)=0
const_1 = ...# Bestimmung der Konstanten
#
sol_x1 = ...
# 
#
# Graphische Auswertung
sol = ...
p_sol = ...
#
t = ...# Liste mit Zeitpunkten
y = ...# Funktionswerte

fig = plt.figure()
#
...
#
plt.show()
