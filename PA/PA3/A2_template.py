#
# Duffing angetrieben
#
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from tqdm import tqdm # Fortschrittsbalken für Schleifen
#
logfile = ...
logfile.write('Programm gestartet\n')
#
logfile.write('Funktionen werden definiert\n')
# Duffing-Schwinger als System 1. Ordnung
def duffing(t, state, xi, alpha, gamma, f, Omega):
    x, v = state
    dx = v
    dv = -xi*v - alpha*x - gamma*x**3 + f*np.cos(Omega*t) 
    return [dx, dv]
#
# Event zum Auffinden der Umkehrpunkte mit positivem VZ
def v_event(t, state, xi, alpha, gamma, f, Omega):
    return state[1]
v_event.direction = -1
#
logfile.write('Parameter werden definiert\n')
# Parameter
xi=0.1
alpha= 1
gamma=-0.1
f=0.2
#
eps = ...
logfile.write('Schwellenwert für Einschwingen ist: ' + str(eps) + '\n')
#
# Werte der Erregerfrequenz, für die Lösungen berechnet werden
Omegas_up = ...
Omegas_down = ...
#
# Felder für Ergebnisse anlegen
x_val_up = [] # Omega Werte
y_val_up = [] # maximale Auslenkung
x_val_down = []
y_val_down = []
#
# Initiale Anfangsbedingungen
[t0,x0,v0] = [0,0.1,0]
#
logfile.write('Sweep-up wird gestartet\n')
# Sweep up
for Omega in tqdm(Omegas_up):
    ...
#
logfile.write('Sweep-up erfolgreich beendet\n\n\n')
#
#
logfile.write('Sweep-down wird gestartet\n')
# Sweep up
for Omega in tqdm(Omegas_down):
    ...
#
logfile.write('Sweep-down erfolgreich beendet\n\n\n')
#
#
fig = plt.figure()
#
ax = fig.add_subplot(1, 1, 1)
ax.plot(x_val_up,y_val_up,marker = 'o',markerfacecolor='red',markeredgecolor='red',linestyle = 'None')
ax.plot(x_val_down,y_val_down,marker = 'o',markerfacecolor='blue',markeredgecolor='blue',linestyle = 'None')
ax.set_ylabel('Maximale Auslenkung')
ax.set_xlabel('Erregerkreisfrequenz $\Omega$')
ax.set_title('Sweep Duffing-Schwinger')
#
plt.show()
#
logfile.close()
#