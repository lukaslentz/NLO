#
# Duffing angetrieben
#
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from tqdm import tqdm # Fortschrittsbalken für Schleifen
import os
#
dir_script = os.path.dirname(os.path.abspath(__file__))
#
logfile = open(dir_script + '\logfile.txt','w')
logfile.write('*' * 50) 
logfile.write('\nProgramm gestartet\n')
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
# 
# Parameter
xi=0.1
alpha= 1
gamma=-0.1
f=0.2
#
eps = 1e-5
logfile.write('Schwellenwert fuer Einschwingen ist: ' + str(eps) + '\n')
#
# Werte der Erregerfrequenz, für die Lösungen berechnet werden
Omegas_up = np.arange(0.3,1.3,0.01)
Omegas_down = np.flip(Omegas_up)
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
     logfile.write(f'\t Omega = {Omega:.2f} \n')
     T = 2*np.pi/Omega # aktuelle Periodendauer
     finished = False
     cnt = 0
     while not finished and cnt<50:
        cnt += 1
        results_ivp_duffing = solve_ivp(duffing, [t0, t0 + 2*T], [x0,v0], max_step = 0.01, events = v_event, args = [xi, alpha, gamma, f, Omega])
        if results_ivp_duffing.success : # Testen ob Integration erfolgreich war
            t0 = results_ivp_duffing.t_events[0][-1]%T # Startzeit für die nächste Integration
            [x0,v0] = results_ivp_duffing.y_events[0][-1] # Anfangswerte für die nächste Integration
            if len(results_ivp_duffing.y_events[0])>=2: 
                delta = np.abs(x0 - results_ivp_duffing.y_events[0][-2][0])
                logfile.write(f'\t\tDurchlauf: {cnt}, aktuelle Delta : {delta} \n')
                if delta < eps:
                    finished = True
                    x_val_up.append(Omega) # Ergebnisse speichern 
                    y_val_up.append(x0)
                    logfile.write(f'\t\t\t {delta} < {eps} ---> Durchlauf beendet\n')
            else:
                logfile.write(f'\t\t-> Im Durchlauf: {cnt} wurde nur ein Wendepunkt gefunden \n')
        else: 
            break
    # 
#
logfile.write('\n\nSweep-up erfolgreich beendet\n\n\n')
#
#
logfile.write('Sweep-down wird gestartet\n')
# Sweep down
for Omega in tqdm(Omegas_down):
     logfile.write(f'\t Omega = {Omega:.2f} \n')
     T = 2*np.pi/Omega # aktuelle Periodendauer
     finished = False
     cnt = 0
     while not finished and cnt<50:
        cnt += 1
        results_ivp_duffing = solve_ivp(duffing, [t0, t0 + 2*T], [x0,v0], max_step = 0.01, events = v_event, args = [xi, alpha, gamma, f, Omega])
        if results_ivp_duffing.success : # Testen ob Integration erfolgreich war
            t0 = results_ivp_duffing.t_events[0][-1]%T # Startzeit für die nächste Integration
            [x0,v0] = results_ivp_duffing.y_events[0][-1] # Anfangswerte für die nächste Integration
            if len(results_ivp_duffing.y_events[0])>=2: 
                delta = np.abs(x0 - results_ivp_duffing.y_events[0][-2][0])
                logfile.write(f'\t\tDurchlauf: {cnt}, aktuelle Delta : {delta} \n')
                if delta < eps:
                    finished = True
                    x_val_down.append(Omega) # Ergebnisse speichern 
                    y_val_down.append(x0)
                    logfile.write(f'\t\t\t {delta} < {eps} ---> Durchlauf beendet\n')
            else:
                logfile.write(f'\t\t-> Im Durchlauf: {cnt} wurde nur ein Wendepunkt gefunden \n')
        else: 
            break
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