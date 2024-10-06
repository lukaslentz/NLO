#
# Duffing angetrieben
#
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
#
logfile = open("C:/SeafileContainer/Seafile/Meine Bibliothek/Lehre/NLO/NLO_23_24/Programme/Aufgabe 3/log_file.txt","a")
logfile.write('Programm started\n')
#
logfile.write('Definition of functions\n')
def duffing(t, state, xi, alpha, gamma, f, Omega):
    x, v = state
    dx = v
    dv = -xi*v - alpha*x - gamma*x**3 + f*np.cos(Omega*t) 
    return [dx, dv]
#
def v_event(t, state, xi, alpha, gamma, f, Omega):
    return state[1]
v_event.direction = -1
#
logfile.write('Set Parameters\n')
xi=0.1
alpha= 1
gamma=-0.1
f=0.2
#
eps = 1e-3
#
Omegas_up = np.arange(0.3,1.3,0.01)
Omegas_down = np.flip(Omegas_up)
#
[t0,x0,v0] = [0,0.1,0.1]
x_val_up = []
y_val_up = []
x_val_down = []
y_val_down = []
# 
from tqdm import tqdm
#
logfile.write('Start sweep-up\n')
for Omega in tqdm(Omegas_up):
    logfile.write(f'\tOmega={Omega}\tt0={t0}\tx0={x0}\tv0={v0}\n')
    T = 2*np.pi/Omega
    finished = False
    cnt = 1
    while not finished: 
        results_ivp_duffing = solve_ivp(duffing, [t0,t0+10*T], [x0,v0], max_step = 0.01, events = v_event, args = [xi, alpha, gamma, f, Omega])
        if results_ivp_duffing.success == True:
            t0 = (results_ivp_duffing.t_events[0][0])%T
            [x0,v0] = results_ivp_duffing.y_events[0][-1]
            delta = min([abs(x0-event[0]) for event in results_ivp_duffing.y_events[0][0:-2]])
            print(delta)
            if delta < eps:
                finished = True
                x_val_up.append(Omega)
                y_val_up.append(x0) 
                logfile.write(f'\tfinished computation after {cnt} integration loops with a delta of {delta}\n')
            cnt += 1
        else:
            break        
#
logfile.write('Start sweep-down\n')
for Omega in tqdm(Omegas_down):
    logfile.write(f'\tOmega={Omega}\n')
    T = 2*np.pi/Omega
    finished = False
    cnt = 1
    while not finished: 
        results_ivp_duffing = solve_ivp(duffing, [t0,t0+10*T], [x0,v0], max_step = 0.01, events = v_event, args = [xi, alpha, gamma, f, Omega])
        if results_ivp_duffing.success == True:
            t0 = (results_ivp_duffing.t_events[0][0])%T
            [x0,v0] = results_ivp_duffing.y_events[0][-1]
            if abs(results_ivp_duffing.y_events[0][-1,0]-results_ivp_duffing.y_events[0][-2,0]) < eps:
                finished = True
                x_val_down.append(Omega)
                y_val_down.append(x0)
                logfile.write(f'\tfinished computation after {cnt} integration loops\n')
            cnt += 1
        else:
            break 
#
logfile.write('Computation finished\n')
logfile.write('Start graphical evaluation\n')
#
fig = plt.figure()
#
ax = fig.add_subplot(1, 1, 1)
ax.plot(x_val_up,y_val_up,marker = 'o',markerfacecolor='red',markeredgecolor='red',linestyle = 'None')
ax.plot(x_val_down,y_val_down,marker = 'o',markerfacecolor='blue',markeredgecolor='blue',linestyle = 'None')
ax.set_ylabel('Maximale Auslenkung')
ax.set_xlabel('Erregerkreisfrequenz $\Omega$')
ax.set_title('Sweep Duffing-Schwinger')
ax.legend(['sweep-up','sweep-down'],loc='upper right')
#
plt.show()
#
logfile.write('Programm finished\n\n\n')
logfile.close()