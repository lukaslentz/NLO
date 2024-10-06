#
# Duffing angetrieben
#
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
#
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
xi=0.1
alpha= 1
gamma=-0.1
f=0.2
#
Omegas_up = np.arange(0.3,1.3,0.1)
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
for Omega in tqdm(Omegas_up):
    T = 2*np.pi/Omega
    results_ivp_duffing = solve_ivp(duffing, [t0,t0+20*T], [x0,v0], max_step = 0.01, events = v_event, args = [xi, alpha, gamma, f, Omega])
    if results_ivp_duffing.success == True:
        t0 = (results_ivp_duffing.t_events[0][0])%T
        [x0,v0] = results_ivp_duffing.y_events[0][-1]
        x_val_up.append(Omega)
        y_val_up.append(x0)
#
for Omega in tqdm(Omegas_down):
    T = 2*np.pi/Omega
    results_ivp_duffing = solve_ivp(duffing, [t0,t0+20*T], [x0,v0], max_step = 0.01, events = v_event, args = [xi, alpha, gamma, f, Omega])
    if results_ivp_duffing.success == True:
        t0 = (results_ivp_duffing.t_events[0][0])%T
        [x0,v0] = results_ivp_duffing.y_events[0][-1]
        x_val_down.append(Omega)
        y_val_down.append(x0)
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