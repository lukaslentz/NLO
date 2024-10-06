#
# mathematische Pendel
#
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
#
g = 9.81
l = 1
#
def time_set(Omega):
    T = 2*np.pi/Omega
    t_end = 5*T
    delta = T/100
    t_span = (0.0, t_end)
    t = np.arange(0.0, t_end, delta)
    return [t_span, t, T]
#
def pendel(t, state):
    x, v = state
    dx = v
    dv = -0.*v-g/l*np.sin(x) 
    return [dx, dv]
#
def v_event(t, state):
    return state[1]
v_event.direction = -1
#
[t_span, t, T] = time_set(np.sqrt(g/l))
[x0,v0] = [0.2,0]
#
results_ivp_pendel = solve_ivp(pendel, t_span, [x0,v0], method='RK45', t_eval=t, max_step = 0.01, events = v_event )
t_event = results_ivp_pendel.t_events[0]
x_event = results_ivp_pendel.y_events[0][:,0]
#
Ts = np.diff(t_event)
T = Ts.mean()
freq = 1/T
print(T)
print(freq)
#
fig = plt.figure()
#
ax = fig.add_subplot(1, 1, 1)
for pt in t_event:
    plt.axvline(pt, lw=0.5, color='black', alpha=0.5)
ax.plot(t,results_ivp_pendel.y[0],color='blue')
ax.plot(t_event,x_event,marker = 'o',markerfacecolor='red',markeredgecolor='red',linestyle = 'None')
ax.set_xlabel("Zeit")
ax.set_ylabel("Auslenkung")
ax.set_title('Wendepunkte')
#
plt.show()
#
