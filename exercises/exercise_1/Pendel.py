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
# returns a time set for numerical integration
# t_span :  tuple wit start and end value
def time_set(omega):
    T = 2*np.pi/omega
    t_end = 5*T
    delta = T/100
    t_span = (0.0, t_end)
    t = np.arange(0.0, t_end, delta)
    return [t_span, t, T]
#
# returns the rhs for the penduulum as 1st order system
def pendel(t, state):
    x, v = state
    dx = v
    dv = -g/l*np.sin(x) 
    return [dx, dv]
#
# returns the rhs for the linearized penduulum as 1st order system
def pendel_lin(t, state):
    x, v = state
    dx = v
    dv = -g/l*x
    return [dx, dv]
#
[t_span, t, T] = time_set(np.sqrt(g/l)) # get timeset
[x0,v0] = [1,0] # initial values
#
results_ivp_pendel = solve_ivp(pendel, t_span, [x0,v0], method='RK45', t_eval=t, max_step = 0.01) # numerical integration of penduulum
x_p = results_ivp_pendel.y[0] # displacement
v_p = results_ivp_pendel.y[1] # velocity
#
results_ivp_lin = solve_ivp(pendel_lin, t_span, [x0,v0], method='RK45', t_eval=t, max_step = 0.01) # numerical integration of linearized penduulum
x_lin = results_ivp_lin.y[0] # displacement
v_lin = results_ivp_lin.y[1] # velocity
#
# plot results
fig = plt.figure()
#
ax = fig.add_subplot(1, 1, 1)
ax.plot(x_p,v_p,color='blue')
ax.plot(x_lin,v_lin,color='red')
ax.set_xlabel("Auslenkung")
ax.set_ylabel("Geschwindigkeit")
ax.legend(['Pendel','Linearisierung'],loc='upper right')
ax.set_title(f'Lösung für Anfangsbedingungen x0 = {x0:.2f} und v0 = {v0:.2f}')
#
plt.show()
#