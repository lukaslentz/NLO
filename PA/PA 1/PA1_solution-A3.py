# #
# # mathematische Pendel
# #
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
# #
g = 9.81
l = 1
# #
# returns a time set for numerical integration
# t_span :  tuple with start and end value
# t : numpy ndarray with values for time
# T : periodicity
def time_set(omega):
    T = 2*np.pi/omega
    t_end = 15*T
    delta = T/499
    t_span = (0, t_end)
    t = np.arange(0, t_end, delta) # [0,t_end)
    return [t_span, t, T]
# 
# 
# returns the rhs for the penduulum as 1st order system
def pendel(t, state):
    [x, v] = state
    dx = v
    dv = -g/l*np.sin(x)
    return [dx, dv]
#
# returns the rhs for the linearized penduulum as 1st order system
def duffing(t, state):
    [x, v] = state
    dx = v
    dv = g/l*(-x + 1/6*x**3)
    return [dx, dv]
#
[t_span, t, T] = time_set(np.sqrt(g/l)) # get timeset
#
def v_event(t, state):
    return state[1]
v_event.direction = -1 
#
x0s = [0.1,0.2,0.3,0.5,0.8,1,1.2,1.5,2,2*np.pi-5]
freqs = [0]*len(x0s)
freqs_duff = [0]*len(x0s)
freqs_ref = [0]*len(x0s)

om0=np.sqrt(g/l)
mu=-om0**2/6

for i,x0 in enumerate(x0s):
    results_ivp_pendel = solve_ivp(pendel, t_span, [x0,0], t_eval = t, events = v_event) # numerical integration of penduulum
    t_event_pendel = results_ivp_pendel.t_events[0] # event times
    freqs[i] = 2*np.pi/np.mean(np.diff(t_event_pendel))
    #
    results_ivp_duff = solve_ivp(duffing, t_span, [x0,0], t_eval = t, events = v_event) # numerical integration of penduulum
    t_event_duff = results_ivp_duff.t_events[0] # event times
    freqs_duff[i] = 2*np.pi/np.mean(np.diff(t_event_duff))
    #
    freqs_ref[i] = om0*np.sqrt(1+3/4*mu/om0**2*x0**2)

# 
fig = plt.figure()

plt.plot(x0s, freqs, 'o',label='Pendel')
plt.plot(x0s, freqs_duff, 'o',label='Duffing')
plt.plot(x0s, freqs_ref, 'o',label='Näherung')
plt.legend()

plt.show()
# 