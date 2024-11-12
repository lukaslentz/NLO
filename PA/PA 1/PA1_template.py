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
    t_end = 5*T
    delta = T/99
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
def pendel_lin(t, state):
    [x, v] = state
    dx = v
    dv = -g/l*x
    return [dx, dv]
#
[t_span, t, T] = time_set(np.sqrt(g/l)) # get timeset
# 
[x0,v0] = [0,15] # initial values
# 
results_ivp_pendel = solve_ivp(pendel, t_span, [x0,v0], t_eval = t) # numerical integration of penduulum
x_pendel = results_ivp_pendel.y[0] # displacement
v_pendel = results_ivp_pendel.y[1] # velocity
# 
#
results_ivp_lin = solve_ivp(pendel_lin, t_span, [x0,v0], t_eval = t) # numerical integration of linearized penduulum
x_lin = results_ivp_lin.y[0] # displacement
v_lin = results_ivp_lin.y[1] # velocity
#
# plot results
fig = plt.figure()
#
plt.plot(x_pendel, v_pendel)
plt.plot(x_lin, v_lin)
#
plt.show()
# #