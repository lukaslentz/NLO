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
# t_span :  tuple with start and end value
# t : numpy ndarray with values for time
# T : periodicity
def time_set(omega):
    T = ...
    t_end = ...
    delta = ...
    t_span = ...
    t = ...
    return [t_span, t, T]
#
#
# returns the rhs for the penduulum as 1st order system
def pendel(t, state):
    [x, v] = ...
    dx = ...
    dv = ...
    return [dx, dv]
# 
# returns the rhs for the linearized penduulum as 1st order system
def pendel_lin(t, state):
    [x, v] = ...
    dx = ...
    dv = ...
    return [dx, dv]
#
[t_span, t, T] = time_set(np.sqrt(g/l)) # get timeset
[x0,v0] = [.5,0] # initial values
# #
results_ivp_pendel = solve_ivp(...) # numerical integration of penduulum
x_p = ... # displacement
v_p = ... # velocity
#
results_ivp_lin = solve_ivp(...) # numerical integration of linearized penduulum
x_lin = ... # displacement
v_lin = ... # velocity
#
# plot results
fig = plt.figure()
#
ax = fig.add_subplot(1, 1, 1)
ax.plot(t,x_p,color='blue')
# ax.plot(x_lin,v_lin,color='red')
ax.set_xlabel("Displacement")
ax.set_ylabel("Velocity")
ax.legend(['Full Equation','Linearization'],loc='upper right')
ax.set_title(f'Solution for Initial Conditions x0 = {x0:.2f} and v0 = {v0:.2f}')
#
plt.show()
#