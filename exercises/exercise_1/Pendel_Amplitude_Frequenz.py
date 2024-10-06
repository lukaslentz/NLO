#
# mathematische Pendel
#
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
#
omega0 = 1
mu = -1/6*omega0**2
#
def time_set(omega):
    T = 2*np.pi/omega
    t_end = 5*T
    delta = T/100
    t_span = (0.0, t_end)
    t = np.arange(0.0, t_end, delta)
    return [t_span, t, T]
#
def pendel(t, state):
    x, v = state
    dx = v
    dv = -omega0**2*np.sin(x) 
    return [dx, dv]
#
def duffing(t, state):
    x, v = state
    dx = v
    dv = -omega0**2*x-mu*x**3 
    return [dx, dv]
#
def v_event(t, state):
    return state[1]
v_event.direction = -1
#
[t_span, t, T] = time_set(omega0)
#
#x0s = [0.1,0.2,0.3,0.5,0.8,1.0,1.2,1.5,2.0,5.0,2*np.pi-5]
x0s = np.arange(0.1,2.1,0.1)
freqs = [0]*len(x0s) 
#
refs = [ omega0 *np.sqrt(1+3/4*mu/omega0**2*x0**2)for x0 in x0s]
#
for i,x0 in enumerate(x0s):
#
    results_ivp_pendel = solve_ivp(duffing, t_span, [x0,0], method='RK45', t_eval=t, max_step = 0.01, events = v_event )
    t_event = results_ivp_pendel.t_events[0]
    x_event = results_ivp_pendel.y_events[0][:,0]
    #
    Ts = np.diff(t_event)
    T = Ts.mean()
    freqs[i] = 2*np.pi/T
#
fig = plt.figure()
#
ax = fig.add_subplot(1, 1, 1)
ax.plot(refs,x0s,color='blue')
ax.plot(freqs,x0s,marker = 'o',markerfacecolor='red',markeredgecolor='red',linestyle = 'None')
ax.set_ylabel("maximale Auslenkung")
ax.set_xlabel("Frequenz")
ax.set_title('Abhängigkeit der Frequenz von der Auslenkung')
#
plt.show()
#
# write data to file
import csv
# open the file in the write mode
with open('C:/SeafileContainer/Seafile/Meine Bibliothek/Lehre/NLO/Programme/Aufgabe 1/freq.csv', 'w', newline='') as f:
    # create the csv writer
    writer = csv.writer(f)
    writer.writerow(['amp\tfreq'])
    # write a row to the csv file
    for i, freq in enumerate(freqs):
        writer.writerow([str(x0s[i])+'\t'+str(freq)])