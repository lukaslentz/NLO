
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter

reds = [
    "lightcoral", "salmon", "darksalmon", "indianred", "crimson",
    "firebrick", "darkred", "red", "tomato", "orangered"
]

blues = [
    "lightblue", "skyblue", "deepskyblue", "dodgerblue", "cornflowerblue",
    "steelblue", "royalblue", "blue", "mediumblue", "darkblue", "navy"
]

g=9.81
l=1
m=20

om0 = np.sqrt(g/l)
theta = m*l**2

def solu(range,E0):
    return om0*np.sqrt(2*(E0+np.cos(range)))

x = np.linspace(-2.1*np.pi,2.1*np.pi,2999)
v0s = np.arange(0,8.1,1)
E0s = [1/2*(v0/om0)**2-1 for v0 in v0s]
E0s.append(1)

plt.figure()
for i,E0 in enumerate(E0s):
    y = solu(x,E0)
    if E0 > 1:
        color = reds.pop()
    elif E0 == 1:
        color = "black"
    else:
        color = blues.pop()
    plt.plot(x,y,color=color,label=f'{E0:.2f}')
    plt.plot(x,-y,color=color)
plt.scatter([-np.pi,np.pi],[0,0], color='gray', edgecolor='black', s=100,label='instabiles Gleichgewicht')
plt.scatter([-2*np.pi,0,2*np.pi],[0,0,0], color='orange', edgecolor='black', s=100,label='stabiles Gleichgewicht')
plt.legend(title="$E_0/(\Theta\omega_0^2)$", loc="upper left", bbox_to_anchor=(1, 1))
plt.grid(True)
plt.gca().xaxis.set_major_locator(MultipleLocator(base=np.pi))
plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda val, pos: 
    f"{int(val / np.pi)}π" if val != 0 else "0"))
plt.xlim(-2.1*np.pi,2.1*np.pi)
plt.ylim(-8,8)
plt.xlabel("Winkel x")
plt.ylabel("Winkelgeschwindigkeit $\dot{x}$")
plt.title("Phasenportrait Pendel")
plt.show()