import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

xi, alpha, gamma, f = 0.1, 1, -0.1, 0.2

def hb_equations():
    Omega, t = sp.symbols('Omega t')
    a, b = sp.symbols('a b')
    x_hb = a*sp.cos(Omega*t) + b*sp.sin(Omega*t)

    exp_res = (sp.diff(x_hb, t, 2) + xi*sp.diff(x_hb, t, 1) + alpha*x_hb + gamma*x_hb**3 - f*sp.cos(Omega*t)).rewrite(sp.cos,sp.exp).rewrite(sp.sin,sp.exp).expand()
    res = exp_res.rewrite(sp.exp,sp.cos).rewrite(sp.exp,sp.sin).simplify()

    eq_1 = res.coeff(sp.cos(Omega*t)) 
    eq_2 = res.coeff(sp.sin(Omega*t)) 

    return eq_1, eq_2, x_hb


def solve_hb_equations(eq_1, eq_2, x_hb, cur_Omega):
    Omega, t = sp.symbols('Omega t')
    a, b = sp.symbols('a b')

    cur_eq_1 = eq_1.subs({Omega:cur_Omega})
    cur_eq_2 = eq_2.subs({Omega:cur_Omega})

    hb_sol_full = sp.solve([cur_eq_1, cur_eq_2],[a,b])

    hb_sol = []
    eps = 1e-19

    for sol in hb_sol_full:
            if abs(sp.im(sol[0])) < eps and abs(sp.im(sol[1])) < eps:
                 hb_sol.append([sp.re(sol[0]), sp.re(sol[1])]) 

    return [sp.lambdify(t, x_hb.subs({a:sol[0], b:sol[1], Omega:cur_Omega})) for sol in hb_sol]


def plot_time_series(funcs, t_0, t_end):
     
    t_vals = np.linspace(t_0, t_end, 1000)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for func in funcs:
        ax.plot(t_vals, func(t_vals))
    ax.set_xlabel("Zeit $t$")
    ax.set_ylabel("Amplitude $C$")

    plt.show()

def main():

    cur_Omega = 0.8
    cur_T = 2*np.pi/cur_Omega
    eq_1, eq_2, x_hb = hb_equations()

    cur_sol = solve_hb_equations(eq_1, eq_2, x_hb, cur_Omega)

    plot_time_series(cur_sol,0,2*cur_T)


if __name__ == "__main__":
     main()