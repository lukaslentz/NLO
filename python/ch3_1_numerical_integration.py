"""
Chapter 3.1 -- Numerical Integration of Nonlinear Oscillators

This script compares three explicit one-step integrators on a *conservative*
nonlinear oscillator: the explicit (forward) Euler method, the symplectic
(semi-implicit) Euler method, and the classical fourth-order Runge-Kutta
method RK4.  It measures their convergence order, their long-term energy
behaviour, and it illustrates adaptive step size control with
scipy.integrate.solve_ivp on the van der Pol relaxation oscillator.

Test problem (plain text) -- the undamped pendulum of chapter 2.3:

    theta'' + omega0^2 * sin(theta) = 0 ,   state z = (theta, theta')

    z' = f(z) = ( z2 , -omega0^2 * sin(z1) )

It is used deliberately: it is genuinely nonlinear, it conserves the energy

    E = 0.5 * theta'^2 + omega0^2 * (1 - cos(theta)) ,

and its solution is known in closed form through Jacobi elliptic functions,

    theta(t)  = 2 * arcsin( k * sn(omega0 t + K(k), k) ) ,   k = sin(theta0/2),
    theta'(t) = 2 * omega0 * k * cn(omega0 t + K(k), k) ,

the phase shift K(k) placing the start at rest in the turning point,

so the *global* error can be measured against an exact reference instead of
against another numerical solution.

Methods (update from step n to n+1, h = step size):
    explicit Euler     z_{n+1} = z_n + h f(z_n)                      order 1
    symplectic Euler   v_{n+1} = v_n + h a(x_n)                      order 1
                       x_{n+1} = x_n + h v_{n+1}
    RK4                the usual four-stage weighted average          order 4

What the figure shows
---------------------
 (a) Global error (largest state error on the interval) versus step size on
     log-log axes.  The slope is the convergence order; the fitted values are
     printed and annotated (about 1, 1 and 4).
 (b) Relative energy error versus time at a fixed step size.  Explicit Euler
     grows exponentially (each step multiplies the energy of a harmonic
     oscillator by 1 + h^2 omega0^2), symplectic Euler oscillates within a
     bounded band without drifting, RK4 drifts slowly but stays tiny.
 (c) The same statement in the phase plane: the Euler orbit spirals outwards,
     the symplectic orbit stays on a closed curve.
 (d) Adaptive RK45 on the van der Pol oscillator: the accepted step size drops
     by more than an order of magnitude at every fast transition.
 (e) The absolute stability regions |R(h*lambda)| <= 1 of Euler and RK4.  An
     undamped oscillator sits on the imaginary axis, which Euler's unit disk
     only touches at the origin, while RK4 contains it up to |h omega0| = 2.83.
 (f) The stroboscopic Poincare section of the damped, forced Duffing
     oscillator (the example of the slide listing): the transient spirals in
     and the steady state collapses onto a single point, i.e. a period-T
     response.

Runtime is a few seconds.  The expensive part is the convergence study; making
the step-size list longer or its smallest step much smaller (STEPS_PER_PERIOD)
increases the runtime roughly in proportion to the number of steps taken.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import ellipk, ellipj

# =====================================================================
# PARAMETERS -- change these to experiment
# =====================================================================
OMEGA0 = 1.0             # natural frequency of the pendulum [rad/s]
THETA0_DEG = 60.0        # libration amplitude of the reference orbit [deg]

# Convergence study: number of steps per exact period, from coarse to fine.
# The finest value dominates the runtime (about 5*20000 steps per method).
STEPS_PER_PERIOD = np.array([20, 40, 80, 160, 320, 640, 1280,
                             2560, 5120, 10240, 20480])
N_PERIODS_ERROR = 5.0    # length of the integration interval, in periods

# Energy study: fixed step size and a long interval.
STEPS_PER_PERIOD_ENERGY = 50
N_PERIODS_ENERGY = 200.0

# Phase-plane illustration: deliberately coarse steps, few periods.
STEPS_PER_PERIOD_PHASE = 50
N_PERIODS_PHASE = 8.0
PHASE_LIMIT = 3.0        # axis range of the phase-plane panel [rad, rad/s]

# Van der Pol demonstration of adaptive step size control:
# x'' - mu (1 - x^2) x' + x = 0
VDP_MU = 5.0
VDP_T_END = 40.0
VDP_RTOL = 1.0e-8
VDP_ATOL = 1.0e-10

# Stroboscopic Poincare section of the forced Duffing oscillator
# x'' + 2 D om0 x' + om0^2 (x + eps x^3) = fhat cos(Omega t)
DUF_D = 0.05
DUF_OM0 = 1.0
DUF_EPS = 1.0
DUF_FHAT = 0.3
DUF_OMEGA = 1.2
DUF_N_PERIODS = 600      # excitation periods to integrate
DUF_N_TRANSIENT = 100    # periods discarded before recording section points

THETA0 = np.radians(THETA0_DEG)
K_MOD = np.sin(0.5 * THETA0)                       # elliptic modulus
U_QUARTER = ellipk(K_MOD ** 2)                     # K(k); quarter period in u
T_EXACT = 4.0 / OMEGA0 * U_QUARTER                 # exact period (m = k^2!)


# =====================================================================
# Test problem: right-hand side, exact solution, energy
# =====================================================================
def accel(theta):
    """Acceleration of the undamped pendulum, theta'' = -omega0^2 sin(theta)."""
    return -OMEGA0 ** 2 * np.sin(theta)


def rhs(t, z):
    """State-space right-hand side f(z) = (z2, -omega0^2 sin(z1))."""
    return np.array([z[1], accel(z[0])])


def energy(theta, theta_dot):
    """Conserved energy E = 0.5 theta'^2 + omega0^2 (1 - cos theta)."""
    return 0.5 * theta_dot ** 2 + OMEGA0 ** 2 * (1.0 - np.cos(theta))


def exact_state(t):
    """
    Exact solution of the test problem at time t (array or scalar).

    The reference orbit is released from rest at the amplitude THETA0.  Since
    sn(K) = 1 and cn(K) = 0, shifting the argument by the quarter period K(k)
    puts the turning point at t = 0, exactly as in chapter 2.3.
    """
    sn, cn, _dn, _ph = ellipj(OMEGA0 * t + U_QUARTER, K_MOD ** 2)
    return np.array([2.0 * np.arcsin(K_MOD * sn), 2.0 * OMEGA0 * K_MOD * cn])


Z0 = exact_state(0.0)            # initial state: at rest at theta = theta0
E0 = energy(Z0[0], Z0[1])        # exact energy level of the reference orbit


# =====================================================================
# Integrators (fixed step size, written out explicitly for teaching)
# =====================================================================
def integrate_euler(z0, h, n_steps):
    """
    Explicit (forward) Euler:  z_{n+1} = z_n + h f(z_n).

    One function evaluation per step, global error O(h).  For an oscillator
    the amplification factor per step is |1 + i h omega0| = sqrt(1+h^2 omega0^2)
    > 1 for every h > 0, so the method is unconditionally unstable here -- the
    amplitude and the energy grow no matter how small the step is chosen.
    """
    z = np.empty((n_steps + 1, 2))
    z[0] = z0
    for n in range(n_steps):
        x, v = z[n]
        z[n + 1] = (x + h * v, v + h * accel(x))
    return z


def integrate_symplectic_euler(z0, h, n_steps):
    """
    Symplectic (semi-implicit) Euler:

        v_{n+1} = v_n + h a(x_n)
        x_{n+1} = x_n + h v_{n+1}     <-- uses the *new* velocity

    Still only first order and still one force evaluation per step, but the
    map is area preserving in the phase plane.  It therefore conserves a
    slightly perturbed ("shadow") energy exactly, which keeps the true energy
    error bounded and oscillatory for arbitrarily long times.
    """
    z = np.empty((n_steps + 1, 2))
    z[0] = z0
    x, v = z0
    for n in range(n_steps):
        v = v + h * accel(x)
        x = x + h * v
        z[n + 1] = (x, v)
    return z


def integrate_rk4(z0, h, n_steps):
    """
    Classical Runge-Kutta of fourth order: four evaluations per step,
    global error O(h^4).  Halving h reduces the error by a factor 16 while
    doubling the cost -- the reason RK4 beats Euler by orders of magnitude at
    equal work.  RK4 is not symplectic, so its energy still drifts, but the
    drift is proportional to h^4 and invisible over moderate intervals.
    """
    z = np.empty((n_steps + 1, 2))
    z[0] = z0
    for n in range(n_steps):
        y = z[n]
        k1 = rhs(0.0, y)
        k2 = rhs(0.0, y + 0.5 * h * k1)
        k3 = rhs(0.0, y + 0.5 * h * k2)
        k4 = rhs(0.0, y + h * k3)
        z[n + 1] = y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    return z


METHODS = [("explicit Euler", integrate_euler, "tab:red"),
           ("symplectic Euler", integrate_symplectic_euler, "tab:green"),
           ("RK4", integrate_rk4, "tab:blue")]


# =====================================================================
# Study 1: convergence order
# =====================================================================
def convergence_study():
    """
    Integrate the test problem over N_PERIODS_ERROR periods with several fixed
    step sizes and measure the global error as the largest state error that
    occurs anywhere on the interval,

        err(h) = max_n || z_num(t_n) - z_exact(t_n) ||_2 .

    The maximum, not the value at the final time, is the honest measure: the
    error of a nonlinear oscillator is dominated by a phase (frequency) error
    that grows with time, and sampling it only at t_end can accidentally hit a
    turning point, where a phase error hardly shows up in the position.

    Returns (h_values, {method name: error array}).
    """
    t_end = N_PERIODS_ERROR * T_EXACT
    h_values = t_end / (N_PERIODS_ERROR * STEPS_PER_PERIOD)
    errors = {}
    for name, method, _color in METHODS:
        err = []
        for n_per_period in STEPS_PER_PERIOD:
            n_steps = int(round(N_PERIODS_ERROR * n_per_period))
            h = t_end / n_steps
            z = method(Z0, h, n_steps)
            z_ref = exact_state(np.arange(n_steps + 1) * h)
            err.append(np.max(np.hypot(z[:, 0] - z_ref[0], z[:, 1] - z_ref[1])))
        errors[name] = np.array(err)
    return h_values, errors


def fitted_order(h_values, err, err_floor=1.0e-12, err_ceiling=1.0, n_fit=4):
    """
    Estimate the convergence order as the slope of log(err) over log(h).

    Only points that are neither saturated (error as large as the solution
    itself, i.e. no digits left) nor dominated by round-off are usable, and
    of those only the n_fit finest steps are fitted: the order statement is
    an asymptotic one for h -> 0, and coarse steps are still pre-asymptotic.
    """
    mask = (err > err_floor) & (err < err_ceiling)
    h_ok, e_ok = h_values[mask], err[mask]
    if h_ok.size < 2:
        return np.nan
    order = np.argsort(h_ok)                 # finest steps first
    idx = order[:min(n_fit, h_ok.size)]
    return np.polyfit(np.log(h_ok[idx]), np.log(e_ok[idx]), 1)[0]


# =====================================================================
# Study 2: energy behaviour
# =====================================================================
def energy_study():
    """
    Integrate over N_PERIODS_ENERGY periods with a fixed, moderate step size
    and return (t, {method name: relative energy error}).
    """
    h = T_EXACT / STEPS_PER_PERIOD_ENERGY
    n_steps = int(round(N_PERIODS_ENERGY * STEPS_PER_PERIOD_ENERGY))
    t = np.arange(n_steps + 1) * h
    drift = {}
    for name, method, _color in METHODS:
        z = method(Z0, h, n_steps)
        drift[name] = (energy(z[:, 0], z[:, 1]) - E0) / E0
    return t, drift, h


def phase_study():
    """Coarse-step orbits in the phase plane, for the qualitative picture."""
    h = T_EXACT / STEPS_PER_PERIOD_PHASE
    n_steps = int(round(N_PERIODS_PHASE * STEPS_PER_PERIOD_PHASE))
    return {name: method(Z0, h, n_steps) for name, method, _c in METHODS}, h


# =====================================================================
# Study 3: adaptive step size control (van der Pol, RK45)
# =====================================================================
def van_der_pol(t, z):
    """Van der Pol oscillator x'' - mu (1 - x^2) x' + x = 0 in state form."""
    return [z[1], VDP_MU * (1.0 - z[0] ** 2) * z[1] - z[0]]


def adaptive_study():
    """
    Solve the van der Pol oscillator with the adaptive RK45 pair and return
    the accepted step points.  sol.t contains exactly the accepted steps, so
    np.diff(sol.t) is the sequence of step sizes actually used.
    """
    sol = solve_ivp(van_der_pol, (0.0, VDP_T_END), [2.0, 0.0], method="RK45",
                    rtol=VDP_RTOL, atol=VDP_ATOL, dense_output=True)
    if not sol.success:
        raise RuntimeError("van der Pol integration failed: " + sol.message)
    return sol


# =====================================================================
# Study 4: absolute stability regions
# =====================================================================
def stability_functions(z):
    """
    Amplification factors R(h*lambda) of Euler and RK4 for the scalar test
    equation z' = lambda z.  A method is absolutely stable where |R| <= 1.

    Euler: R = 1 + z.
    RK4:   R = 1 + z + z^2/2 + z^3/6 + z^4/24 (the Taylor series up to h^4).
    """
    r_euler = 1.0 + z
    r_rk4 = 1.0 + z + z ** 2 / 2.0 + z ** 3 / 6.0 + z ** 4 / 24.0
    return r_euler, r_rk4


# =====================================================================
# Study 5: stroboscopic Poincare section of the forced Duffing oscillator
# =====================================================================
def duffing(t, z):
    """Forced, damped Duffing oscillator in state-space form."""
    return [z[1],
            -2.0 * DUF_D * DUF_OM0 * z[1]
            - DUF_OM0 ** 2 * (z[0] + DUF_EPS * z[0] ** 3)
            + DUF_FHAT * np.cos(DUF_OMEGA * t)]


def poincare_study():
    """
    Integrate the forced Duffing oscillator and sample the state once per
    excitation period T = 2 pi / Omega (stroboscopic map).

    dense_output=True lets us evaluate the interpolant exactly at t_n = n T:
    the adaptive solver never steps onto those instants by itself, and
    sampling the nearest stored point instead would blur the section.
    The first DUF_N_TRANSIENT periods are discarded as transient.
    """
    t_exc = 2.0 * np.pi / DUF_OMEGA
    sol = solve_ivp(duffing, (0.0, DUF_N_PERIODS * t_exc), [0.0, 0.0],
                    method="RK45", rtol=1.0e-9, atol=1.0e-11,
                    max_step=t_exc / 50.0, dense_output=True)
    if not sol.success:
        raise RuntimeError("Duffing integration failed: " + sol.message)
    t_poin = np.arange(DUF_N_TRANSIENT, DUF_N_PERIODS) * t_exc
    t_trans = np.arange(0, DUF_N_TRANSIENT) * t_exc
    return sol.sol(t_poin), sol.sol(t_trans), t_exc


# =====================================================================
# Console output
# =====================================================================
def print_results(h_values, errors, t_energy, drift, h_energy, sol,
                  poincare, t_exc):
    """Print the quantitative results behind the six panels."""
    print(f"Test problem: pendulum, omega0 = {OMEGA0}, "
          f"theta0 = {THETA0_DEG:.0f} deg")
    print(f"exact period T = {T_EXACT:.6f} s, energy E0 = {E0:.6f}")
    print()

    print("Global error over "
          f"{N_PERIODS_ERROR:.0f} periods (max Euclidean state error)")
    header = "  h        steps/period " + "".join(f"{n:>18s}" for n, _, _ in METHODS)
    print(header)
    print("-" * len(header))
    for i, h in enumerate(h_values):
        row = f"{h:9.2e} {STEPS_PER_PERIOD[i]:10d}   "
        row += "".join(f"{errors[n][i]:18.3e}" for n, _, _ in METHODS)
        print(row)
    print()
    print("measured convergence order (slope of the log-log fit):")
    for name, _m, _c in METHODS:
        print(f"  {name:18s} p = {fitted_order(h_values, errors[name]):.2f}")
    print()

    print(f"Energy behaviour over {N_PERIODS_ENERGY:.0f} periods "
          f"at h = T/{STEPS_PER_PERIOD_ENERGY} = {h_energy:.4f} s:")
    for name, _m, _c in METHODS:
        d = drift[name]
        print(f"  {name:18s} final (E-E0)/E0 = {d[-1]:12.4e}   "
              f"max |(E-E0)/E0| = {np.max(np.abs(d)):.4e}")
    # For the *linearised* oscillator explicit Euler multiplies the energy by
    # exactly (1 + h^2 omega0^2) per step, i.e. by the factor below per period.
    # That closed form explains the initial exponential rise of the red curve;
    # later the growth slows down only because the pendulum has left the
    # small-angle regime (the restoring term sin(theta) stays bounded).
    per_period = (1.0 + (h_energy * OMEGA0) ** 2) ** STEPS_PER_PERIOD_ENERGY
    print(f"  linear-theory Euler energy growth per period "
          f"(1+h^2 omega0^2)^(T/h) = {per_period:.3f}")
    print()

    dt = np.diff(sol.t)
    print(f"Adaptive RK45 on van der Pol (mu = {VDP_MU}): "
          f"{sol.t.size - 1} accepted steps, {sol.nfev} function evaluations")
    print(f"  step size: min {dt.min():.2e}, median {np.median(dt):.2e}, "
          f"max {dt.max():.2e}  (ratio {dt.max() / dt.min():.1f})")
    print(f"  RK4 stability limit on the imaginary axis: "
          f"h*omega0 <= 2*sqrt(2) = {2.0 * np.sqrt(2.0):.4f}; "
          f"explicit Euler: none")
    print()

    # The stroboscopic samples of a periodic steady state all collapse onto a
    # single point; the scatter below is therefore a measure of how close the
    # response is to a plain period-T oscillation.
    x_p, v_p = poincare
    print(f"Duffing Poincare section (Omega = {DUF_OMEGA}, "
          f"f_hat = {DUF_FHAT}, D = {DUF_D}): "
          f"{x_p.size} points, excitation period T = {t_exc:.4f} s")
    print(f"  mean point (x, x') = ({x_p.mean():.6f}, {v_p.mean():.6f})")
    print(f"  spread: std(x) = {x_p.std():.3e}, std(x') = {v_p.std():.3e} "
          f"-> {'period-T steady state' if x_p.std() < 1e-6 else 'not a single point'}")
    print()


# =====================================================================
# Figure
# =====================================================================
def make_figure(h_values, errors, t_energy, drift, h_energy, orbits, h_phase,
                sol, poincare, transient, t_exc):
    """Assemble the six-panel summary figure."""
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 8.8))
    (ax_conv, ax_energy, ax_phase, ax_adapt, ax_stab, ax_poin) = axes.flat

    # ---------------- (a) convergence order -------------------------
    for name, _m, color in METHODS:
        p = fitted_order(h_values, errors[name])
        ax_conv.loglog(h_values, errors[name], "o-", ms=4, color=color,
                       label=f"{name}  (fitted $p$ = {p:.2f})")
    # Reference slopes 1 and 4 for the eye.
    href = np.array([h_values[-1], h_values[0]])
    ax_conv.loglog(href, errors["explicit Euler"][-1] * (href / href[0]) ** 1,
                   "k:", lw=1.0)
    ax_conv.loglog(href, errors["RK4"][-1] * (href / href[0]) ** 4,
                   "k--", lw=1.0)
    ax_conv.text(href[0] * 3, errors["RK4"][-1] * 30, r"slope 4", fontsize=9)
    ax_conv.text(href[0] * 3, errors["explicit Euler"][-1] * 4, r"slope 1",
                 fontsize=9)
    ax_conv.set_xlabel(r"step size $h$ [s]")
    ax_conv.set_ylabel(r"global error $\max_n\|z_h(t_n)-z(t_n)\|$")
    ax_conv.set_title(f"(a) Global error over {N_PERIODS_ERROR:.0f} periods")
    ax_conv.legend(loc="lower right", fontsize=9)
    ax_conv.grid(True, which="both", alpha=0.3)

    # ---------------- (b) energy drift ------------------------------
    for name, _m, color in METHODS:
        ax_energy.semilogy(t_energy / T_EXACT, np.abs(drift[name]) + 1.0e-18,
                           color=color, lw=1.0, label=name)
    ax_energy.set_ylim(1.0e-8, 1.0e3)
    ax_energy.set_xlabel(r"time $t/T$ [periods]")
    ax_energy.set_ylabel(r"$|E-E_0|/E_0$")
    ax_energy.set_title(f"(b) Energy error at $h=T/{STEPS_PER_PERIOD_ENERGY}$"
                        f" = {h_energy:.3f} s")
    ax_energy.legend(loc="lower right", fontsize=9)
    ax_energy.grid(True, which="both", alpha=0.3)

    # ---------------- (c) phase plane -------------------------------
    t_ref = np.linspace(0.0, T_EXACT, 600)
    z_ref = exact_state(t_ref)
    ax_phase.plot(z_ref[0], z_ref[1], "k-", lw=1.5, label="exact orbit")
    for name, _m, color in METHODS:
        z = orbits[name]
        ax_phase.plot(z[:, 0], z[:, 1], color=color, lw=0.8, alpha=0.85,
                      label=name)
    # Fixed limits: the explicit Euler orbit gains so much energy that it
    # crosses the separatrix and runs off as a whirling motion -- it is meant
    # to leave the frame.
    ax_phase.set_xlim(-PHASE_LIMIT, PHASE_LIMIT)
    ax_phase.set_ylim(-PHASE_LIMIT, PHASE_LIMIT)
    ax_phase.set_xlabel(r"$\theta$ [rad]")
    ax_phase.set_ylabel(r"$\dot\theta$ [rad/s]")
    ax_phase.set_title(f"(c) {N_PERIODS_PHASE:.0f} periods at the coarse step "
                       f"$h=T/{STEPS_PER_PERIOD_PHASE}$")
    ax_phase.legend(loc="upper right", fontsize=9)
    ax_phase.grid(alpha=0.3)
    ax_phase.set_aspect("equal")

    # ---------------- (d) adaptive step size ------------------------
    t_dense = np.linspace(0.0, VDP_T_END, 4000)
    ax_adapt.plot(t_dense, sol.sol(t_dense)[0], color="tab:purple", lw=1.0,
                  label=r"$x(t)$, van der Pol")
    ax_adapt.set_xlabel("time $t$ [s]")
    ax_adapt.set_ylabel(r"$x$", color="tab:purple")
    ax_adapt.tick_params(axis="y", labelcolor="tab:purple")
    ax_adapt.set_title(f"(d) Adaptive RK45, van der Pol $\\mu={VDP_MU:.0f}$: "
                       f"{sol.t.size - 1} accepted steps")
    ax_adapt.grid(alpha=0.3)

    ax_h = ax_adapt.twinx()
    # Each accepted step is plotted at its own start time; the step size
    # collapses exactly where the relaxation oscillation jumps.
    ax_h.semilogy(sol.t[:-1], np.diff(sol.t), ".", ms=3, color="tab:orange")
    ax_h.set_ylabel("accepted step size $h$ [s]", color="tab:orange")
    ax_h.tick_params(axis="y", labelcolor="tab:orange")

    # ---------------- (e) absolute stability regions ----------------
    re = np.linspace(-3.5, 1.5, 500)
    im = np.linspace(-3.5, 3.5, 700)
    zz = re[None, :] + 1j * im[:, None]
    r_euler, r_rk4 = stability_functions(zz)
    ax_stab.contourf(re, im, np.abs(r_rk4), levels=[0.0, 1.0],
                     colors=["tab:blue"], alpha=0.30)
    ax_stab.contour(re, im, np.abs(r_rk4), levels=[1.0], colors="tab:blue")
    ax_stab.contourf(re, im, np.abs(r_euler), levels=[0.0, 1.0],
                     colors=["tab:red"], alpha=0.30)
    ax_stab.contour(re, im, np.abs(r_euler), levels=[1.0], colors="tab:red")
    # An undamped oscillator has h*lambda = +-i h omega0, i.e. it lives on the
    # imaginary axis: Euler's disk only touches it at the origin, whereas RK4
    # contains it up to |h omega0| = 2 sqrt(2) = 2.83.
    ax_stab.axvline(0.0, color="k", lw=0.8)
    ax_stab.plot([0.0, 0.0], [-2.0 * np.sqrt(2.0), 2.0 * np.sqrt(2.0)],
                 "k-", lw=2.5, alpha=0.6)
    ax_stab.annotate(r"RK4 covers $|h\omega_0|\leq 2\sqrt{2}$",
                     xy=(-3.3, -3.2), fontsize=9)
    ax_stab.plot([], [], color="tab:red", label="explicit Euler")
    ax_stab.plot([], [], color="tab:blue", label="RK4")
    ax_stab.set_xlabel(r"$\mathrm{Re}(h\lambda)$")
    ax_stab.set_ylabel(r"$\mathrm{Im}(h\lambda)$")
    ax_stab.set_title("(e) Absolute stability regions")
    ax_stab.legend(loc="upper left", fontsize=9)
    ax_stab.set_aspect("equal")
    ax_stab.grid(alpha=0.3)

    # ---------------- (f) stroboscopic Poincare section -------------
    # The transient samples (light) spiral into the steady-state point; after
    # the transient all 500 stroboscopic samples fall on one single point,
    # which is the signature of a period-T (harmonic) steady state.
    ax_poin.plot(transient[0], transient[1], ".", ms=3, color="0.7",
                 label=f"first {DUF_N_TRANSIENT} periods (transient)")
    ax_poin.plot(poincare[0], poincare[1], "o", ms=6, color="steelblue",
                 label=f"periods {DUF_N_TRANSIENT}-{DUF_N_PERIODS}: one point")
    ax_poin.legend(loc="lower left", fontsize=8)
    ax_poin.set_xlabel(r"$x(nT)$")
    ax_poin.set_ylabel(r"$\dot x(nT)$")
    ax_poin.set_title(f"(f) Duffing Poincare section, "
                      f"$\\Omega={DUF_OMEGA}$, $\\hat f={DUF_FHAT}$")
    ax_poin.grid(alpha=0.3)

    fig.suptitle("3.1 Numerical integration: convergence order, energy "
                 "behaviour, stability, adaptive steps and Poincare sections")
    fig.tight_layout()


def main():
    h_values, errors = convergence_study()
    t_energy, drift, h_energy = energy_study()
    orbits, h_phase = phase_study()
    sol = adaptive_study()
    poincare, transient, t_exc = poincare_study()

    print_results(h_values, errors, t_energy, drift, h_energy, sol,
                  poincare, t_exc)
    make_figure(h_values, errors, t_energy, drift, h_energy, orbits, h_phase,
                sol, poincare, transient, t_exc)
    plt.show()


if __name__ == "__main__":
    main()
