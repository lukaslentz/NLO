"""
Chapter 1.1 -- Nonlinear Systems

Two hallmark phenomena of nonlinear oscillators that have no counterpart in a
linear model are computed and displayed in one figure:

1. Self-excited oscillation and limit cycle (van der Pol oscillator).
   Governing equation (plain text):

       x'' - mu * (1 - x**2) * x' + x = 0 ,      mu > 0

   For |x| < 1 the effective damping coefficient -mu*(1-x^2) is negative, so
   energy is fed into the system; for |x| > 1 it is positive and energy is
   dissipated.  The balance of the two produces an isolated closed orbit, the
   limit cycle.  Trajectories started inside spiral outward, trajectories
   started outside spiral inward, and all of them end on the same orbit: the
   steady-state amplitude is independent of the initial conditions.
   The figure shows the time histories x(t) for three initial conditions and
   the corresponding phase portrait (x, dx/dt) on top of the vector field.

2. Amplitude-dependent natural frequency (backbone curve) of the free,
   undamped hardening Duffing oscillator:

       x'' + x + eps * x**3 = 0 ,                eps > 0 (hardening)

   The Lindstedt-Poincare perturbation method gives the first-order result
   omega(A) = 1 + (3/8)*eps*A**2.  The script measures omega(A) directly from
   numerical time integration (period taken from successive zero crossings of
   the velocity, located with a root-finding event) and compares it with the
   perturbation result and with the exact period expressed through the complete
   elliptic integral of the first kind.  The linear oscillator would give the
   vertical line omega = 1.

Only numpy, scipy and matplotlib are used.  Total runtime is a few seconds.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import ellipk

# =====================================================================
# PARAMETERS -- change these to experiment
# =====================================================================

# --- van der Pol oscillator -------------------------------------------------
MU = 1.0                       # nonlinear damping parameter mu (slides: mu = 1)
X0_LIST = [0.2, 1.0, 2.5]      # three initial displacements, all with v0 = 0
T_END = 40.0                   # integration time (long enough to settle)
N_SAMPLES = 4000               # number of output samples for the time history

# --- free Duffing oscillator (backbone curve) -------------------------------
EPS = 0.5                      # cubic coefficient eps > 0 -> hardening spring
A_MIN, A_MAX = 0.05, 2.0       # amplitude range of the backbone curve
N_AMPLITUDES = 40              # number of amplitudes evaluated numerically

# --- numerics ---------------------------------------------------------------
RTOL, ATOL = 1e-10, 1e-12      # tight tolerances: the period must be accurate


# =====================================================================
# VAN DER POL OSCILLATOR
# =====================================================================

def van_der_pol(t, y, mu):
    """Right-hand side of the van der Pol oscillator in state-space form.

    State vector y = [x, v] with v = dx/dt, so that

        dx/dt = v
        dv/dt = mu * (1 - x**2) * v - x

    The damping term changes sign at |x| = 1, which is the mechanism that
    creates the limit cycle.
    """
    x, v = y
    return [v, mu * (1.0 - x**2) * v - x]


def integrate_vdp(x0, mu=MU, t_end=T_END, n=N_SAMPLES):
    """Integrate the van der Pol oscillator from (x0, 0) and return (t, x, v)."""
    t_eval = np.linspace(0.0, t_end, n)
    sol = solve_ivp(van_der_pol, (0.0, t_end), [x0, 0.0], args=(mu,),
                    t_eval=t_eval, rtol=RTOL, atol=ATOL, dense_output=False)
    return sol.t, sol.y[0], sol.y[1]


def limit_cycle_amplitude(x, v, t, t_settle):
    """Peak displacement of the (already settled) limit cycle.

    Only the tail of the record, t > t_settle, is used so that the initial
    transient does not contaminate the amplitude estimate.
    """
    mask = t > t_settle
    return np.max(np.abs(x[mask]))


# =====================================================================
# FREE DUFFING OSCILLATOR -- BACKBONE CURVE omega(A)
# =====================================================================

def duffing_free(t, y, eps):
    """Right-hand side of the free, undamped Duffing oscillator.

        dx/dt = v
        dv/dt = -x - eps * x**3

    There is no damping and no forcing, so the total energy
    E = v^2/2 + x^2/2 + eps*x^4/4 is conserved and every solution is periodic.
    """
    x, v = y
    return [v, -x - eps * x**3]


def _velocity_zero(t, y, eps):
    """Event function: zero crossing of the velocity v = y[1]."""
    return y[1]


# The motion starts at a turning point (v = 0).  Half a period later the
# velocity vanishes again at the opposite turning point, so the spacing of two
# successive velocity zeros is exactly half a period.  Direction is left free,
# so every zero of the velocity is recorded.
_velocity_zero.direction = 0.0
_velocity_zero.terminal = False


def duffing_frequency_numeric(amplitude, eps=EPS):
    """Measure the oscillation frequency of the free Duffing oscillator.

    The oscillator is released from rest at x(0) = A, i.e. at a turning point.
    Successive zeros of the velocity are half a period apart; the event root
    finder of solve_ivp locates them to integration accuracy, which is far more
    accurate than reading peaks off a sampled time series.
    """
    # Rough linear guess for the period, generously extended so that at least
    # two velocity zeros are certain to be contained in the interval.
    t_end = 3.0 * 2.0 * np.pi
    sol = solve_ivp(duffing_free, (0.0, t_end), [amplitude, 0.0], args=(eps,),
                    events=_velocity_zero, rtol=RTOL, atol=ATOL)
    zeros = sol.t_events[0]
    zeros = zeros[zeros > 1e-9]          # discard the trivial zero at t = 0
    if len(zeros) < 2:
        raise RuntimeError("velocity zeros not found -- increase t_end")
    period = 2.0 * (zeros[1] - zeros[0])  # one half period -> full period
    return 2.0 * np.pi / period


def duffing_frequency_lindstedt(amplitude, eps=EPS):
    """First-order Lindstedt-Poincare result omega(A) = 1 + (3/8)*eps*A^2."""
    return 1.0 + 0.375 * eps * amplitude**2


def duffing_frequency_exact(amplitude, eps=EPS):
    """Exact frequency of x'' + x + eps*x^3 = 0 via elliptic integrals.

    Energy conservation gives a quadrature for the quarter period which
    evaluates to the complete elliptic integral of the first kind K(m):

        omega(A) = pi/2 * sqrt(1 + eps*A^2) / K(m),
        m = eps*A^2 / (2*(1 + eps*A^2)).

    scipy's ellipk takes the parameter m = k^2 (not the modulus k).
    """
    q = eps * amplitude**2
    m = q / (2.0 * (1.0 + q))
    return 0.5 * np.pi * np.sqrt(1.0 + q) / ellipk(m)


# =====================================================================
# PLOTTING
# =====================================================================

def plot_vdp_time_history(ax, results):
    """Time histories x(t) for the three initial conditions."""
    for x0, (t, x, _v) in results.items():
        ax.plot(t, x, lw=1.2, label=f"$x_0 = {x0}$")
    ax.set_xlabel("time $t$")
    ax.set_ylabel("displacement $x(t)$")
    ax.set_title(f"van der Pol time history ($\\mu = {MU}$)")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.3)


def plot_vdp_phase_portrait(ax, results):
    """Phase portrait with the vector field of the van der Pol oscillator."""
    # Background vector field (normalised arrows: only the direction matters).
    xg, vg = np.meshgrid(np.linspace(-3.0, 3.0, 21), np.linspace(-4.0, 4.0, 21))
    dx = vg
    dv = MU * (1.0 - xg**2) * vg - xg
    speed = np.hypot(dx, dv)
    speed_safe = np.where(speed > 0.0, speed, 1.0)   # avoid 0/0 at the equilibrium
    dx, dv = dx / speed_safe, dv / speed_safe
    ax.quiver(xg, vg, dx, dv, color="0.65", pivot="mid", scale=32,
              width=0.003)

    for x0, (_t, x, v) in results.items():
        ax.plot(x, v, lw=1.2, label=f"$x_0 = {x0}$")
        ax.plot(x[0], v[0], "o", ms=4, color="k")

    ax.plot(0.0, 0.0, "x", color="crimson", ms=8, mew=2,
            label="unstable equilibrium")
    ax.set_xlabel("$x$")
    ax.set_ylabel("$\\dot{x}$")
    ax.set_title("Phase portrait: all orbits reach one limit cycle")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.3)


def plot_backbone(ax, amps, omega_num, omega_lp, omega_ex):
    """Backbone curve: oscillation frequency versus amplitude."""
    ax.plot(omega_num, amps, "o", ms=4, color="tab:blue",
            label="numerical integration")
    ax.plot(omega_ex, amps, "-", lw=1.6, color="tab:green",
            label="exact (elliptic integral)")
    ax.plot(omega_lp, amps, "--", lw=1.6, color="tab:red",
            label="Lindstedt-Poincare  $1 + \\frac{3}{8}\\varepsilon A^2$")
    ax.axvline(1.0, color="k", lw=1.2, ls=":",
               label="linear oscillator $\\omega_0 = 1$")
    ax.set_xlabel("oscillation frequency $\\omega$")
    ax.set_ylabel("amplitude $A$")
    ax.set_title(f"Backbone curve of $\\ddot{{x}} + x + \\varepsilon x^3 = 0$ "
                 f"($\\varepsilon = {EPS}$, hardening)")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(alpha=0.3)


# =====================================================================
# MAIN
# =====================================================================

def main():
    # ---------------- van der Pol limit cycle ----------------
    results = {x0: integrate_vdp(x0) for x0 in X0_LIST}

    print("van der Pol oscillator, mu = %.2f" % MU)
    for x0, (t, x, v) in results.items():
        amp = limit_cycle_amplitude(x, v, t, t_settle=0.75 * T_END)
        print("  x0 = %4.1f  ->  settled amplitude = %.4f" % (x0, amp))
    print("  (all three agree: the limit cycle amplitude does not depend on x0)")

    # ---------------- Duffing backbone curve ----------------
    amps = np.linspace(A_MIN, A_MAX, N_AMPLITUDES)
    omega_num = np.array([duffing_frequency_numeric(a) for a in amps])
    omega_lp = duffing_frequency_lindstedt(amps)
    omega_ex = duffing_frequency_exact(amps)

    print("\nFree Duffing oscillator, eps = %.2f" % EPS)
    print("     A      omega_num   omega_exact   omega_Lindstedt")
    for a in (0.25, 0.5, 1.0, 1.5, 2.0):
        print("  %5.2f    %9.5f   %9.5f     %9.5f"
              % (a, duffing_frequency_numeric(a),
                 duffing_frequency_exact(a), duffing_frequency_lindstedt(a)))
    err = np.max(np.abs(omega_num - omega_ex))
    print("  max |numerical - exact| over the whole range: %.2e" % err)

    # ---------------- figure ----------------
    fig = plt.figure(figsize=(11.5, 8.0))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.25)
    plot_vdp_time_history(fig.add_subplot(gs[0, 0]), results)
    plot_vdp_phase_portrait(fig.add_subplot(gs[0, 1]), results)
    plot_backbone(fig.add_subplot(gs[1, :]), amps, omega_num, omega_lp, omega_ex)
    fig.suptitle("1.1 Nonlinear Systems: limit cycle and "
                 "amplitude-dependent natural frequency", fontsize=13)
    plt.show()


if __name__ == "__main__":
    main()
