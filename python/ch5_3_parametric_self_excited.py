"""
5.3 Parametric and Self-Excited Oscillations
============================================

This script covers the two mechanisms of chapter 5.3 that produce oscillations
without an additive periodic force.

(1) Parametric excitation.  A system parameter (here the stiffness) is modulated
    periodically.  After normalising the excitation period to 2*pi the equation
    of motion is the damped Mathieu equation

        x'' + 2*D*x' + (delta + eps*cos(t))*x = 0,

    with delta = (2*omega_0/Omega)**2 the squared frequency ratio and eps the
    modulation depth.  Because the equation is linear with 2*pi-periodic
    coefficients, Floquet theory applies: the monodromy matrix M (the map that
    advances the state by one excitation period) decides stability.  Its
    eigenvalues, the Floquet multipliers rho, give the Floquet exponents
    lambda = ln(rho)/(2*pi).  Motion is bounded when every |rho| <= 1.  Scanning
    the (delta, eps) plane produces the Ince-Strutt chart with its instability
    tongues emanating from delta = n**2/4, i.e. from Omega = 2*omega_0/n.
    Viscous damping lifts the tongues off the eps = 0 axis: the principal tongue
    only opens for eps > eps_min = 2*D, its boundaries being
    delta = 1/4 +- 0.5*sqrt(eps**2 - 4*D**2).

(2) Self-excitation.  The Van der Pol oscillator

        x'' - mu*(1 - x**2)*x' + x = 0

    has amplitude-dependent damping d(x) = -mu*(1 - x**2): energy is pumped in
    for |x| < 1 and dissipated for |x| > 1.  The energy balance over one cycle of
    x(t) ~ A*cos(t) gives  mu*pi*A**2*(1 - A**2/4) = 0,  hence the limit cycle
    amplitude A = 2, independent of mu.  The limit cycle is an isolated
    attracting orbit born in a supercritical Hopf bifurcation at mu = 0.

The figure has six panels:
  (a) Ince-Strutt chart of the undamped Mathieu equation, coloured by the Floquet
      growth rate and computed from the monodromy matrix on a (delta, eps) grid;
  (b) the same chart with damping D, showing the closed tongue tips;
  (c) time responses just inside and just outside the principal tongue;
  (d) the Van der Pol damping coefficient d(x) and the energy input per cycle
      W(A), whose zero fixes A = 2;
  (e) the Van der Pol phase plane with trajectories spiralling onto the cycle
      from inside and from outside;
  (f) two time histories that settle on the same amplitude A = 2.

Runtime is roughly 15-25 s, dominated by the two stability charts.
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# =====================================================================
# PARAMETERS  (change these to explore)
# =====================================================================

# --- Ince-Strutt charts ------------------------------------------------
DELTA_MIN, DELTA_MAX = -0.5, 4.5     # squared frequency ratio (2*w0/Omega)**2
EPS_MIN_PLOT, EPS_MAX = 0.0, 2.5     # modulation depth of the stiffness
N_DELTA = 320                        # grid points along delta
N_EPS = 220                          # grid points along eps
N_STEPS = 600                        # RK4 steps per excitation period 2*pi
D_DAMP = 0.10                        # damping ratio for the second chart
# Cost scales as N_DELTA * N_EPS * N_STEPS.  Doubling N_DELTA and N_EPS
# quadruples the runtime of a chart; N_STEPS = 600 already gives a monodromy
# matrix accurate to ~1e-9, so refine the grid rather than the time step.

# --- Time responses inside / outside the principal tongue --------------
EPS_RESP = 0.30                      # modulation depth used for panel (c)
DELTA_IN = 0.25                      # centre of the principal tongue -> unstable
DELTA_OUT = 0.60                     # detuned away from it -> bounded
T_RESP = 80.0                        # integration time (excitation periods*2pi)

# --- Van der Pol -------------------------------------------------------
MU = 1.0                             # nonlinearity / negative damping strength
MU_RELAX = 5.0                       # large mu, relaxation oscillation check
VDP_IC_INSIDE = (0.10, 0.0)          # start close to the unstable focus
VDP_IC_OUTSIDE = (3.50, 0.0)         # start well outside the limit cycle
T_VDP = 60.0                         # total integration time
T_TRANSIENT = 40.0                   # discard this much before measuring

RTOL, ATOL = 1e-10, 1e-12            # tolerances for the solve_ivp runs


# =====================================================================
# Floquet analysis of the Mathieu equation
# =====================================================================

def monodromy_grid(delta, eps, damping=0.0, n_steps=N_STEPS):
    """Monodromy matrix of the damped Mathieu equation on a whole grid.

    The equation x'' + 2*D*x' + (delta + eps*cos(t))*x = 0 is linear, so the
    state after one excitation period T = 2*pi is M @ x(0) with the monodromy
    matrix M.  Its columns are the solutions started from the unit initial
    conditions (1, 0) and (0, 1).  Both columns are integrated here with a
    fixed-step classical RK4 that is vectorised over every (delta, eps) pair of
    the grid at once - orders of magnitude faster than calling an adaptive
    solver once per grid point, and the smooth periodic coefficient makes a
    fixed step perfectly adequate.

    Parameters
    ----------
    delta, eps : ndarray
        Broadcastable arrays of parameter values.
    damping : float
        Damping ratio D in the term 2*D*x'.
    n_steps : int
        RK4 steps per period.

    Returns
    -------
    trace, det : ndarray
        Trace and determinant of M at every grid point.
    """
    delta, eps = np.broadcast_arrays(np.asarray(delta, float),
                                     np.asarray(eps, float))
    # Columns of the fundamental matrix: (x1, v1) from (1,0), (x2, v2) from (0,1)
    x1 = np.ones_like(delta)
    v1 = np.zeros_like(delta)
    x2 = np.zeros_like(delta)
    v2 = np.ones_like(delta)

    def rhs(t, x1, v1, x2, v2):
        """Right-hand side of the first-order form, applied to both columns."""
        k = delta + eps * np.cos(t)          # instantaneous stiffness
        return v1, -k * x1 - 2.0 * damping * v1, \
               v2, -k * x2 - 2.0 * damping * v2

    h = 2.0 * np.pi / n_steps
    t = 0.0
    for _ in range(n_steps):
        a = rhs(t, x1, v1, x2, v2)
        b = rhs(t + 0.5 * h, *(s + 0.5 * h * k for s, k in
                               zip((x1, v1, x2, v2), a)))
        c = rhs(t + 0.5 * h, *(s + 0.5 * h * k for s, k in
                               zip((x1, v1, x2, v2), b)))
        d = rhs(t + h, *(s + h * k for s, k in zip((x1, v1, x2, v2), c)))
        x1, v1, x2, v2 = (s + (h / 6.0) * (ka + 2 * kb + 2 * kc + kd)
                          for s, ka, kb, kc, kd in
                          zip((x1, v1, x2, v2), a, b, c, d))
        t += h

    trace = x1 + v2
    det = x1 * v2 - x2 * v1
    return trace, det


def max_multiplier(trace, det):
    """Largest Floquet multiplier modulus from trace and determinant of M.

    The multipliers solve rho**2 - trace*rho + det = 0.  Working in complex
    arithmetic covers both the oscillatory case (complex pair on a circle of
    radius sqrt(det)) and the unstable case (two real multipliers).
    """
    disc = np.sqrt(np.asarray(trace, complex) ** 2 - 4.0 * det)
    rho1 = 0.5 * (trace + disc)
    rho2 = 0.5 * (trace - disc)
    return np.maximum(np.abs(rho1), np.abs(rho2))


def growth_rate(trace, det):
    """Floquet growth rate Re(lambda) = ln(max|rho|) / T, T = 2*pi."""
    return np.log(max_multiplier(trace, det)) / (2.0 * np.pi)


def tongue_boundary(eps, damping=0.0, side=+1, n_steps=400):
    """Locate the principal-tongue boundary in delta, for verification only.

    The boundary is the delta at which the largest multiplier modulus crosses 1.
    Starting at the tongue centre delta = 1/4 the search marches outwards in
    small steps until the first stable point is met and then bisects.  Marching
    (instead of bracketing with a far-away guess) matters because the plane
    contains further unstable regions - notably all of delta < 0 - which a wide
    bracket would jump across.  Returns NaN if damping has closed the tongue.
    """
    def unstable(d):
        return max_multiplier(*monodromy_grid(d, eps, damping, n_steps)) > 1.0

    inside = 0.25                                   # tongue centre
    if not unstable(inside):
        return np.nan                               # tongue closed by damping
    step = 0.002 * side
    outside = inside
    for _ in range(500):
        outside += step
        if not unstable(outside):
            break
    else:
        return np.nan
    lo, hi = outside - step, outside               # lo unstable, hi stable
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if unstable(mid):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def mathieu_rhs(t, y, delta, eps, damping):
    """First-order form of the damped Mathieu equation for solve_ivp."""
    x, v = y
    return [v, -(delta + eps * np.cos(t)) * x - 2.0 * damping * v]


# =====================================================================
# Van der Pol oscillator
# =====================================================================

def vdp_rhs(t, y, mu):
    """Van der Pol: x'' - mu*(1 - x**2)*x' + x = 0 as a first-order system."""
    x, v = y
    return [v, mu * (1.0 - x * x) * v - x]


def vdp_hopf_rhs(t, y, mu):
    """Unscaled Van der Pol variant x'' - (mu - x**2)*x' + x = 0.

    This is the form in which the Hopf bifurcation at mu = 0 shows the textbook
    amplitude law.  Energy balance with x ~ A*cos(t) gives
    pi*A**2*(mu - A**2/4) = 0, hence A = 2*sqrt(mu):  the limit cycle grows
    continuously out of the equilibrium, A proportional to sqrt(mu).
    The substitution x = sqrt(mu)*z maps it onto the standard Van der Pol
    equation, which is why the standard form always yields A = 2 - there the
    sqrt(mu) has already been absorbed into the amplitude scale.
    """
    x, v = y
    return [v, (mu - x * x) * v - x]


def vdp_energy_per_cycle(amp, mu):
    """Net work of the nonlinear damping over one cycle of x(t) = A*cos(t).

    W(A) = mu*pi*A**2*(1 - A**2/4).  Positive means net energy input (the
    amplitude grows), negative means net dissipation.  The zero at A = 2 is the
    energy balance that fixes the limit cycle amplitude, independently of mu.
    """
    return mu * np.pi * amp ** 2 * (1.0 - amp ** 2 / 4.0)


def measure_limit_cycle(mu, t_end=T_VDP, t_transient=T_TRANSIENT,
                        rhs=vdp_rhs, y0=VDP_IC_INSIDE):
    """Integrate a Van der Pol type oscillator, measure amplitude and period.

    The period is obtained from successive upward zero crossings of x(t) after
    the transient has decayed, refined by linear interpolation between samples.
    """
    sol = solve_ivp(rhs, [0.0, t_end], list(y0), args=(mu,),
                    rtol=RTOL, atol=ATOL, dense_output=True, max_step=0.05)
    t = np.linspace(t_transient, t_end, 40000)
    x = sol.sol(t)[0]
    amplitude = np.max(np.abs(x))

    # upward zero crossings -> one per period
    sign_change = np.where((x[:-1] < 0.0) & (x[1:] >= 0.0))[0]
    crossings = t[sign_change] - x[sign_change] * (t[sign_change + 1] - t[sign_change]) \
        / (x[sign_change + 1] - x[sign_change])
    period = np.mean(np.diff(crossings)) if crossings.size > 1 else np.nan
    return sol, amplitude, period


# =====================================================================
# Main
# =====================================================================

def main():
    """Compute everything, print the numerical checks and draw the figure."""
    # -----------------------------------------------------------------
    # (a)+(b) Ince-Strutt charts
    # -----------------------------------------------------------------
    delta_axis = np.linspace(DELTA_MIN, DELTA_MAX, N_DELTA)
    eps_axis = np.linspace(EPS_MIN_PLOT, EPS_MAX, N_EPS)
    dd, ee = np.meshgrid(delta_axis, eps_axis)

    tr0, det0 = monodromy_grid(dd, ee, damping=0.0)
    tr1, det1 = monodromy_grid(dd, ee, damping=D_DAMP)
    sigma0 = growth_rate(tr0, det0)
    sigma1 = growth_rate(tr1, det1)

    # Numerical health checks: for the undamped equation the Wronskian is
    # conserved, det(M) = 1; with damping det(M) = exp(-2*D*T) exactly.
    det_err0 = np.max(np.abs(det0 - 1.0))
    det_err1 = np.max(np.abs(det1 - np.exp(-2.0 * D_DAMP * 2.0 * np.pi)))
    print("Floquet analysis of the Mathieu equation")
    print("  grid: %d x %d points, %d RK4 steps per period"
          % (N_DELTA, N_EPS, N_STEPS))
    print("  max |det(M) - 1|              (D = 0)      : %.2e" % det_err0)
    print("  max |det(M) - exp(-2*D*T)|    (D = %.2f)   : %.2e"
          % (D_DAMP, det_err1))

    # Verify the perturbation result delta = 1/4 +- eps/2 for small eps.
    print("  principal tongue, undamped (numerical vs. 1/4 -+ eps/2):")
    for eps_test in (0.05, 0.20, 0.50):
        lo = tongue_boundary(eps_test, 0.0, side=-1)
        hi = tongue_boundary(eps_test, 0.0, side=+1)
        print("    eps = %.2f : delta = [%.5f, %.5f]   theory [%.5f, %.5f]"
              % (eps_test, lo, hi, 0.25 - eps_test / 2, 0.25 + eps_test / 2))
    print("  damped tongue (D = %.2f), theory eps_min = 2D = %.2f:"
          % (D_DAMP, 2 * D_DAMP))
    for eps_test in (0.15, 0.25, 0.50):
        lo = tongue_boundary(eps_test, D_DAMP, side=-1)
        hi = tongue_boundary(eps_test, D_DAMP, side=+1)
        if np.isnan(lo):
            print("    eps = %.2f : tongue closed (no instability)" % eps_test)
        else:
            half = 0.5 * np.sqrt(max(eps_test ** 2 - 4 * D_DAMP ** 2, 0.0))
            print("    eps = %.2f : delta = [%.5f, %.5f]   theory [%.5f, %.5f]"
                  % (eps_test, lo, hi, 0.25 - half, 0.25 + half))

    # -----------------------------------------------------------------
    # (c) response inside and outside the principal tongue
    # -----------------------------------------------------------------
    t_eval = np.linspace(0.0, T_RESP, 4000)
    sol_in = solve_ivp(mathieu_rhs, [0.0, T_RESP], [1e-2, 0.0],
                       args=(DELTA_IN, EPS_RESP, 0.0),
                       t_eval=t_eval, rtol=RTOL, atol=ATOL)
    sol_out = solve_ivp(mathieu_rhs, [0.0, T_RESP], [1e-2, 0.0],
                        args=(DELTA_OUT, EPS_RESP, 0.0),
                        t_eval=t_eval, rtol=RTOL, atol=ATOL)
    sigma_in = growth_rate(*monodromy_grid(DELTA_IN, EPS_RESP, 0.0))
    sigma_out = growth_rate(*monodromy_grid(DELTA_OUT, EPS_RESP, 0.0))
    print("  operating points at eps = %.2f:" % EPS_RESP)
    print("    delta = %.2f (inside) : growth rate Re(lambda) = %+.4f"
          " (theory eps/2 = %.4f)" % (DELTA_IN, sigma_in, EPS_RESP / 2))
    print("    delta = %.2f (outside): growth rate Re(lambda) = %+.4f"
          % (DELTA_OUT, sigma_out))

    # -----------------------------------------------------------------
    # (d)-(f) Van der Pol
    # -----------------------------------------------------------------
    sol_vdp_in, amp_in, period_in = measure_limit_cycle(MU)
    sol_vdp_out = solve_ivp(vdp_rhs, [0.0, T_VDP], list(VDP_IC_OUTSIDE),
                            args=(MU,), rtol=RTOL, atol=ATOL,
                            dense_output=True, max_step=0.05)
    _, amp_relax, period_relax = measure_limit_cycle(MU_RELAX, t_end=120.0,
                                                     t_transient=60.0)
    print("Van der Pol oscillator")
    print("  mu = %.2f : limit cycle amplitude = %.4f (energy balance: 2)"
          % (MU, amp_in))
    print("  mu = %.2f : period = %.4f (harmonic limit 2*pi = %.4f)"
          % (MU, period_in, 2 * np.pi))
    # Leading relaxation term (3 - 2*ln2)*mu is only the mu -> infinity limit;
    # the next terms of the asymptotic expansion still matter at mu = 5.
    t_lead = (3.0 - 2.0 * np.log(2.0)) * MU_RELAX
    t_next = t_lead + 7.0143 * MU_RELAX ** (-1.0 / 3.0)
    print("  mu = %.2f : amplitude = %.4f, period = %.4f"
          % (MU_RELAX, amp_relax, period_relax))
    print("             leading term (3-2ln2)*mu = %.4f,"
          " with next term = %.4f" % (t_lead, t_next))

    # The standard Van der Pol equation gives A = 2 for every mu > 0: the
    # sqrt(mu) of the Hopf amplitude law is already absorbed in its scaling.
    # The unscaled variant x'' - (mu - x**2)x' + x = 0 shows the law directly.
    print("  supercritical Hopf, unscaled form  x'' - (mu - x^2)x' + x = 0:")
    for mu_small in (0.05, 0.10, 0.20):
        _, amp_s, _ = measure_limit_cycle(mu_small, t_end=400.0,
                                          t_transient=350.0,
                                          rhs=vdp_hopf_rhs, y0=(0.01, 0.0))
        print("    mu = %.2f : A = %.4f   2*sqrt(mu) = %.4f"
              % (mu_small, amp_s, 2 * np.sqrt(mu_small)))

    # =================================================================
    # Figure
    # =================================================================
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.0))
    fig.suptitle("5.3 Parametric excitation (Mathieu / Ince-Strutt) and "
                 "self-excitation (Van der Pol)", fontsize=14)

    # --- (a) undamped Ince-Strutt ------------------------------------
    ax = axes[0, 0]
    mesh = _plot_chart(ax, delta_axis, eps_axis, sigma0, 0.0)
    ax.set_title("(a) Ince-Strutt chart, undamped")
    fig.colorbar(mesh, ax=ax, label=r"growth rate Re$\lambda$")

    # --- (b) damped Ince-Strutt --------------------------------------
    ax = axes[0, 1]
    mesh = _plot_chart(ax, delta_axis, eps_axis, sigma1, D_DAMP)
    ax.set_title("(b) Ince-Strutt chart, D = %.2f" % D_DAMP)
    fig.colorbar(mesh, ax=ax, label=r"growth rate Re$\lambda$")

    # --- (c) responses -----------------------------------------------
    ax = axes[0, 2]
    ax.plot(sol_in.t, sol_in.y[0], color="crimson", lw=1.0,
            label=r"inside tongue: $\delta$ = %.2f, Re$\lambda$ = %+.3f"
                  % (DELTA_IN, sigma_in))
    ax.plot(sol_out.t, sol_out.y[0], color="tab:blue", lw=1.0,
            label=r"outside: $\delta$ = %.2f, Re$\lambda$ = %+.3f"
                  % (DELTA_OUT, sigma_out))
    # exponential envelope predicted by the Floquet exponent
    env = 1e-2 * np.exp(sigma_in * sol_in.t)
    ax.plot(sol_in.t, env, "k--", lw=1.0, label=r"envelope $e^{\sigma t}$")
    ax.plot(sol_in.t, -env, "k--", lw=1.0)
    # symmetric log scale: the unstable solution grows by five orders of
    # magnitude, the bounded one stays of order the initial amplitude
    ax.set_yscale("symlog", linthresh=1e-2)
    ax.set_xlabel("normalised time $t$ (excitation period $2\\pi$)")
    ax.set_ylabel("$x$ (symlog scale)")
    ax.set_title(r"(c) Response at $\varepsilon$ = %.2f" % EPS_RESP)
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)

    # --- (d) amplitude-dependent damping and energy balance ----------
    ax = axes[1, 0]
    xs = np.linspace(-2.6, 2.6, 400)
    ax.plot(xs, -MU * (1.0 - xs ** 2), color="tab:blue", lw=2,
            label=r"damping $d(x) = -\mu(1-x^2)$")
    amps = np.linspace(0.0, 2.8, 400)
    ax.plot(amps, vdp_energy_per_cycle(amps, MU) / np.pi, color="tab:orange",
            lw=2, label=r"energy per cycle $W(A)/\pi$")
    ax.axhline(0.0, color="k", lw=0.8)
    ax.axvline(2.0, color="crimson", ls=":", lw=1.5, label="$A = 2$ (balance)")
    ax.axvspan(-1.0, 1.0, color="tab:green", alpha=0.10)
    ax.text(0.0, -2.6, "energy input\n$|x|<1$", ha="center", fontsize=8,
            color="tab:green")
    ax.set_xlabel("$x$ resp. amplitude $A$")
    ax.set_ylabel("damping coefficient / energy")
    ax.set_title(r"(d) Amplitude-dependent damping, $\mu$ = %.1f" % MU)
    ax.set_ylim(-3.5, 3.0)
    ax.legend(loc="upper center", fontsize=8)
    ax.grid(alpha=0.3)

    # --- (e) phase plane ---------------------------------------------
    ax = axes[1, 1]
    t_fine = np.linspace(0.0, T_VDP, 20000)
    y_in = sol_vdp_in.sol(t_fine)
    y_out = sol_vdp_out.sol(t_fine)
    ax.plot(y_in[0], y_in[1], color="tab:blue", lw=0.7,
            label="from inside (%.2f, %.1f)" % VDP_IC_INSIDE)
    ax.plot(y_out[0], y_out[1], color="tab:orange", lw=0.7,
            label="from outside (%.2f, %.1f)" % VDP_IC_OUTSIDE)
    mask = t_fine > T_TRANSIENT                    # converged limit cycle
    ax.plot(y_in[0][mask], y_in[1][mask], color="crimson", lw=2.0,
            label="limit cycle")
    ax.plot(0.0, 0.0, "ko", ms=5)
    ax.text(0.15, 0.15, "unstable focus", fontsize=8)
    ax.set_xlabel("$x$")
    ax.set_ylabel(r"$\dot{x}$")
    ax.set_title(r"(e) Van der Pol phase plane, $\mu$ = %.1f" % MU)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.3)

    # --- (f) time histories ------------------------------------------
    ax = axes[1, 2]
    ax.plot(t_fine, y_in[0], color="tab:blue", lw=1.0,
            label="$x_0$ = %.2f" % VDP_IC_INSIDE[0])
    ax.plot(t_fine, y_out[0], color="tab:orange", lw=1.0,
            label="$x_0$ = %.2f" % VDP_IC_OUTSIDE[0])
    ax.axhline(2.0, color="crimson", ls=":", lw=1.5, label="$A = 2$")
    ax.axhline(-2.0, color="crimson", ls=":", lw=1.5)
    ax.set_xlabel("time $t$")
    ax.set_ylabel("$x$")
    ax.set_title("(f) Both initial states settle on $A$ = 2\n"
                 "measured $A$ = %.3f, $T$ = %.3f" % (amp_in, period_in))
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    plt.show()


def _plot_chart(ax, delta_axis, eps_axis, sigma, damping):
    """Draw one Ince-Strutt chart: growth rate map plus analytic boundaries.

    Stable points (Re lambda = 0) are masked and stay white; unstable points are
    shaded by their growth rate.  Returns the mesh so a colorbar can be added.
    """
    unstable = np.ma.masked_where(sigma <= 1e-6, sigma)
    mesh = ax.pcolormesh(delta_axis, eps_axis, unstable, cmap="inferno_r",
                         shading="auto", vmin=0.0, vmax=0.35)
    # first-order transition curves of the principal tongue (n = 1)
    # only drawn up to eps = 1: it is a first-order result in eps
    eps_line = np.linspace(0.0, min(1.0, eps_axis[-1]), 400)
    half = 0.5 * np.sqrt(np.maximum(eps_line ** 2 - 4.0 * damping ** 2, 0.0))
    valid = eps_line >= 2.0 * damping
    ax.plot(0.25 + half[valid], eps_line[valid], "b--", lw=1.2,
            label=r"$\delta = \frac{1}{4} \pm "
                  r"\frac{1}{2}\sqrt{\varepsilon^2-4D^2}$")
    ax.plot(0.25 - half[valid], eps_line[valid], "b--", lw=1.2)
    # tongue origins delta = n**2/4, i.e. Omega = 2*w0/n
    for n in range(1, 5):
        d0 = n ** 2 / 4.0
        if delta_axis[0] <= d0 <= delta_axis[-1]:
            ax.plot(d0, 0.0, "ko", ms=4, clip_on=False)
            ax.annotate("$n$=%d" % n, (d0, 0.0), textcoords="offset points",
                        xytext=(2, 6), fontsize=8)
    ax.set_xlabel(r"$\delta = (2\omega_0/\Omega)^2$")
    ax.set_ylabel(r"modulation depth $\varepsilon$")
    ax.set_xlim(delta_axis[0], delta_axis[-1])
    ax.set_ylim(eps_axis[0], eps_axis[-1])
    ax.legend(loc="upper left", fontsize=8)
    ax.text(0.98, 0.97, "white = stable (bounded)\ncoloured = unstable tongue",
            transform=ax.transAxes, ha="right", va="top", fontsize=8)
    return mesh


if __name__ == "__main__":
    main()
