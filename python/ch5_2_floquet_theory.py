"""
Chapter 5.2 -- Floquet Theory

Monodromy matrix, Floquet multipliers and the Ince-Strutt stability chart of the
damped Mathieu equation.

Governing equation (plain text):

    x'' + 2*zeta*x' + (delta + eps*cos(t))*x = 0,     period T = 2*pi

which is Hill's equation  x'' + p(t) x = 0  with  p(t) = delta + eps*cos(t)
when the damping zeta vanishes.  In first-order form, x = (x, x'),

    x' = A(t) x,   A(t) = [[0, 1], [-(delta + eps*cos t), -2*zeta]],  A(t+T)=A(t).

Floquet theory: the fundamental matrix X(t) with X(0) = I factorises as
X(t) = P(t) exp(B t) with P(t+T) = P(t); the monodromy matrix M = X(T) maps a
perturbation over exactly one period, x(kT) = M^k x(0).  Its eigenvalues are the
Floquet multipliers mu_i, and the trivial solution is asymptotically stable iff
|mu_i| < 1 for all i.  Liouville's formula gives det M = exp(int_0^T tr A dt)
= exp(-4*pi*zeta), so for the undamped Mathieu equation det M = 1, mu_1*mu_2 = 1,
and boundedness is equivalent to |tr M| < 2 (the Hill criterion).

What this script computes
-------------------------
1. M by integrating the augmented system  d/dt vec(X) = A(t) vec(X),  X(0) = I,
   over one period with solve_ivp, and the multipliers mu = eig(M).
   The Liouville identity det M = exp(-4*pi*zeta) is used as a numerical check.
2. The Ince-Strutt chart: the stable region of the (delta, eps) plane, for
   zeta = 0 (Hill criterion |tr M| < 2) and for a damped case (max|mu| < 1).
   The tongue tips are compared with the theoretical values delta = n**2/4.
3. Two scans of the modulation depth eps across the principal tongue (which sits
   at delta = 1/4, i.e. excitation at twice the natural frequency).  At a
   slightly detuned delta = 0.35 the multipliers start as a complex pair on the
   circle |mu| = sqrt(det M) = exp(-2*pi*zeta), collide on the negative real axis
   and one of them leaves the unit circle through mu = -1: the period-doubling
   (flip) bifurcation of the lecture.  At delta = 1/4 the spectral radius
   max|mu| is plotted against eps for three damping levels, which shows the
   minimum modulation depth eps_crit that damping imposes before parametric
   resonance can start.
4. Two time histories, one inside and one outside the principal tongue, with the
   Floquet envelope exp(Re(lambda) t), lambda = ln(mu)/T, drawn on top.

Figure: four panels -- Ince-Strutt chart, multipliers versus the unit circle,
spectral radius versus eps, and the two time histories.

Numerics note
-------------
Single monodromy matrices are computed with scipy's solve_ivp.  The chart needs
about 8*10**4 of them, which is far too slow that way, so it uses a fixed-step
RK4 integrator vectorised over the whole parameter grid (all grid points are
advanced simultaneously with numpy array operations).  Both routines are
compared against each other at a few points; they agree to ~1e-9.

Runtime: about 10 s.  N_DELTA/N_EPS (chart resolution) and RK4_STEPS are the
expensive knobs -- the cost is proportional to N_DELTA*N_EPS*RK4_STEPS.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ----------------------------------------------------------------------
# PARAMETERS  (edit here)
# ----------------------------------------------------------------------
T_PERIOD = 2.0 * np.pi     # period of the Mathieu coefficient p(t)

# Ince-Strutt chart
DELTA_RANGE = (-0.5, 4.5)  # horizontal axis of the chart
EPS_RANGE = (0.0, 2.0)     # vertical axis (modulation depth)
N_DELTA = 251              # grid resolution in delta
N_EPS = 161                # grid resolution in eps
RK4_STEPS = 800            # fixed RK4 steps per period for the grid
ZETA_CHART = 0.05          # damping used for the second (damped) boundary

# Scan across the principal tongue
DELTA_SCAN = 0.25          # tip of the first tongue, delta = n**2/4 with n = 1
EPS_SCAN = np.linspace(0.0, 0.6, 121)
# Path in the complex plane: slightly detuned from the tip, so that the two
# multipliers start as a genuine complex pair and can be seen to collide.
DELTA_PATH = 0.35
EPS_PATH = np.linspace(0.0, 0.45, 91)
ZETA_SCAN = (0.0, 0.02, 0.05)      # undamped and two damping levels

# Time histories (inside / outside the principal tongue)
ZETA_TIME = 0.02
CASE_UNSTABLE = (0.25, 0.40)       # (delta, eps) inside the tongue
CASE_STABLE = (0.70, 0.40)         # same modulation depth, but off resonance
N_PERIODS_TIME = 12                # forcing periods shown

RTOL, ATOL = 1e-11, 1e-13          # tolerances for the solve_ivp monodromy


# ----------------------------------------------------------------------
# Monodromy matrix by direct integration (reference implementation)
# ----------------------------------------------------------------------
def coefficient(t, delta, eps):
    """Periodic stiffness p(t) = delta + eps*cos(t) of the Mathieu equation."""
    return delta + eps * np.cos(t)


def monodromy(delta, eps, zeta=0.0, period=T_PERIOD):
    """Monodromy matrix M = X(T) of the damped Mathieu equation.

    The four entries of the fundamental matrix X are carried as the state
    vector: column j of X is the solution starting from the j-th unit vector, so
    integrating X' = A(t) X from X(0) = I over one period yields M directly.
    The physical state itself never has to be integrated separately, because
    x(t) = X(t) x(0) for every initial condition.
    """
    def rhs(t, y):
        p = coefficient(t, delta, eps)
        x_row, v_row = y[:2], y[2:]          # rows of X: displacements, velocities
        return np.concatenate([v_row, -p * x_row - 2.0 * zeta * v_row])

    y0 = np.array([1.0, 0.0, 0.0, 1.0])      # X(0) = I, stored row-wise
    sol = solve_ivp(rhs, [0.0, period], y0, method="DOP853", rtol=RTOL, atol=ATOL)
    if not sol.success:
        raise RuntimeError("monodromy integration failed: " + sol.message)
    return sol.y[:, -1].reshape(2, 2)


def multipliers(matrix):
    """Floquet multipliers = eigenvalues of the monodromy matrix."""
    return np.linalg.eigvals(matrix)


def floquet_exponents(mu, period=T_PERIOD):
    """Floquet exponents lambda = Ln(mu)/T (principal branch).
    Re(lambda) = ln|mu|/T is unique, Im(lambda) only modulo 2*pi/T."""
    return np.log(mu.astype(complex)) / period


# ----------------------------------------------------------------------
# Monodromy on a whole parameter grid (vectorised fixed-step RK4)
# ----------------------------------------------------------------------
def monodromy_grid(delta, eps, zeta=0.0, period=T_PERIOD, n_step=RK4_STEPS):
    """Trace and determinant of M for every point of a (delta, eps) grid.

    delta and eps are broadcastable arrays.  The state holds both columns of X
    for all grid points at once: y = [x1, v1, x2, v2] with X = [[x1, x2],
    [v1, v2]].  Classical RK4 with a fixed step is used because the number of
    parameter points is large and every point needs exactly the same amount of
    work -- an adaptive per-point solver would cost hundreds of times more.

    Returns (trace M, det M), both with the shape of the broadcast grid.
    """
    delta, eps = np.broadcast_arrays(np.asarray(delta, dtype=float),
                                     np.asarray(eps, dtype=float))
    shape = delta.shape
    y = np.zeros((4,) + shape)
    y[0] = 1.0        # first column of X(0) = e1
    y[3] = 1.0        # second column of X(0) = e2

    def rhs(t, y):
        p = coefficient(t, delta, eps)
        return np.array([y[1], -p * y[0] - 2.0 * zeta * y[1],
                         y[3], -p * y[2] - 2.0 * zeta * y[3]])

    h = period / n_step
    t = 0.0
    for _ in range(n_step):
        k1 = rhs(t, y)
        k2 = rhs(t + 0.5 * h, y + 0.5 * h * k1)
        k3 = rhs(t + 0.5 * h, y + 0.5 * h * k2)
        k4 = rhs(t + h, y + h * k3)
        y = y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        t += h
    trace = y[0] + y[3]
    det = y[0] * y[3] - y[2] * y[1]
    return trace, det


def spectral_radius_from_invariants(trace, det):
    """max|mu| of a real 2x2 matrix from its trace and determinant.

    mu = (tr +- sqrt(tr**2 - 4 det))/2.  For a complex pair (tr**2 < 4 det) both
    multipliers have modulus sqrt(det); otherwise the larger real root decides.
    Working with the invariants avoids assembling 10**5 small matrices.
    """
    disc = trace ** 2 - 4.0 * det
    complex_pair = disc < 0.0
    root = np.sqrt(np.abs(disc))
    real_max = 0.5 * (np.abs(trace) + root)          # larger |real root|
    return np.where(complex_pair, np.sqrt(np.abs(det)), real_max)


# ----------------------------------------------------------------------
# Transition curves of the undamped chart
# ----------------------------------------------------------------------
def tongue_boundaries(eps, n, half_width=0.4, n_scan=4001):
    """Left and right boundary of the n-th instability tongue at modulation eps.

    The transition curves are the level sets |tr M| = 2 of the undamped problem.
    Inside a tongue |tr M| > 2, outside |tr M| < 2, so the boundaries are found
    by walking outwards from the tip delta = n**2/4 until the sign of
    |tr M| - 2 changes, and refining that bracket with brentq.  Higher tongues
    are extremely narrow (width ~ eps**n), hence the fine scan.  Returns
    (delta_left, delta_right); NaN if no crossing is found in the window.
    """
    tip = 0.25 * n ** 2

    def gap(d):
        tr, _ = monodromy_grid(np.asarray(d, dtype=float),
                               np.full(np.shape(d), float(eps)))
        return np.abs(tr) - 2.0

    deltas = np.linspace(tip - half_width, tip + half_width, n_scan)
    values = gap(deltas)
    crossings = np.flatnonzero(np.sign(values[:-1]) != np.sign(values[1:]))
    if crossings.size < 2:
        return (np.nan, np.nan)

    # keep the unstable band (|tr M| > 2) whose midpoint is closest to the tip
    bands = [(crossings[i], crossings[i + 1]) for i in range(len(crossings) - 1)
             if values[crossings[i] + 1] > 0.0]
    if not bands:
        return (np.nan, np.nan)
    i0, i1 = min(bands, key=lambda b: abs(0.5 * (deltas[b[0]] + deltas[b[1]]) - tip))
    left = brentq(lambda d: float(gap(d)), deltas[i0], deltas[i0 + 1], xtol=1e-13)
    right = brentq(lambda d: float(gap(d)), deltas[i1], deltas[i1 + 1], xtol=1e-13)
    return left, right


def critical_eps(delta, zeta, eps_max=2.0):
    """Smallest modulation depth eps at which the orbit loses stability
    (max|mu| = 1) at the given delta and damping."""
    def margin(e):
        tr, det = monodromy_grid(np.array(delta), np.array(e), zeta=zeta)
        return float(spectral_radius_from_invariants(tr, det)) - 1.0

    if margin(eps_max) < 0.0:
        return np.nan
    return brentq(margin, 1e-9, eps_max, xtol=1e-8)


# ----------------------------------------------------------------------
# Time histories
# ----------------------------------------------------------------------
def time_history(delta, eps, zeta, n_periods, x0=(1.0, 0.0), n_out=2000):
    """Solve the Mathieu equation over n_periods for a plot of x(t)."""
    def rhs(t, y):
        return [y[1], -coefficient(t, delta, eps) * y[0] - 2.0 * zeta * y[1]]

    t_end = n_periods * T_PERIOD
    t_eval = np.linspace(0.0, t_end, n_out)
    sol = solve_ivp(rhs, [0.0, t_end], list(x0), t_eval=t_eval,
                    method="DOP853", rtol=1e-10, atol=1e-12)
    return sol.t, sol.y[0]


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    print("Damped Mathieu equation:  x'' + 2*zeta*x' + (delta + eps*cos t)*x = 0,"
          "  T = 2*pi\n")

    # --- one monodromy matrix, printed in full -------------------------
    d0, e0, z0 = 1.0, 0.5, 0.0
    m0 = monodromy(d0, e0, z0)
    mu0 = multipliers(m0)
    print("delta = %.2f, eps = %.2f, zeta = %.2f" % (d0, e0, z0))
    print("  M = [[%+.8f, %+.8f], [%+.8f, %+.8f]]"
          % (m0[0, 0], m0[0, 1], m0[1, 0], m0[1, 1]))
    print("  tr M = %+.8f, det M = %.10f (Liouville: exp(-4*pi*zeta) = %.10f)"
          % (np.trace(m0), np.linalg.det(m0), np.exp(-4.0 * np.pi * z0)))
    print("  multipliers: %s" % np.array2string(mu0, precision=6))
    print("  |mu| = %.6f, %.6f  ->  %s (Hill criterion: |tr M| = %.4f %s 2)"
          % (abs(mu0[0]), abs(mu0[1]),
             "stable" if max(abs(mu0)) < 1.0 + 1e-9 else "unstable",
             abs(np.trace(m0)), "<" if abs(np.trace(m0)) < 2.0 else ">"))

    # damped case: multipliers must shrink by exp(-2*pi*zeta) in modulus
    m1 = monodromy(d0, e0, ZETA_CHART)
    mu1 = multipliers(m1)
    print("\nsame point with zeta = %.2f:" % ZETA_CHART)
    print("  det M = %.10f, exp(-4*pi*zeta) = %.10f"
          % (np.linalg.det(m1), np.exp(-4.0 * np.pi * ZETA_CHART)))
    print("  |mu| = %.6f, %.6f (both = sqrt(det M) = %.6f for a complex pair)"
          % (abs(mu1[0]), abs(mu1[1]), np.sqrt(np.linalg.det(m1))))

    # --- cross-check of the vectorised RK4 against solve_ivp -----------
    check = [(0.25, 1.0), (1.0, 0.5), (2.3, 1.5), (3.9, 1.9)]
    err = 0.0
    for d, e in check:
        tr_rk4, det_rk4 = monodromy_grid(np.array(d), np.array(e))
        m = monodromy(d, e)
        err = max(err, abs(float(tr_rk4) - np.trace(m)),
                  abs(float(det_rk4) - np.linalg.det(m)))
    print("\nvectorised RK4 vs solve_ivp: max deviation of (tr M, det M) = %.2e" % err)

    # --- Ince-Strutt chart ---------------------------------------------
    delta_axis = np.linspace(*DELTA_RANGE, N_DELTA)
    eps_axis = np.linspace(*EPS_RANGE, N_EPS)
    dd, ee = np.meshgrid(delta_axis, eps_axis)

    tr_u, det_u = monodromy_grid(dd, ee, zeta=0.0)
    rho_undamped = spectral_radius_from_invariants(tr_u, det_u)
    tr_d, det_d = monodromy_grid(dd, ee, zeta=ZETA_CHART)
    rho_damped = spectral_radius_from_invariants(tr_d, det_d)

    print("\nInce-Strutt chart on a %d x %d grid" % (N_DELTA, N_EPS))
    eps_probe = 0.2
    for n in (1, 2, 3):
        left, right = tongue_boundaries(eps_probe, n)
        print("  tongue n = %d at eps = %.1f: delta in [%.6f, %.6f]  "
              "(tip theory n**2/4 = %.4f, width %.2e)"
              % (n, eps_probe, left, right, 0.25 * n ** 2, right - left))

    # --- multiplier path across the principal tongue -------------------
    # At delta = 0.35 (just off the tip) the two multipliers start as a complex
    # conjugate pair on the circle |mu| = sqrt(det M) = exp(-2*pi*zeta); raising
    # eps drives them together on the negative real axis, where they split and
    # one of them leaves the unit circle through mu = -1.
    zeta_path = ZETA_SCAN[-1]
    eps_path = EPS_PATH
    mu_path = np.array([multipliers(monodromy(DELTA_PATH, e, zeta_path))
                        for e in eps_path])
    e_flip = critical_eps(DELTA_PATH, zeta_path)
    print("\nmultiplier path at delta = %.2f, zeta = %.2f:" % (DELTA_PATH, zeta_path))
    print("  eps = 0: mu = %s (complex pair, |mu| = %.4f = exp(-2*pi*zeta) = %.4f)"
          % (np.array2string(mu_path[0], precision=4), abs(mu_path[0][0]),
             np.exp(-2.0 * np.pi * zeta_path)))
    print("  unit-circle crossing at eps = %.4f, mu = %s"
          % (e_flip, np.array2string(
              np.sort_complex(multipliers(monodromy(DELTA_PATH, e_flip, zeta_path))),
              precision=4)))

    rho_scan = {}
    for zeta in ZETA_SCAN:
        tr, det = monodromy_grid(DELTA_SCAN, EPS_SCAN, zeta=zeta)
        rho_scan[zeta] = spectral_radius_from_invariants(tr, det)
        e_c = critical_eps(DELTA_SCAN, zeta)
        print("  delta = %.2f, zeta = %.2f: instability threshold eps_crit = %s"
              % (DELTA_SCAN, zeta,
                 "%.4f" % e_c if np.isfinite(e_c) else "none below eps = 2"))

    # multiplier at the flip point: it should pass through mu = -1
    e_c = critical_eps(DELTA_SCAN, zeta_path)
    mu_c = multipliers(monodromy(DELTA_SCAN, e_c, zeta_path))
    print("  at eps_crit (zeta = %.2f) the multipliers are %s  ->  crossing at "
          "mu = -1 (period doubling, response period 2T)"
          % (zeta_path, np.array2string(np.sort_complex(mu_c), precision=4)))

    # --- time histories -------------------------------------------------
    t_un, x_un = time_history(*CASE_UNSTABLE, ZETA_TIME, N_PERIODS_TIME)
    t_st, x_st = time_history(*CASE_STABLE, ZETA_TIME, N_PERIODS_TIME)
    mu_un = multipliers(monodromy(CASE_UNSTABLE[0], CASE_UNSTABLE[1], ZETA_TIME))
    lam_un = floquet_exponents(mu_un)
    growth = np.max(lam_un.real)
    print("\ntime histories at zeta = %.2f:" % ZETA_TIME)
    print("  inside tongue  (delta=%.2f, eps=%.2f): max|mu| = %.4f, "
          "Re(lambda) = %+.4f  -> growth" %
          (CASE_UNSTABLE[0], CASE_UNSTABLE[1], np.max(np.abs(mu_un)), growth))
    mu_st = multipliers(monodromy(CASE_STABLE[0], CASE_STABLE[1], ZETA_TIME))
    print("  outside tongue (delta=%.2f, eps=%.2f): max|mu| = %.4f, "
          "Re(lambda) = %+.4f  -> decay" %
          (CASE_STABLE[0], CASE_STABLE[1], np.max(np.abs(mu_st)),
           np.max(floquet_exponents(mu_st).real)))

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.0), constrained_layout=True)
    ax_chart, ax_circle, ax_scan, ax_time = axes.ravel()

    # (a) Ince-Strutt chart
    ax_chart.contourf(dd, ee, (rho_undamped > 1.0 + 1e-9).astype(float),
                      levels=[0.5, 1.5], colors=["#c6d9f0"])
    ax_chart.contour(dd, ee, rho_undamped, levels=[1.0 + 1e-9],
                     colors="#1f4e79", linewidths=1.6)
    ax_chart.contour(dd, ee, rho_damped, levels=[1.0],
                     colors="#d62728", linewidths=1.6, linestyles="--")
    ax_chart.plot([], [], color="#1f4e79", lw=1.6,
                  label=r"transition curves, $\zeta=0$ ($|\mathrm{tr}\,M|=2$)")
    ax_chart.plot([], [], color="#d62728", lw=1.6, ls="--",
                  label=r"stability boundary, $\zeta=%.2f$ ($\max|\mu|=1$)" % ZETA_CHART)
    ax_chart.fill_between([], [], color="#c6d9f0", label="unstable (undamped)")
    for n in (1, 2, 3, 4):
        tip = 0.25 * n ** 2
        if DELTA_RANGE[0] <= tip <= DELTA_RANGE[1]:
            ax_chart.plot(tip, 0.0, "kv", ms=6)
            ax_chart.annotate(r"$n=%d$" % n, (tip, 0.0), textcoords="offset points",
                              xytext=(3, 8), fontsize=9)
    ax_chart.axvline(0.0, color="0.5", lw=0.8)
    ax_chart.set_xlim(*DELTA_RANGE)
    ax_chart.set_ylim(*EPS_RANGE)
    ax_chart.set_xlabel(r"$\delta$")
    ax_chart.set_ylabel(r"$\varepsilon$")
    ax_chart.set_title(r"Ince--Strutt chart: tongue tips at $\delta=n^2/4$")
    ax_chart.legend(loc="upper right", fontsize=8)

    # (b) multipliers versus the unit circle
    theta = np.linspace(0.0, 2.0 * np.pi, 400)
    ax_circle.plot(np.cos(theta), np.sin(theta), "k-", lw=1.2, label="unit circle")
    ax_circle.plot(np.exp(-2.0 * np.pi * zeta_path) * np.cos(theta),
                   np.exp(-2.0 * np.pi * zeta_path) * np.sin(theta),
                   color="0.6", ls=":", lw=1.2,
                   label=r"$|\mu|=\sqrt{\det M}=e^{-2\pi\zeta}$")
    sc = ax_circle.scatter(mu_path.real.ravel(), mu_path.imag.ravel(),
                           c=np.repeat(eps_path, 2), cmap="viridis", s=18,
                           zorder=3, label=r"multipliers $\mu_{1,2}(\varepsilon)$")
    ax_circle.plot(-1.0, 0.0, "r*", ms=14, zorder=4,
                   label=r"flip point $\mu=-1$")
    fig.colorbar(sc, ax=ax_circle, label=r"$\varepsilon$")
    ax_circle.axhline(0.0, color="0.7", lw=0.8)
    ax_circle.axvline(0.0, color="0.7", lw=0.8)
    ax_circle.set_aspect("equal")
    ax_circle.set_xlim(-2.8, 1.6)
    ax_circle.set_ylim(-1.5, 1.5)
    ax_circle.set_xlabel(r"$\mathrm{Re}\,\mu$")
    ax_circle.set_ylabel(r"$\mathrm{Im}\,\mu$")
    ax_circle.set_title(r"Multipliers at $\delta=%.2f$, $\zeta=%.2f$: exit through $\mu=-1$"
                        % (DELTA_PATH, zeta_path))
    ax_circle.legend(loc="lower left", fontsize=8)

    # (c) spectral radius along the scan
    for zeta, colour in zip(ZETA_SCAN, ["#1f77b4", "#2ca02c", "#d62728"]):
        ax_scan.plot(EPS_SCAN, rho_scan[zeta], color=colour, lw=1.6,
                     label=r"$\zeta=%.2f$" % zeta)
        e_c = critical_eps(DELTA_SCAN, zeta)
        if np.isfinite(e_c):
            ax_scan.plot(e_c, 1.0, "o", color=colour, ms=6)
    ax_scan.axhline(1.0, color="k", ls="--", lw=1.0, label=r"$|\mu|=1$ (stability limit)")
    ax_scan.set_xlabel(r"$\varepsilon$")
    ax_scan.set_ylabel(r"$\max_i|\mu_i|$")
    ax_scan.set_title(r"Crossing the unit circle at $\delta=%.2f$ "
                      r"(principal resonance $\Omega=2\omega_0$)" % DELTA_SCAN)
    ax_scan.set_ylim(0.0, 3.0)
    ax_scan.grid(alpha=0.3)
    ax_scan.legend(loc="upper left", fontsize=9)

    # (d) time histories with the Floquet envelopes.
    # |x| is plotted logarithmically: growth and decay differ by many decades,
    # and an exponential envelope becomes a straight line.
    decay = np.max(floquet_exponents(mu_st).real)
    ax_time.semilogy(t_un, np.abs(x_un), color="#d62728", lw=1.0,
                     label=r"inside tongue: $\delta=%.2f$, $\varepsilon=%.2f$"
                           % CASE_UNSTABLE)
    ax_time.semilogy(t_un, np.abs(x_un[0]) * np.exp(growth * t_un), "k--", lw=1.2,
                     label=r"envelope $e^{\mathrm{Re}(\lambda)t}$, "
                           r"$\mathrm{Re}\,\lambda=%+.3f$" % growth)
    ax_time.semilogy(t_st, np.abs(x_st), color="#1f77b4", lw=1.0,
                     label=r"outside tongue: $\delta=%.2f$, $\varepsilon=%.2f$"
                           % CASE_STABLE)
    ax_time.semilogy(t_st, np.abs(x_st[0]) * np.exp(decay * t_st), "k:", lw=1.2,
                     label=r"envelope $e^{\mathrm{Re}(\lambda)t}$, "
                           r"$\mathrm{Re}\,\lambda=%+.3f$" % decay)
    ax_time.set_ylim(1e-3, 1e6)
    ax_time.set_xlabel(r"$t$")
    ax_time.set_ylabel(r"$|x(t)|$")
    ax_time.set_title(r"Parametric resonance vs. bounded response ($\zeta=%.2f$)"
                      % ZETA_TIME)
    ax_time.grid(alpha=0.3)
    ax_time.legend(loc="upper left", fontsize=9)

    fig.suptitle(r"Floquet analysis of the Mathieu equation "
                 r"$\ddot{x}+2\zeta\dot{x}+(\delta+\varepsilon\cos t)x=0$: "
                 r"monodromy matrix, multipliers and stability chart")
    plt.show()


if __name__ == "__main__":
    main()
