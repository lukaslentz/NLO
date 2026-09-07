"""
6.1 Bifurcations
================

A bifurcation is a qualitative change of the solution structure of

    xdot = f(x, mu)

when the parameter mu passes a critical value mu_c.  This script builds the
bifurcation diagrams of the local normal forms of chapter 6.1 numerically, in
exactly the way the chapter describes: by *equilibrium continuation*.  Starting
from one known equilibrium, a tangent predictor steps along the solution branch
and Newton iterations correct back onto f = 0; the stability of each computed
point follows from the Jacobian df/dx (stable if df/dx < 0, unstable if > 0),
and a sign change of df/dx marks a bifurcation point.  The continuation is
formulated in pseudo-arclength form, so it walks around fold points where the
naive parametrisation by mu breaks down (the saddle-node and the subcritical
pitchfork both need this).

Normal forms treated:
    saddle-node (fold)        xdot = mu - x**2
    transcritical             xdot = mu*x - x**2
    pitchfork, supercritical  xdot = mu*x - x**3
    pitchfork, subcritical    xdot = mu*x + x**3 - x**5   (with saturation)
    Hopf, supercritical       rdot = mu*r - r**3,      thetadot = omega
    Hopf, subcritical         rdot = mu*r + r**3 - r**5, thetadot = omega
Finally the period-doubling cascade of the logistic map x -> r*x*(1-x) is
computed as the universal route to chaos, together with a numerical estimate of
the Feigenbaum constant delta_F = 4.669201...

The figure has six panels: the four one-dimensional diagrams (with solid stable
and dashed unstable branches and the bifurcation points marked), the Hopf
amplitude diagram verified against time integration of the planar system, and
the orbit diagram of the logistic map.  Panel (d) additionally shows a
quasi-static forward and backward sweep that traces the hysteresis loop of the
subcritical pitchfork.

Runtime is a few seconds.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib.pyplot as plt

# =====================================================================
# PARAMETERS  (change these to explore)
# =====================================================================

MU_MIN, MU_MAX = -1.5, 2.0      # parameter window of the 1-D diagrams
DS = 0.01                       # pseudo-arclength step of the continuation
MAX_STEPS = 4000                # safety limit per branch
NEWTON_TOL = 1e-12              # corrector tolerance
NEWTON_MAX = 40                 # corrector iterations

# Subcritical pitchfork: xdot = mu*x + x**3 - x**5
# non-trivial equilibria satisfy mu = x**4 - x**2, folds at x**2 = 1/2,
# hence mu_c = -1/4.
SWEEP_MU = np.linspace(-0.6, 0.4, 81)   # quasi-static hysteresis sweep
SWEEP_SETTLE = 300.0                    # relaxation time per sweep point
SWEEP_SEED = 1e-2                       # symmetry-breaking perturbation
# The upward jump is observed slightly above mu = 0: right at the threshold the
# growth rate mu is so small that the state needs longer than SWEEP_SETTLE to
# leave the neighbourhood of x = 0 - the same effect delays the jump in a real
# slow sweep (critical slowing down).

OMEGA_H = 2.0                   # Hopf frequency omega_H
HOPF_MU_CHECK = (0.1, 0.4, 0.9, 1.6)    # mu values verified by integration

# Logistic map / period doubling
R_MIN, R_MAX = 2.8, 4.0         # parameter window of the orbit diagram
N_R = 1200                      # parameter resolution of the orbit diagram
N_TRANSIENT = 600               # iterations discarded per r value
N_PLOT = 200                    # iterations plotted per r value
N_CASCADE = 8                   # number of superstable points used for delta_F
# The orbit diagram costs N_R * (N_TRANSIENT + N_PLOT) map iterations; it is
# vectorised over r, so even N_R = 5000 stays well below a second.


# =====================================================================
# Normal forms:  f(x, mu), df/dx, df/dmu
# =====================================================================

def saddle_node(x, mu):
    """Saddle-node (fold) normal form: xdot = mu - x**2."""
    return mu - x ** 2


def saddle_node_jac(x, mu):
    """df/dx and df/dmu of the saddle-node normal form."""
    return -2.0 * x, 1.0


def transcritical(x, mu):
    """Transcritical normal form: xdot = mu*x - x**2."""
    return mu * x - x ** 2


def transcritical_jac(x, mu):
    """df/dx and df/dmu of the transcritical normal form."""
    return mu - 2.0 * x, x


def pitchfork_super(x, mu):
    """Supercritical pitchfork: xdot = mu*x - x**3 (Z2 symmetric)."""
    return mu * x - x ** 3


def pitchfork_super_jac(x, mu):
    """df/dx and df/dmu of the supercritical pitchfork."""
    return mu - 3.0 * x ** 2, x


def pitchfork_sub(x, mu):
    """Subcritical pitchfork with saturation: xdot = mu*x + x**3 - x**5."""
    return mu * x + x ** 3 - x ** 5


def pitchfork_sub_jac(x, mu):
    """df/dx and df/dmu of the subcritical pitchfork."""
    return mu + 3.0 * x ** 2 - 5.0 * x ** 4, x


def hopf_super(r, mu):
    """Radial equation of the supercritical Hopf normal form: rdot = mu*r-r**3."""
    return mu * r - r ** 3


def hopf_super_jac(r, mu):
    """d/dr and d/dmu of the supercritical Hopf radial equation."""
    return mu - 3.0 * r ** 2, r


def hopf_sub(r, mu):
    """Radial equation of the subcritical Hopf: rdot = mu*r + r**3 - r**5."""
    return mu * r + r ** 3 - r ** 5


def hopf_sub_jac(r, mu):
    """d/dr and d/dmu of the subcritical Hopf radial equation."""
    return mu + 3.0 * r ** 2 - 5.0 * r ** 4, r


# =====================================================================
# Pseudo-arclength continuation
# =====================================================================

def continue_branch(f, jac, x0, mu0, ds=DS, direction=+1,
                    mu_range=(MU_MIN, MU_MAX), x_limit=3.0,
                    max_steps=MAX_STEPS):
    """Trace one equilibrium branch of f(x, mu) = 0 by continuation.

    Algorithm (the one sketched on the 'Bifurcation Diagrams' slide, in its
    fold-safe pseudo-arclength variant):

    1. Tangent predictor.  The tangent (tx, tmu) of the solution curve satisfies
       f_x*tx + f_mu*tmu = 0 and is normalised to unit length; its sign is kept
       consistent with the previous step so the branch is traversed in one
       direction.  The predictor is  (x, mu) + ds*(tx, tmu).
    2. Newton corrector.  Two equations are solved for the two unknowns
       (x, mu):  f(x, mu) = 0 together with the arclength condition
       (x - x_p)*tx + (mu - mu_p)*tmu = ds, which pins the new point onto the
       hyperplane normal to the tangent.  Because mu is an unknown, the method
       turns around fold points instead of diverging there.
    3. Stability.  df/dx is evaluated at every accepted point: negative means
       stable, positive unstable, and a sign change flags a bifurcation point.

    Parameters
    ----------
    f : callable
        f(x, mu).
    jac : callable
        Returns the pair (df/dx, df/dmu).
    x0, mu0 : float
        A point on the branch (need not be exact, it is corrected first).
    ds : float
        Arclength step.
    direction : int
        +1 to start walking towards larger mu, -1 towards smaller mu.  (At a
        fold, where the tangent is vertical, the sign of dx decides instead.)
    mu_range, x_limit : tuple, float
        The continuation stops when it leaves this window.

    Returns
    -------
    mus, xs, stable : ndarray
        Parameter values, equilibrium values and a boolean stability flag.
    """
    x, mu = float(x0), float(mu0)
    # Correct the starting point onto the branch (mu fixed).
    for _ in range(NEWTON_MAX):
        fx, _ = jac(x, mu)
        if abs(fx) < 1e-14:
            break
        dx = -f(x, mu) / fx
        x += dx
        if abs(dx) < NEWTON_TOL:
            break

    xs, mus = [x], [mu]
    fx, fmu = jac(x, mu)
    # initial tangent, normalised
    tx, tmu = -fmu, fx
    norm = np.hypot(tx, tmu)
    tx, tmu = tx / norm, tmu / norm
    # walk towards larger (direction = +1) or smaller (-1) mu
    lead = tmu if abs(tmu) > 1e-8 else tx
    if lead * direction < 0.0:
        tx, tmu = -tx, -tmu

    for _ in range(max_steps):
        # ---- predictor ----
        xp, mup = x + ds * tx, mu + ds * tmu
        xn, mun = xp, mup
        # ---- corrector: Newton on the 2x2 extended system ----
        converged = False
        for _ in range(NEWTON_MAX):
            fx, fmu = jac(xn, mun)
            r1 = f(xn, mun)
            r2 = (xn - x) * tx + (mun - mu) * tmu - ds
            det = fx * tmu - fmu * tx
            if abs(det) < 1e-14:
                break
            dxn = (-r1 * tmu + r2 * fmu) / det
            dmun = (-fx * r2 + tx * r1) / det
            xn += dxn
            mun += dmun
            if abs(dxn) + abs(dmun) < NEWTON_TOL:
                converged = True
                break
        if not converged:
            break

        # ---- new tangent, sign chosen to keep the direction of travel ----
        fx, fmu = jac(xn, mun)
        ntx, ntmu = -fmu, fx
        norm = np.hypot(ntx, ntmu)
        if norm < 1e-14:
            break
        ntx, ntmu = ntx / norm, ntmu / norm
        if ntx * tx + ntmu * tmu < 0.0:
            ntx, ntmu = -ntx, -ntmu
        x, mu, tx, tmu = xn, mun, ntx, ntmu

        xs.append(x)
        mus.append(mu)
        if not (mu_range[0] <= mu <= mu_range[1]) or abs(x) > x_limit:
            break

    xs = np.array(xs)
    mus = np.array(mus)
    stable = np.array([jac(xi, mi)[0] < 0.0 for xi, mi in zip(xs, mus)])
    return mus, xs, stable


def plot_branch(ax, mus, xs, stable, color, label_stable=None,
                label_unstable=None):
    """Plot a continued branch: solid where stable, dashed where unstable.

    The branch is split at every stability change so that the line style really
    follows the sign of df/dx along the curve.
    """
    change = np.flatnonzero(stable[1:] != stable[:-1]) + 1
    segments = np.split(np.arange(xs.size), change)
    used_s, used_u = False, False
    for seg in segments:
        if seg.size < 2:
            continue
        is_stable = stable[seg[0]]
        style = "-" if is_stable else "--"
        lab = None
        if is_stable and not used_s:
            lab, used_s = label_stable, True
        elif not is_stable and not used_u:
            lab, used_u = label_unstable, True
        ax.plot(mus[seg], xs[seg], style, color=color, lw=2, label=lab)


def bifurcation_points(mus, xs, stable):
    """Parameter values at which df/dx changes sign along a branch."""
    change = np.flatnonzero(stable[1:] != stable[:-1])
    return mus[change], xs[change]


# =====================================================================
# Hysteresis sweep of the subcritical pitchfork
# =====================================================================

def quasi_static_sweep(f, mu_values, x_start):
    """Follow the attractor while mu is changed slowly.

    For every mu the previous state is relaxed to its steady state by
    integrating xdot = f(x, mu) long enough.  This mimics an experiment in
    which the parameter is ramped slowly, and it reproduces the jumps and the
    hysteresis loop that the continuation curve only hints at.
    """
    x = x_start
    out = np.empty_like(mu_values)
    for i, mu in enumerate(mu_values):
        sol = solve_ivp(lambda t, y: [f(y[0], mu)], [0.0, SWEEP_SETTLE], [x],
                        rtol=1e-10, atol=1e-12)
        x = float(sol.y[0, -1])
        out[i] = x
    return out


# =====================================================================
# Hopf normal form in the plane
# =====================================================================

def hopf_planar_rhs(t, y, mu, omega, c3, c5):
    """Hopf normal form in Cartesian coordinates.

    In polar form rdot = mu*r + c3*r**3 + c5*r**5 and thetadot = omega, with
    (c3, c5) = (-1, 0) for the supercritical case (first Lyapunov coefficient
    l1 < 0) and (+1, -1) for the subcritical case with saturation.  Written out
    in x, y this is

        xdot = mu*x - omega*y + (c3*R + c5*R**2)*x
        ydot = omega*x + mu*y + (c3*R + c5*R**2)*y ,   R = x**2 + y**2.

    The Jacobian at the origin is [[mu, -omega], [omega, mu]] with eigenvalues
    mu +- i*omega: the pair crosses the imaginary axis at mu = 0, which is the
    Hopf condition.
    """
    x, y = y
    rr = x * x + y * y
    common = mu + c3 * rr + c5 * rr * rr
    return [common * x - omega * y, omega * x + common * y]


def measure_hopf_radius(mu, omega=OMEGA_H, c3=-1.0, c5=0.0, r0=0.3,
                        t_end=200.0):
    """Integrate the planar Hopf system and measure the limit cycle radius."""
    sol = solve_ivp(hopf_planar_rhs, [0.0, t_end], [r0, 0.0],
                    args=(mu, omega, c3, c5), rtol=1e-10, atol=1e-12,
                    dense_output=True)
    t = np.linspace(0.8 * t_end, t_end, 4000)
    x, y = sol.sol(t)
    return float(np.mean(np.hypot(x, y)))


# =====================================================================
# Logistic map: period-doubling cascade
# =====================================================================

def logistic_orbit_diagram(r_values, n_transient=N_TRANSIENT, n_plot=N_PLOT):
    """Attractor of x -> r*x*(1-x) for many r, iterated in parallel.

    Returns flattened (r, x) pairs of the plotted points.  The transient is
    discarded so that only the attractor - fixed point, 2^n cycle or chaotic
    band - remains.
    """
    x = 0.5 * np.ones_like(r_values)
    for _ in range(n_transient):
        x = r_values * x * (1.0 - x)
    r_out = np.repeat(r_values, n_plot)
    x_out = np.empty((r_values.size, n_plot))
    for k in range(n_plot):
        x = r_values * x * (1.0 - x)
        x_out[:, k] = x
    return r_out, x_out.ravel()


def logistic_iterate(r, n, x0=0.5):
    """Apply the logistic map n times to x0."""
    x = x0
    for _ in range(n):
        x = r * x * (1.0 - x)
    return x


def superstable_parameters(n_max=N_CASCADE):
    """Superstable parameter values R_n of the period-doubling cascade.

    The 2**n cycle is superstable when it contains the critical point x = 0.5,
    i.e. when g_n(r) = F^(2**n)(0.5) - 0.5 vanishes.  These points are much
    easier to compute than the bifurcation points themselves and give the same
    Feigenbaum ratio in the limit.  R_n is found by marching upwards from
    R_{n-1} until g_n changes sign and then bracketing the root.
    """
    roots = [2.0]                      # R_0: superstable fixed point at x = 0.5
    gap = 1.0
    for n in range(1, n_max + 1):
        period = 2 ** n

        def g(r, period=period):
            return logistic_iterate(r, period) - 0.5

        r_lo = roots[-1] + 1e-12
        step = max(gap / 25.0, 1e-12)
        g_lo = g(r_lo)
        r_hi = r_lo
        found = False
        for _ in range(4000):
            r_hi += step
            if r_hi > 4.0:
                break
            if g(r_hi) * g_lo < 0.0:
                found = True
                break
            r_lo, g_lo = r_hi, g(r_hi)
        if not found:
            break
        root = brentq(g, r_lo, r_hi, xtol=1e-15, rtol=8.9e-16)
        gap = root - roots[-1]
        roots.append(root)
    return np.array(roots)


def feigenbaum_ratios(roots):
    """delta_n = (R_{n-1} - R_{n-2}) / (R_n - R_{n-1}), converging to 4.6692."""
    gaps = np.diff(roots)
    return gaps[:-1] / gaps[1:]


# =====================================================================
# Main
# =====================================================================

def main():
    """Continue all branches, print the numerical checks and draw the figure."""
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.0))
    fig.suptitle("6.1 Bifurcations: normal forms by equilibrium continuation, "
                 "and the period-doubling route to chaos", fontsize=14)

    # -----------------------------------------------------------------
    # (a) Saddle-node:  xdot = mu - x**2
    #     One single branch x = +-sqrt(mu) that folds at (mu, x) = (0, 0);
    #     continuation started on the upper part and walked around the fold.
    # -----------------------------------------------------------------
    ax = axes[0, 0]
    mus, xs, st = continue_branch(saddle_node, saddle_node_jac,
                                  x0=np.sqrt(MU_MAX), mu0=MU_MAX,
                                  direction=-1)
    plot_branch(ax, mus, xs, st, "tab:blue", "stable", "unstable")
    mu_b, x_b = bifurcation_points(mus, xs, st)
    ax.plot(mu_b, x_b, "ro", ms=7, zorder=5, label="bifurcation point")
    ax.set_title(r"(a) Saddle-node:  $\dot{x} = \mu - x^2$")
    print("Saddle-node   : fold detected at mu = %s (exact 0), x = %s"
          % (np.round(mu_b, 4), np.round(x_b, 4)))
    _finish_1d(ax)

    # -----------------------------------------------------------------
    # (b) Transcritical:  xdot = mu*x - x**2, branches x = 0 and x = mu
    # -----------------------------------------------------------------
    ax = axes[0, 1]
    for x0, mu0, dirn, lab in ((0.0, MU_MIN, +1, ("stable", "unstable")),
                               (MU_MIN, MU_MIN, +1, (None, None))):
        mus, xs, st = continue_branch(transcritical, transcritical_jac,
                                      x0=x0, mu0=mu0, direction=dirn)
        plot_branch(ax, mus, xs, st, "tab:blue", *lab)
    ax.plot(0.0, 0.0, "ro", ms=7, zorder=5, label="bifurcation point")
    ax.set_title(r"(b) Transcritical:  $\dot{x} = \mu x - x^2$")
    print("Transcritical : branches x = 0 and x = mu exchange stability at mu = 0")
    _finish_1d(ax)

    # -----------------------------------------------------------------
    # (c) Supercritical pitchfork:  xdot = mu*x - x**3
    #     Trivial branch plus the two symmetric branches x = +-sqrt(mu).
    # -----------------------------------------------------------------
    ax = axes[0, 2]
    mus, xs, st = continue_branch(pitchfork_super, pitchfork_super_jac,
                                  x0=0.0, mu0=MU_MIN, direction=+1)
    plot_branch(ax, mus, xs, st, "tab:blue", "stable", "unstable")
    for sign in (+1, -1):
        mus2, xs2, st2 = continue_branch(pitchfork_super, pitchfork_super_jac,
                                         x0=sign * np.sqrt(MU_MAX), mu0=MU_MAX,
                                         direction=-1)
        plot_branch(ax, mus2, xs2, st2, "tab:green")
        err = np.max(np.abs(np.abs(xs2) - np.sqrt(np.maximum(mus2, 0.0))))
        print("Pitchfork sup.: branch error vs. sqrt(mu) = %.2e (sign %+d)"
              % (err, sign))
    ax.plot(0.0, 0.0, "ro", ms=7, zorder=5, label="bifurcation point")
    ax.set_title(r"(c) Pitchfork, supercritical:  $\dot{x} = \mu x - x^3$")
    _finish_1d(ax)

    # -----------------------------------------------------------------
    # (d) Subcritical pitchfork:  xdot = mu*x + x**3 - x**5
    #     The non-trivial branches leave x = 0 backwards (unstable), fold at
    #     mu_c = -1/4 and become stable: bistability and hysteresis.
    # -----------------------------------------------------------------
    ax = axes[1, 0]
    mus, xs, st = continue_branch(pitchfork_sub, pitchfork_sub_jac,
                                  x0=0.0, mu0=MU_MIN, direction=+1,
                                  mu_range=(-0.6, 0.4))
    plot_branch(ax, mus, xs, st, "tab:blue", "stable", "unstable")
    fold_mu = []
    for sign in (+1, -1):
        # start on the outer stable part and continue back around the fold
        x_start = sign * np.sqrt(0.5 * (1.0 + np.sqrt(1.0 + 4.0 * 0.4)))
        mus2, xs2, st2 = continue_branch(pitchfork_sub, pitchfork_sub_jac,
                                         x0=x_start, mu0=0.4, direction=-1,
                                         mu_range=(-0.6, 0.4))
        plot_branch(ax, mus2, xs2, st2, "tab:green")
        mb, xb = bifurcation_points(mus2, xs2, st2)
        fold_mu.extend(mb)
        ax.plot(mb, xb, "ro", ms=7, zorder=5)
    ax.plot(0.0, 0.0, "ro", ms=7, zorder=5, label="bifurcation point")
    # Both continuations traverse the complete S-curve (outer stable branch,
    # fold, unstable branch through the origin into the mirrored branch), so
    # each of them detects both folds; the values are collapsed here.
    print("Pitchfork sub.: folds found at mu = %s (exact -0.25)"
          % np.unique(np.round(np.array(fold_mu), 3)))

    # quasi-static sweeps: forward jumps at mu = 0, backward at mu = mu_c
    up = quasi_static_sweep(pitchfork_sub, SWEEP_MU, SWEEP_SEED)
    down = quasi_static_sweep(pitchfork_sub, SWEEP_MU[::-1], up[-1])
    ax.plot(SWEEP_MU, up, ".", color="crimson", ms=4,
            label="sweep $\\mu$ up")
    ax.plot(SWEEP_MU[::-1], down, ".", color="darkorange", ms=4,
            label="sweep $\\mu$ down")
    jump_up = SWEEP_MU[np.argmax(np.abs(np.diff(up, prepend=up[0])) > 0.1)]
    jump_dn = SWEEP_MU[::-1][np.argmax(np.abs(np.diff(down, prepend=down[0]))
                                       > 0.1)]
    print("               hysteresis: jump up at mu = %.3f, "
          "jump back at mu = %.3f" % (jump_up, jump_dn))
    ax.set_title(r"(d) Pitchfork, subcritical:  $\dot{x} = \mu x + x^3 - x^5$")
    _finish_1d(ax, xlim=(-0.6, 0.4))

    # -----------------------------------------------------------------
    # (e) Hopf bifurcation: amplitude of the limit cycle
    # -----------------------------------------------------------------
    ax = axes[1, 1]
    # supercritical: r = sqrt(mu), stable cycle for mu > 0
    mus_h, rs_h, st_h = continue_branch(hopf_super, hopf_super_jac,
                                        x0=np.sqrt(MU_MAX), mu0=MU_MAX,
                                        direction=-1, mu_range=(-0.6, MU_MAX))
    plot_branch(ax, mus_h, rs_h, st_h, "tab:green", "stable cycle",
                "unstable cycle")
    # trivial solution r = 0 (the equilibrium), stable for mu < 0
    mu_line = np.linspace(-0.6, MU_MAX, 200)
    ax.plot(mu_line[mu_line <= 0], np.zeros_like(mu_line[mu_line <= 0]),
            "-", color="tab:blue", lw=2, label="stable equilibrium")
    ax.plot(mu_line[mu_line > 0], np.zeros_like(mu_line[mu_line > 0]),
            "--", color="tab:blue", lw=2, label="unstable equilibrium")
    # subcritical: unstable cycle for mu < 0, folding at mu = -1/4
    mus_s, rs_s, st_s = continue_branch(hopf_sub, hopf_sub_jac,
                                        x0=1.2, mu0=0.4, direction=-1,
                                        mu_range=(-0.6, 0.4))
    plot_branch(ax, mus_s, rs_s, st_s, "tab:red", "subcritical (stable)",
                "subcritical (unstable)")
    # verification by time integration of the planar system
    radii = [measure_hopf_radius(m) for m in HOPF_MU_CHECK]
    ax.plot(HOPF_MU_CHECK, radii, "ko", ms=6, label="time integration")
    print("Hopf          : eigenvalues at the origin are mu +- i*omega,"
          " omega = %.1f" % OMEGA_H)
    for m, r in zip(HOPF_MU_CHECK, radii):
        print("                mu = %+.2f : measured radius %.5f,"
              " sqrt(mu) = %.5f" % (m, r, np.sqrt(m)))
    ax.plot(0.0, 0.0, "ro", ms=7, zorder=5)
    ax.set_title(r"(e) Hopf:  $\dot{r} = \mu r \mp r^3\;(- r^5)$, "
                 r"$\dot{\theta} = \omega_H$")
    ax.set_xlabel(r"parameter $\mu$")
    ax.set_ylabel(r"limit cycle amplitude $r$")
    ax.set_xlim(-0.6, MU_MAX)
    ax.set_ylim(-0.05, 1.6)
    ax.axhline(0.0, color="0.6", lw=0.5)
    ax.axvline(0.0, color="0.6", lw=0.5, ls=":")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=7, ncol=2)

    # -----------------------------------------------------------------
    # (f) Period-doubling cascade of the logistic map
    # -----------------------------------------------------------------
    ax = axes[1, 2]
    r_values = np.linspace(R_MIN, R_MAX, N_R)
    r_pts, x_pts = logistic_orbit_diagram(r_values)
    ax.plot(r_pts, x_pts, ",", color="k", alpha=0.35)
    roots = superstable_parameters()
    ratios = feigenbaum_ratios(roots)
    for r_s in roots[1:]:
        ax.axvline(r_s, color="tab:red", lw=0.6, alpha=0.7)
    ax.set_title("(f) Period doubling, logistic map "
                 r"$x_{n+1} = r\,x_n(1-x_n)$")
    ax.set_xlabel("parameter $r$")
    ax.set_ylabel("attractor $x_n$")
    ax.set_xlim(R_MIN, R_MAX)
    ax.set_ylim(0.0, 1.0)
    ax.text(0.03, 0.05,
            "red: superstable $2^n$ cycles\n"
            r"$\delta_F \approx %.4f$ (exact 4.6692)" % ratios[-1],
            transform=ax.transAxes, fontsize=9,
            bbox=dict(fc="white", alpha=0.8, ec="0.7"))
    print("Period doubling of the logistic map")
    for n, r_s in enumerate(roots):
        print("  superstable R_%d = %.10f" % (n, r_s))
    print("  Feigenbaum ratios : %s" % np.round(ratios, 4))
    print("  best estimate delta_F = %.5f (exact 4.66920)" % ratios[-1])

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    plt.show()


def _finish_1d(ax, xlim=(MU_MIN, MU_MAX)):
    """Common cosmetics of the four one-dimensional bifurcation diagrams."""
    ax.set_xlabel(r"parameter $\mu$")
    ax.set_ylabel(r"equilibrium $x^*$")
    ax.axhline(0.0, color="0.6", lw=0.5)
    ax.axvline(0.0, color="0.6", lw=0.5, ls=":")
    ax.set_xlim(*xlim)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)


if __name__ == "__main__":
    main()
