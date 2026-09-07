"""
Chapter 4.4 -- Forced Duffing Oscillator
========================================

This script demonstrates the two signature phenomena of the harmonically forced
Duffing oscillator with a hardening cubic spring:

1.  The *jump phenomenon and hysteresis*.  A slow (quasi-static) sweep of the
    excitation frequency upwards and then downwards through the primary
    resonance produces two different amplitude curves.  Sweeping up, the
    response climbs the bent-over resonance peak and drops abruptly at the
    jump-down frequency Omega_jd; sweeping down, it stays on the small-amplitude
    branch and jumps up at Omega_ju < Omega_jd.  Between the two jump frequencies
    three steady states coexist (two stable, one unstable) and the observed
    amplitude depends on the loading history.  The simulated sweep is overlaid
    on the frequency-response curve predicted by single-term harmonic balance.

2.  The *route to chaos*, visualised with Poincare sections.  At a fixed
    excitation frequency the forcing amplitude F is raised through a
    period-doubling cascade.  The trajectory is sampled once per excitation
    period T = 2*pi/Omega; a period-n orbit then shows exactly n points, while a
    chaotic attractor produces a fractal cloud of points.  The script verifies
    the character of each response by counting the distinct Poincare points.

Governing equation (plain text):

    x'' + 2*zeta*omega0*x' + omega0**2*x + alpha*x**3 = F*cos(Omega*t)

with omega0 the linear natural frequency, zeta the viscous damping ratio,
alpha > 0 a hardening cubic stiffness, F the forcing amplitude and Omega the
excitation frequency.

The resulting figure has four panels: the frequency sweep with hysteresis (top,
full width) and three Poincare sections for a period-1, a period-2 and a chaotic
response (bottom row).

Runtime: roughly 35-50 s.  The two dominant costs are the two slow sweeps
(SWEEP_DURATION) and the three long Poincare runs (N_PERIODS).  Increasing
either of them raises the runtime proportionally.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

# =====================================================================
# PARAMETERS -- change these to explore the system
# =====================================================================

# --- oscillator ------------------------------------------------------
ZETA = 0.05        # viscous damping ratio
OMEGA0 = 1.0       # linear natural frequency
ALPHA = 1.0        # cubic stiffness, > 0 = hardening spring

# --- part 1: slow frequency sweep (weak forcing, primary resonance) ---
F_SWEEP = 0.3          # forcing amplitude used for the sweep
OMEGA_MIN = 0.5        # start of the up-sweep / end of the down-sweep
OMEGA_MAX = 2.5        # end of the up-sweep / start of the down-sweep
SWEEP_DURATION = 12000.0   # physical time for one sweep; must be long enough
                           # that the sweep is quasi-static (see comment below)
N_SWEEP_SAMPLES = 40000    # output samples per sweep (only for peak picking)

# --- part 2: Poincare sections (strong forcing, chaotic regime) -------
# These values are taken from the lecture slides and are verified below.
OMEGA_POINC = 1.4                     # excitation frequency
F_POINCARE = (15.0, 20.0, 30.0)       # period-1, period-2, chaos
N_PERIODS = 1500                      # total number of excitation periods
N_TRANSIENT = 300                     # periods discarded as transient
Y0_POINCARE = (0.1, 0.0)              # initial condition (x, x')

# --- numerics --------------------------------------------------------
RTOL_SWEEP = 1e-8
ATOL_SWEEP = 1e-10
RTOL_POINC = 1e-9      # a tight tolerance is essential in the chaotic case:
ATOL_POINC = 1e-11     # errors grow exponentially, so the *shape* of the
                       # attractor, not the individual orbit, is the result
CLUSTER_TOL = 2e-3     # two Poincare points closer than this count as one


# =====================================================================
# Model and analytical reference curves
# =====================================================================

def duffing_rhs(t, y, force_amplitude, omega_exc):
    """Right-hand side of the forced Duffing oscillator.

    State y = (x, v).  Returns (x', v') for

        x'  = v
        v'  = -2*zeta*omega0*v - omega0**2*x - alpha*x**3 + F*cos(Omega*t).
    """
    x, v = y
    accel = (-2.0 * ZETA * OMEGA0 * v
             - OMEGA0**2 * x
             - ALPHA * x**3
             + force_amplitude * np.cos(omega_exc * t))
    return [v, accel]


def harmonic_balance_branch(amplitudes, force_amplitude):
    """Frequency-response curve from single-term harmonic balance.

    With the ansatz x(t) ~ A*cos(Omega*t - phi) the slides give the implicit
    amplitude equation

        (omega0**2 + 3*alpha*A**2/4 - Omega**2)**2 + (2*zeta*omega0*Omega)**2
            = (F/A)**2.

    It is *quadratic in u = Omega**2* for a prescribed amplitude A, which is the
    convenient way to trace the folded curve: each A generally yields two
    frequencies (the low-frequency and the high-frequency intersection).

        u**2 + (4*zeta**2*omega0**2 - 2*R)*u + (R**2 - F**2/A**2) = 0,
        R = omega0**2 + 3*alpha*A**2/4.

    Returns two arrays (omega_low, omega_high); entries where no real solution
    exists are set to NaN so that matplotlib simply breaks the line.
    """
    a = np.asarray(amplitudes, dtype=float)
    r = OMEGA0**2 + 0.75 * ALPHA * a**2
    b = 4.0 * ZETA**2 * OMEGA0**2 - 2.0 * r
    c = r**2 - (force_amplitude / a)**2

    disc = b**2 - 4.0 * c
    ok = disc >= 0.0
    sqrt_disc = np.sqrt(np.where(ok, disc, 0.0))

    u_low = 0.5 * (-b - sqrt_disc)
    u_high = 0.5 * (-b + sqrt_disc)

    def to_omega(u):
        u = np.where(ok & (u > 0.0), u, np.nan)
        return np.sqrt(u)

    return to_omega(u_low), to_omega(u_high)


def harmonic_balance_stable(amplitude, omega_exc):
    """Stability of a harmonic-balance solution (A, Omega).

    Averaging / slow-flow analysis of the primary resonance gives the classical
    criterion: the periodic solution is stable as long as

        (omega0**2 + 9*alpha*A**2/4 - Omega**2)
      * (omega0**2 + 3*alpha*A**2/4 - Omega**2)
      + (2*zeta*omega0*Omega)**2  >  0.

    A negative value marks the middle branch between the two saddle-node folds,
    which is never observed in an experiment.
    """
    term1 = OMEGA0**2 + 2.25 * ALPHA * amplitude**2 - omega_exc**2
    term2 = OMEGA0**2 + 0.75 * ALPHA * amplitude**2 - omega_exc**2
    return term1 * term2 + (2.0 * ZETA * OMEGA0 * omega_exc)**2 > 0.0


def backbone_curve(amplitudes):
    """Backbone curve Omega_bb(A) = omega0*sqrt(1 + 3*alpha*A**2/(4*omega0**2)).

    Locus of the undamped, unforced nonlinear natural frequency.  It passes
    through the tips of the resonance peaks as zeta -> 0 and is independent of
    both F and zeta.
    """
    a = np.asarray(amplitudes, dtype=float)
    return OMEGA0 * np.sqrt(1.0 + 0.75 * ALPHA * a**2 / OMEGA0**2)


# =====================================================================
# Part 1: slow frequency sweep
# =====================================================================

def swept_response(omega_start, omega_end):
    """Integrate the Duffing oscillator with a slowly drifting frequency.

    The excitation phase is carried as a third state variable,

        phi' = Omega(t),   Omega(t) = omega_start + rate*t,

    so that the forcing stays phase-continuous while the frequency drifts --
    simply writing cos(Omega(t)*t) would be wrong, because that expression has
    instantaneous frequency Omega + t*dOmega/dt.

    The sweep is quasi-static if the frequency changes by much less than the
    resonance bandwidth (~2*zeta*omega0) during one settling time (~1/(zeta*omega0)).
    Here rate ~ 1.7e-4 while 2*zeta**2*omega0**2 ~ 5e-3, so the sweep is slow
    enough for the response to track the stable branch.

    Returns (omega_at_peak, peak_amplitude): the instantaneous excitation
    frequency and the amplitude at every local maximum of x(t).
    """
    rate = (omega_end - omega_start) / SWEEP_DURATION

    def rhs(t, y):
        x, v, phi = y
        accel = (-2.0 * ZETA * OMEGA0 * v
                 - OMEGA0**2 * x
                 - ALPHA * x**3
                 + F_SWEEP * np.cos(phi))
        return [v, accel, omega_start + rate * t]

    t_eval = np.linspace(0.0, SWEEP_DURATION, N_SWEEP_SAMPLES)
    sol = solve_ivp(rhs, [0.0, SWEEP_DURATION], [0.0, 0.0, 0.0],
                    t_eval=t_eval, method="DOP853",
                    rtol=RTOL_SWEEP, atol=ATOL_SWEEP)

    x = sol.y[0]
    peaks, _ = find_peaks(x)
    omega_inst = omega_start + rate * t_eval[peaks]
    return omega_inst, x[peaks]


def jump_frequency(omega, amplitude, sweeping_up):
    """Locate the jump by the largest step between consecutive response peaks.

    During a sweep the envelope changes smoothly except at the saddle-node,
    where the branch it was following ceases to exist and the response falls
    (or climbs) to the other branch within a few excitation periods.
    """
    d = np.diff(amplitude)
    idx = int(np.argmin(d)) if sweeping_up else int(np.argmax(d))
    return omega[idx]


# =====================================================================
# Part 2: Poincare sections
# =====================================================================

def poincare_section(force_amplitude):
    """Stroboscopic sampling of the steady-state attractor.

    The trajectory is evaluated at t = n*T with T = 2*pi/Omega, i.e. exactly
    once per excitation period, which is the Poincare section
    Sigma = {(x, v) : Omega*t mod 2*pi = 0} of the slides.  The first
    N_TRANSIENT samples are discarded so that only the attractor remains.
    """
    period = 2.0 * np.pi / OMEGA_POINC
    t_end = N_PERIODS * period
    t_strobe = np.arange(N_TRANSIENT, N_PERIODS) * period

    sol = solve_ivp(duffing_rhs, [0.0, t_end], list(Y0_POINCARE),
                    t_eval=t_strobe, args=(force_amplitude, OMEGA_POINC),
                    method="DOP853", rtol=RTOL_POINC, atol=ATOL_POINC)
    return sol.y[0], sol.y[1]


def count_distinct_points(x_values, tol=CLUSTER_TOL):
    """Number of distinct Poincare points, i.e. the period of the orbit.

    Sorting the x-coordinates and counting the gaps larger than `tol` is a
    cheap and robust classifier: 1 -> period-1, 2 -> period-2, a large number
    -> chaotic (or a very high-period orbit).
    """
    x_sorted = np.sort(np.asarray(x_values))
    if x_sorted.size == 0:
        return 0
    return 1 + int(np.count_nonzero(np.diff(x_sorted) > tol))


def classify(n_points):
    """Turn the point count into a human-readable label."""
    if n_points > 20:
        return "chaotic"
    return "period-%d" % n_points


# =====================================================================
# Main
# =====================================================================

def main():
    # ---------------- part 1: hysteresis -----------------------------
    print("Slow frequency sweep, F = %.2f ..." % F_SWEEP)
    om_up, amp_up = swept_response(OMEGA_MIN, OMEGA_MAX)
    om_dn, amp_dn = swept_response(OMEGA_MAX, OMEGA_MIN)

    omega_jd = jump_frequency(om_up, amp_up, sweeping_up=True)
    omega_ju = jump_frequency(om_dn, amp_dn, sweeping_up=False)
    print("  jump-down frequency (up-sweep)   Omega_jd = %.3f" % omega_jd)
    print("  jump-up   frequency (down-sweep) Omega_ju = %.3f" % omega_ju)
    print("  maximum amplitude reached        A_max    = %.3f" % amp_up.max())

    # analytical reference curve, traced by amplitude
    a_grid = np.linspace(0.01, 1.05 * amp_up.max(), 1200)
    om_lo, om_hi = harmonic_balance_branch(a_grid, F_SWEEP)

    # ---------------- part 2: Poincare sections ----------------------
    sections = []
    for force in F_POINCARE:
        print("Poincare section, F = %.1f ..." % force)
        xs, vs = poincare_section(force)
        n_pts = count_distinct_points(xs)
        label = classify(n_pts)
        print("    distinct Poincare points: %4d  ->  %s" % (n_pts, label))
        print("    x in [%.3f, %.3f]" % (xs.min(), xs.max()))
        sections.append((force, xs, vs, label))

    # ---------------- figure -----------------------------------------
    fig = plt.figure(figsize=(11.5, 7.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1.0], hspace=0.48,
                          wspace=0.30, top=0.90, bottom=0.08,
                          left=0.07, right=0.97)

    # -- top panel: sweep + harmonic balance --------------------------
    ax = fig.add_subplot(gs[0, :])

    # split the harmonic-balance curve into stable and unstable parts so that
    # the middle branch can be dashed
    for omega_branch in (om_lo, om_hi):
        stable = harmonic_balance_stable(a_grid, omega_branch)
        ax.plot(np.where(stable, omega_branch, np.nan), a_grid,
                color="0.35", lw=1.6, zorder=2)
        ax.plot(np.where(~stable, omega_branch, np.nan), a_grid,
                color="0.35", lw=1.6, ls="--", zorder=2)
    # proxy handles for the legend
    ax.plot([], [], color="0.35", lw=1.6, label="harmonic balance (stable)")
    ax.plot([], [], color="0.35", lw=1.6, ls="--",
            label="harmonic balance (unstable)")

    ax.plot(backbone_curve(a_grid), a_grid, color="0.6", lw=1.2, ls=":",
            label=r"backbone $\Omega_{bb}(A)$")
    ax.plot(om_up, amp_up, ".", ms=1.6, color="tab:blue",
            label="sweep up")
    ax.plot(om_dn, amp_dn, ".", ms=1.6, color="tab:red",
            label="sweep down")

    ax.axvline(omega_jd, color="tab:blue", lw=0.9, ls="-.", alpha=0.7)
    ax.axvline(omega_ju, color="tab:red", lw=0.9, ls="-.", alpha=0.7)
    ax.annotate(r"$\Omega_{jd}=%.2f$" % omega_jd,
                xy=(omega_jd, 0.95 * amp_up.max()),
                xytext=(6, 0), textcoords="offset points",
                color="tab:blue", fontsize=9)
    ax.annotate(r"$\Omega_{ju}=%.2f$" % omega_ju,
                xy=(omega_ju, 0.95 * amp_up.max()),
                xytext=(-6, 0), textcoords="offset points",
                ha="right", color="tab:red", fontsize=9)

    ax.set_xlim(OMEGA_MIN, OMEGA_MAX)
    ax.set_ylim(0.0, 1.15 * amp_up.max())
    ax.set_xlabel(r"excitation frequency $\Omega$")
    ax.set_ylabel(r"response amplitude $A$")
    ax.set_title(r"Jump phenomenon and hysteresis "
                 r"($\zeta=%.2f$, $\alpha=%.1f$, $F=%.2f$)"
                 % (ZETA, ALPHA, F_SWEEP))
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    ax.grid(alpha=0.25)

    # -- bottom row: Poincare sections --------------------------------
    poincare_axes = []
    for k, (force, xs, vs, label) in enumerate(sections):
        axp = fig.add_subplot(gs[1, k])
        poincare_axes.append(axp)
        marker_size = 1.0 if label == "chaotic" else 8.0
        axp.plot(xs, vs, ".", ms=marker_size, color="tab:purple")
        axp.set_xlabel(r"$x(nT)$")
        if k == 0:
            axp.set_ylabel(r"$\dot{x}(nT)$")
        axp.set_title(r"$F=%.0f$: %s" % (force, label), fontsize=10, pad=6)
        axp.grid(alpha=0.25)

    # a common window for all three sections makes the comparison meaningful:
    # a single point (period-1), two points (period-2), a fractal cloud (chaos)
    for axp in poincare_axes:
        axp.set_xlim(2.4, 5.6)
        axp.set_ylim(-12.0, 12.0)

    fig.suptitle(r"Forced Duffing oscillator: "
                 r"$\ddot{x}+2\zeta\omega_0\dot{x}+\omega_0^2x+\alpha x^3"
                 r"=F\cos\Omega t$", fontsize=12)
    fig.text(0.5, 0.458,
             "Poincare sections at $\\Omega=%.1f$ "
             "(stroboscopic sampling, one point per excitation period)"
             % OMEGA_POINC,
             ha="center", fontsize=10)
    plt.show()


if __name__ == "__main__":
    main()
