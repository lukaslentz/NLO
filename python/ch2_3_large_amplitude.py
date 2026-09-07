"""
Chapter 2.3 -- Large Amplitude Oscillations (Nonlinear Pendulum)

This script quantifies how the period of a simple pendulum depends on its
swing amplitude, i.e. how the isochronism of the linear oscillator is lost as
soon as sin(theta) may no longer be replaced by theta.

Governing equation (plain text):

    theta'' + omega0^2 * sin(theta) = 0,     omega0 = sqrt(g/L)

Energy first integral (potential normalised to zero at rest):

    E = 0.5 * theta'^2 + omega0^2 * (1 - cos(theta))
    E_separatrix = 2 * omega0^2

For libration with turning-point amplitude theta0 the exact period follows from
the complete elliptic integral of the first kind K(k),

    T = 4/omega0 * K(k),     k = sin(theta0/2),

and the exact time history is given by the Jacobi elliptic sine,

    theta(t) = 2 * arcsin( k * sn(omega0 * t, k) ),
    theta'(t) = 2 * omega0 * k * cn(omega0 * t, k).

Note that SciPy parameterises both ellipk and ellipj by the parameter m = k^2,
not by the modulus k -- this is the single most common source of factor-of-two
confusion when reproducing pendulum tables.

What the script computes
------------------------
1. The exact period T(theta0) from K(k), its two-term series approximation
   T ~ T_lin * (1 + theta0^2/16 + 11*theta0^4/3072), and -- as an independent
   check -- the period measured from a direct numerical integration of the
   nonlinear equation (zero crossings of theta(t)).
2. The frequency backbone omega(theta0)/omega0 = pi / (2 K(k)), which bends to
   the left: the pendulum is the canonical *softening* oscillator.
3. The exact large-amplitude waveform at theta0 = 150 degrees compared with the
   harmonic (small-angle) solution of the same amplitude.

The figure has four panels:
 (a) sin(theta) versus theta with the shaded linearisation error,
 (b) normalised period T/T_lin versus amplitude (exact, series, numerical),
 (c) exact versus harmonic time history at theta0 = 150 degrees,
 (d) phase portrait with libration orbits, the separatrix and a rotating orbit.

A table reproducing the amplitude/period/backbone numbers of the slides is
printed to the console.  Total runtime is a few seconds.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import ellipk, ellipj

# =====================================================================
# PARAMETERS -- change these to experiment
# =====================================================================
G = 9.81                 # gravitational acceleration [m/s^2]
L = 1.0                  # pendulum length [m]

# Amplitudes (in degrees) for the printed table; these are the values of the
# amplitude-period table on the slides.
TABLE_ANGLES_DEG = [10.0, 30.0, 60.0, 90.0, 120.0, 150.0, 170.0]

# Amplitude sweep for the continuous curves (kept away from 180 deg, where the
# period diverges logarithmically).
THETA0_MIN_DEG = 0.5
THETA0_MAX_DEG = 179.0
N_SWEEP = 400

# Amplitude used for the time-history comparison, and how long to plot it
# (measured in linear periods T_lin).
THETA0_HISTORY_DEG = 150.0
N_LINEAR_PERIODS = 3.0

# Numerical integration settings for the independent period measurement
RTOL = 1.0e-11
ATOL = 1.0e-12

OMEGA0 = np.sqrt(G / L)          # linear natural frequency [rad/s]
T_LIN = 2.0 * np.pi / OMEGA0     # linear (small-angle) period [s]


# =====================================================================
# Analytical building blocks
# =====================================================================
def modulus(theta0):
    """Elliptic modulus k = sin(theta0/2) belonging to the amplitude theta0."""
    return np.sin(0.5 * theta0)


def exact_period(theta0):
    """
    Exact libration period T = 4/omega0 * K(k), k = sin(theta0/2).

    SciPy's ellipk takes the parameter m = k^2, hence the square below.
    """
    k = modulus(theta0)
    return 4.0 / OMEGA0 * ellipk(k ** 2)


def series_period(theta0):
    """
    Two-term series approximation of the period,

        T ~ T_lin * (1 + theta0^2/16 + 11/3072 * theta0^4).

    It follows from expanding K(k) in powers of k and k = theta0/2 + ...
    It is accurate only while the amplitude stays moderate (k below ~0.6).
    """
    return T_LIN * (1.0 + theta0 ** 2 / 16.0 + 11.0 * theta0 ** 4 / 3072.0)


def exact_solution(t, theta0):
    """
    Exact libration solution in terms of Jacobi elliptic functions.

    Returns (theta, theta_dot) for the trajectory that passes through
    theta = 0 at t = 0 with the maximum speed belonging to amplitude theta0:

        theta(t)  = 2 arcsin( k sn(omega0 t, k) )
        theta'(t) = 2 omega0 k cn(omega0 t, k)

    The time origin therefore sits at the zero crossing, not at the turning
    point; that is the natural phase for sn, which is odd about u = 0.
    """
    k = modulus(theta0)
    sn, cn, _dn, _ph = ellipj(OMEGA0 * t, k ** 2)   # again: m = k^2
    return 2.0 * np.arcsin(k * sn), 2.0 * OMEGA0 * k * cn


def harmonic_solution(t, theta0):
    """
    Small-angle solution of the *linearised* pendulum with the same amplitude
    and the same phase as exact_solution: theta = theta0 sin(omega0 t).
    """
    return theta0 * np.sin(OMEGA0 * t)


# =====================================================================
# Direct numerical integration (used as an independent check)
# =====================================================================
def pendulum_rhs(t, y):
    """State-space form of theta'' + omega0^2 sin(theta) = 0, y = (theta, theta')."""
    return [y[1], -OMEGA0 ** 2 * np.sin(y[0])]


def numerical_period(theta0):
    """
    Measure the period by integrating the nonlinear equation from rest at
    theta0 and timing successive upward zero crossings of theta(t).

    Half a period elapses between two consecutive zero crossings, so the full
    period is twice the crossing-to-crossing time.  The crossing times are
    located by SciPy's event mechanism, which refines them by root finding and
    is therefore far more accurate than scanning a sampled time series.
    """
    def crossing(t, y):
        return y[0]
    crossing.direction = 0          # count crossings in both directions

    t_guess = exact_period(theta0)  # integrate a bit more than one period
    sol = solve_ivp(pendulum_rhs, (0.0, 1.5 * t_guess), [theta0, 0.0],
                    events=crossing, rtol=RTOL, atol=ATOL, dense_output=False)
    times = sol.t_events[0]
    if len(times) < 2:
        return np.nan
    return 2.0 * (times[1] - times[0])


def integrate_orbit(y0, t_end, n_points=2000):
    """Integrate one pendulum orbit and return (theta, theta') samples."""
    t_eval = np.linspace(0.0, t_end, n_points)
    sol = solve_ivp(pendulum_rhs, (0.0, t_end), y0, t_eval=t_eval,
                    rtol=RTOL, atol=ATOL)
    return sol.y[0], sol.y[1]


# =====================================================================
# Console table
# =====================================================================
def print_table():
    """Print amplitude, modulus, K(k), period ratio and backbone frequency."""
    print(f"Pendulum: g = {G} m/s^2, L = {L} m  ->  omega0 = {OMEGA0:.4f} rad/s, "
          f"T_lin = {T_LIN:.4f} s")
    print()
    print(" theta0      k       K(k)    T_exact   T/T_lin  T_series/T_lin  "
          "T_numeric   omega/omega0   deviation")
    print("-" * 108)
    for deg in TABLE_ANGLES_DEG:
        th0 = np.radians(deg)
        k = modulus(th0)
        kk = ellipk(k ** 2)
        t_ex = exact_period(th0)
        t_se = series_period(th0)
        t_nu = numerical_period(th0)
        ratio = t_ex / T_LIN
        backbone = 1.0 / ratio          # omega/omega0 = T_lin/T
        print(f"{deg:6.0f}  {k:7.4f}  {kk:7.4f}  {t_ex:8.4f}  {ratio:8.4f}  "
              f"{t_se / T_LIN:13.4f}  {t_nu:9.4f}  {backbone:12.4f}  "
              f"{100.0 * (backbone - 1.0):8.1f} %")
    print()

    # Worst-case disagreement between the elliptic-integral formula and the
    # brute-force integration -- a sanity check on both.
    errors = [abs(numerical_period(np.radians(d)) - exact_period(np.radians(d)))
              / exact_period(np.radians(d)) for d in TABLE_ANGLES_DEG]
    print(f"max |T_numeric - T_exact| / T_exact over the table: {max(errors):.2e}")

    # Amplitude at which the linear period is still within 1 %.
    sweep = np.radians(np.linspace(0.1, 90.0, 4000))
    within = sweep[exact_period(sweep) / T_LIN <= 1.01]
    print(f"linear period stays within 1 % up to theta0 = "
          f"{np.degrees(within[-1]):.1f} deg")

    # Separatrix energy and the divergence of the period close to it.
    print(f"separatrix energy E_sep = 2*omega0^2 = {2.0 * OMEGA0 ** 2:.4f}")
    for deg in (179.0, 179.9, 179.99):
        print(f"  theta0 = {deg:7.2f} deg -> T/T_lin = "
              f"{exact_period(np.radians(deg)) / T_LIN:8.4f}")
    print()


# =====================================================================
# Figure
# =====================================================================
def make_figure():
    """Assemble the four-panel summary figure."""
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5))
    ax_sin, ax_period, ax_time, ax_phase = axes.flat

    # ---------------- (a) sin(theta) versus theta -------------------
    th = np.linspace(1.0e-3, np.pi, 500)   # start slightly above 0: the
                                           # relative error below divides by
                                           # sin(theta)
    ax_sin.plot(np.degrees(th), np.sin(th), label=r"exact  $\sin\theta$")
    ax_sin.plot(np.degrees(th), th, "--", label=r"linearised  $\theta$")
    ax_sin.fill_between(np.degrees(th), np.sin(th), th, alpha=0.2,
                        color="tab:red", label="linearisation error")
    # The 5 % error level quoted on the slides.
    rel_err = np.abs(th - np.sin(th)) / np.sin(th)
    idx = np.argmax(rel_err > 0.05)
    ax_sin.axvline(np.degrees(th[idx]), color="0.4", lw=0.8, ls=":")
    ax_sin.annotate(f"5 % error at {np.degrees(th[idx]):.0f}$^\\circ$",
                    xy=(np.degrees(th[idx]), 0.35), xytext=(55, 0.25),
                    fontsize=9, color="0.3")
    ax_sin.set_xlim(0, 180)
    ax_sin.set_ylim(0, 3.2)
    ax_sin.set_xlabel(r"angle $\theta$ [deg]")
    ax_sin.set_ylabel("restoring term")
    ax_sin.set_title("(a) Breakdown of the small-angle approximation")
    ax_sin.legend(loc="upper left", fontsize=9)
    ax_sin.grid(alpha=0.3)

    # ---------------- (b) period versus amplitude -------------------
    th0 = np.radians(np.linspace(THETA0_MIN_DEG, THETA0_MAX_DEG, N_SWEEP))
    ax_period.plot(np.degrees(th0), exact_period(th0) / T_LIN,
                   label=r"exact  $4K(k)/(\omega_0 T_{\rm lin})$")
    ax_period.plot(np.degrees(th0), series_period(th0) / T_LIN, "--",
                   label=r"series  $1+\theta_0^2/16+11\theta_0^4/3072$")
    ax_period.axhline(1.0, color="0.5", lw=0.8, ls=":")
    deg_tab = np.array(TABLE_ANGLES_DEG)
    ax_period.plot(deg_tab,
                   [numerical_period(np.radians(d)) / T_LIN for d in deg_tab],
                   "o", ms=5, mfc="none", color="k",
                   label="numerical integration")
    ax_period.set_xlim(0, 180)
    ax_period.set_ylim(0.95, 3.2)
    ax_period.set_xlabel(r"amplitude $\theta_0$ [deg]")
    ax_period.set_ylabel(r"$T/T_{\rm lin}$")
    ax_period.set_title("(b) Period grows with amplitude (softening)")
    ax_period.legend(loc="center left", fontsize=9)
    ax_period.grid(alpha=0.3)

    # Backbone omega/omega0 on a twin axis: the same information, read as a
    # frequency, is what a forced-response diagram shows as the backbone curve.
    ax_bb = ax_period.twinx()
    ax_bb.plot(np.degrees(th0), T_LIN / exact_period(th0), color="tab:green",
               lw=1.0, alpha=0.8)
    ax_bb.set_ylabel(r"backbone $\omega/\omega_0$", color="tab:green")
    ax_bb.tick_params(axis="y", labelcolor="tab:green")
    ax_bb.set_ylim(0.0, 1.05)

    # ---------------- (c) time history at large amplitude -----------
    th0_h = np.radians(THETA0_HISTORY_DEG)
    t_end = N_LINEAR_PERIODS * T_LIN
    t = np.linspace(0.0, t_end, 3000)
    theta_exact, _ = exact_solution(t, th0_h)
    theta_harm = harmonic_solution(t, th0_h)
    ax_time.plot(t / T_LIN, np.degrees(theta_exact),
                 label=r"exact  $2\arcsin[k\,{\rm sn}(\omega_0 t,k)]$")
    ax_time.plot(t / T_LIN, np.degrees(theta_harm), "--",
                 label=r"harmonic  $\theta_0\sin(\omega_0 t)$")
    ratio = exact_period(th0_h) / T_LIN
    ax_time.set_xlabel(r"time $t/T_{\rm lin}$")
    ax_time.set_ylabel(r"$\theta$ [deg]")
    ax_time.set_title(f"(c) Waveform at "
                      f"$\\theta_0={THETA0_HISTORY_DEG:.0f}^\\circ$ "
                      f"($T={ratio:.2f}\\,T_{{\\rm lin}}$)")
    ax_time.set_ylim(-200.0, 290.0)     # headroom for the legend
    ax_time.legend(loc="upper center", fontsize=9, ncol=2)
    ax_time.grid(alpha=0.3)

    # ---------------- (d) phase portrait ----------------------------
    # Level sets of E = 0.5*thd^2 + omega0^2(1-cos theta); libration inside the
    # separatrix, rotation outside, saddle points at theta = +-pi.
    th_grid = np.linspace(-2.0 * np.pi, 2.0 * np.pi, 800)
    for deg in (30.0, 60.0, 90.0, 120.0, 150.0, 175.0):
        th0_p = np.radians(deg)
        arg = 2.0 * OMEGA0 ** 2 * (np.cos(th_grid) - np.cos(th0_p))
        # Outside the turning points the orbit does not exist: mark those
        # angles as NaN so matplotlib breaks the line instead of connecting
        # neighbouring potential wells with a spurious horizontal segment.
        thd = np.sqrt(np.where(arg > 0.0, arg, np.nan))
        ax_phase.plot(th_grid, thd, color="tab:blue", lw=0.9)
        ax_phase.plot(th_grid, -thd, color="tab:blue", lw=0.9)
    ax_phase.plot([], [], color="tab:blue", lw=0.9, label="libration")

    # Separatrix: thd = +-2 omega0 cos(theta/2)
    sep = 2.0 * OMEGA0 * np.cos(0.5 * th_grid)
    ax_phase.plot(th_grid, sep, "r--", lw=1.3, label="separatrix")
    ax_phase.plot(th_grid, -sep, "r--", lw=1.3)

    # Rotation: total energy above E_sep, integrated numerically.
    for factor in (1.05, 1.25):
        e_tot = factor * 2.0 * OMEGA0 ** 2
        thd0 = np.sqrt(2.0 * e_tot)
        x_rot, v_rot = integrate_orbit([-2.0 * np.pi, thd0],
                                       4.0 * np.pi / (thd0 * 0.5))
        keep = np.abs(x_rot) <= 2.0 * np.pi
        ax_phase.plot(x_rot[keep], v_rot[keep], color="tab:green", lw=1.0)
        ax_phase.plot(x_rot[keep], -v_rot[keep], color="tab:green", lw=1.0)
    ax_phase.plot([], [], color="tab:green", lw=1.0, label="rotation")

    ax_phase.plot([-2 * np.pi, 0, 2 * np.pi], [0, 0, 0], "ko", ms=4)
    ax_phase.plot([-np.pi, np.pi], [0, 0], "rx", ms=7)
    ax_phase.set_xlim(-2.0 * np.pi, 2.0 * np.pi)
    ax_phase.set_ylim(-3.5 * OMEGA0, 3.5 * OMEGA0)
    ax_phase.set_xticks([-2 * np.pi, -np.pi, 0, np.pi, 2 * np.pi])
    ax_phase.set_xticklabels([r"$-2\pi$", r"$-\pi$", "0", r"$\pi$", r"$2\pi$"])
    ax_phase.set_xlabel(r"$\theta$ [rad]")
    ax_phase.set_ylabel(r"$\dot\theta$ [rad/s]")
    ax_phase.set_title("(d) Phase portrait: centres, saddles, separatrix")
    ax_phase.legend(loc="upper right", fontsize=9)
    ax_phase.grid(alpha=0.3)

    fig.suptitle("2.3 Large amplitude oscillations of the pendulum: "
                 r"$\ddot\theta+\omega_0^2\sin\theta=0$")
    fig.tight_layout()


def main():
    print_table()
    make_figure()
    plt.show()


if __name__ == "__main__":
    main()
