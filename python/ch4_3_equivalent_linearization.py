"""
4.3 Equivalent Linearization / Describing Functions
===================================================

Equivalent linearization replaces a nonlinear restoring force f_nl(x) by an
amplitude-dependent linear stiffness.  For an assumed harmonic response
x = A*cos(theta), theta = Omega*t, the coefficient that minimises the
mean-square error <(f_nl - k_eq*x)^2> is the describing function

    N(A) = k_eq(A) = 1/(pi*A) * int_0^{2*pi} f_nl(A*cos(theta))*cos(theta) dth,

i.e. the fundamental Fourier coefficient of the distorted output divided by the
input amplitude.  The nonlinear oscillator

    x'' + 2*zeta*w0*x' + w0^2*x + f_nl(x) = F*cos(Omega*t)

is then replaced by the linear system with stiffness w0^2 + N(A), and the
self-consistency condition A = F/sqrt((w0^2+N(A)-Om^2)^2 + (2*zeta*w0*Om)^2)
gives the amplitude equation -- identical to single-term harmonic balance.

The script computes everything numerically (no hand-derived formulas are used
in the computation, only for comparison) and shows six panels:

  (1) the four nonlinear characteristics treated here: saturation, dead zone,
      cubic stiffness and relay;
  (2) their describing functions N(A) from numerical Fourier integration,
      compared with the closed-form expressions from the lecture;
  (3) frequency response of the Duffing oscillator (f_nl = alpha*x^3,
      hardening) from the describing function -- all three roots of the cubic
      amplitude equation -- against an upward and a downward frequency sweep
      by time integration, which shows the jump phenomenon;
  (4) the same for an oscillator whose *entire* restoring force saturates,
      x'' + 2*zeta*w0*x' + sat(x) = F*cos(Om*t): there N(A) < 1, the resonance
      curve bends to the *left* (softening) and the jump down happens during
      the upward sweep; the describing function used here comes purely from
      the numerical Fourier integral, no closed-form formula is inserted;
  (5) van der Pol oscillator: the equivalent damping c_eq(A) = -mu*(1-A^2/4)
      predicts the limit cycle amplitude A = 2 independent of mu; the panel
      compares this with limit cycle amplitudes from time integration;
  (6) accuracy: relative amplitude error of the describing function against
      time integration at a fixed frequency below resonance, plotted over the
      nonlinearity measure alpha*A^2/w0^2 that the lecture uses to bound the
      error by 5 %.

Runtime: roughly twenty seconds, dominated by the time-integration sweeps
(N_SWEEP and N_PERIODS control this).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ----------------------------------------------------------------------
# Parameters -- change these to explore the model
# ----------------------------------------------------------------------
W0 = 1.0          # linear natural frequency
ZETA = 0.05       # damping ratio
ALPHA = 0.3       # cubic stiffness coefficient of the Duffing oscillator
F_AMP = 0.2       # forcing amplitude

X_SAT = 0.6       # saturation limit of the saturating spring
X_DEAD = 0.5      # half width of the dead zone
M_RELAY = 1.0     # output level of the relay
F_SAT = 0.12      # forcing amplitude for the saturating spring
                  # (must stay below X_SAT, otherwise the spring
                  #  cannot balance the static force at all)

MU_LIST = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]     # van der Pol damping values
ALPHA_LIST = [0.1, 0.3, 0.6]                  # nonlinearity strengths, panel 6
F_ERR_LIST = np.geomspace(0.05, 1.2, 9)       # forcing levels used in panel 6
OM_ERR = 0.9                                  # test frequency for panel 6

N_THETA = 2048    # samples per period for the numerical Fourier integral
N_SWEEP = 26      # frequencies per sweep direction (drives the runtime)
N_PERIODS = 50    # periods discarded as transient in a sweep step
RTOL, ATOL = 1e-9, 1e-11


# ----------------------------------------------------------------------
# Nonlinear characteristics
# ----------------------------------------------------------------------
def f_saturation(x, x_sat=X_SAT):
    """Saturation: linear with unit slope up to x_sat, constant beyond."""
    return np.clip(x, -x_sat, x_sat)


def f_dead_zone(x, x_dead=X_DEAD):
    """Dead zone: no force inside |x| < x_dead, unit slope outside."""
    return np.sign(x) * np.maximum(np.abs(x) - x_dead, 0.0)


def f_cubic(x, alpha=ALPHA):
    """Cubic (Duffing) stiffness."""
    return alpha * x ** 3


def f_relay(x, m=M_RELAY):
    """Ideal relay (sign nonlinearity), discontinuous at x = 0."""
    return m * np.sign(x)


# ----------------------------------------------------------------------
# Describing function by numerical Fourier integration
# ----------------------------------------------------------------------
def describing_function(f_nl, amps):
    """N(A) = (1/(pi*A)) * int_0^{2pi} f_nl(A cos th) cos th dth, numerically.

    The integrand is periodic, so the equidistant trapezoidal rule (here just a
    mean over N_THETA samples) converges spectrally for the smooth
    characteristics.  For the discontinuous relay the convergence drops to
    O(1/N_THETA), which is why a large sample count is used.  Works for a
    single amplitude or for a whole array of amplitudes at once.
    """
    theta = 2.0 * np.pi * np.arange(N_THETA) / N_THETA
    c = np.cos(theta)
    a = np.atleast_1d(amps).astype(float)
    out = 2.0 * np.mean(f_nl(a[:, None] * c[None, :]) * c[None, :], axis=1) / a
    return out if np.ndim(amps) > 0 else float(out[0])


def tabulate_df(f_nl, a_max, n_tab=1500):
    """Tabulate N(A) once and return a fast interpolating function.

    The root finder below evaluates the describing function thousands of times.
    Recomputing the Fourier integral every time would be wasteful, so N(A) is
    sampled on a fine amplitude grid and interpolated linearly afterwards.
    """
    a_tab = np.linspace(1e-4, a_max, n_tab)
    n_val = describing_function(f_nl, a_tab)

    def n_interp(a):
        return np.interp(a, a_tab, n_val)

    return n_interp


# --- closed-form describing functions from the lecture, for comparison ----
def n_saturation_exact(a, x_sat=X_SAT):
    """N(A) of the saturation element (= 1 below the limit)."""
    a = np.atleast_1d(a).astype(float)
    r = np.clip(x_sat / a, 0.0, 1.0)
    n = 2.0 / np.pi * (np.arcsin(r) + r * np.sqrt(1.0 - r ** 2))
    return np.where(a <= x_sat, 1.0, n)


def n_dead_zone_exact(a, x_dead=X_DEAD):
    """N(A) of the dead zone (= 0 below the dead-zone width)."""
    a = np.atleast_1d(a).astype(float)
    r = np.clip(x_dead / a, 0.0, 1.0)
    n = 1.0 - 2.0 / np.pi * (np.arcsin(r) + r * np.sqrt(1.0 - r ** 2))
    return np.where(a <= x_dead, 0.0, n)


def n_cubic_exact(a, alpha=ALPHA):
    """N(A) = 3/4 * alpha * A^2 -- the classical Duffing result."""
    return 0.75 * alpha * np.asarray(a, dtype=float) ** 2


def n_relay_exact(a, m=M_RELAY):
    """N(A) = 4*m/(pi*A) -- the relay gain decays with amplitude."""
    return 4.0 * m / (np.pi * np.asarray(a, dtype=float))


# ----------------------------------------------------------------------
# Amplitude equation of the equivalent linear system
# ----------------------------------------------------------------------
def amplitude_residual(a, om, n_func, force):
    """Self-consistency residual (w0^2 + N(A) - Om^2)^2 + (2 z w0 Om)^2 - (F/A)^2."""
    return ((W0 ** 2 + n_func(a) - om ** 2) ** 2
            + (2.0 * ZETA * W0 * om) ** 2 - (force / a) ** 2)


def df_amplitudes(om, n_func, force, a_min=1e-3, a_max=6.0, n_scan=4000):
    """All positive roots A of the amplitude equation at one frequency.

    The residual is scanned on a fine amplitude grid and every sign change is
    refined with brentq.  Near resonance the equation has three roots (two
    stable, one unstable), elsewhere just one.
    """
    a_grid = np.linspace(a_min, a_max, n_scan)
    r = amplitude_residual(a_grid, om, n_func, force)
    roots = []
    for i in np.where(np.sign(r[:-1]) * np.sign(r[1:]) < 0)[0]:
        roots.append(brentq(amplitude_residual, a_grid[i], a_grid[i + 1],
                            args=(om, n_func, force), xtol=1e-12))
    return np.array(roots)


def df_branches(om_values, n_func, force, a_max=6.0):
    """Amplitude branches over a frequency grid, sorted per frequency.

    Returns a (len(om), 3) array padded with NaN: column 0 is the lowest root,
    column 1 the middle (unstable) one, column 2 the largest.
    """
    out = np.full((len(om_values), 3), np.nan)
    for i, om in enumerate(om_values):
        r = np.sort(df_amplitudes(om, n_func, force, a_max=a_max))
        if len(r) == 1:
            out[i, 0] = r[0]
        elif len(r) >= 3:
            out[i, :3] = r[:3]
        else:                      # 2 roots: numerically at a fold
            out[i, 0], out[i, 2] = r[0], r[-1]
    return out


# ----------------------------------------------------------------------
# Time integration of the full nonlinear system
# ----------------------------------------------------------------------
def make_rhs(f_nl, force):
    """Right-hand side of x'' + 2 z w0 x' + w0^2 x + f_nl(x) = F cos(Om t)."""
    def rhs(t, y, om):
        x, v = y
        return [v, -2.0 * ZETA * W0 * v - W0 ** 2 * x - f_nl(x)
                + force * np.cos(om * t)]
    return rhs


def steady_amplitude(rhs, om, y0, n_periods=N_PERIODS):
    """Fundamental amplitude of the steady state, plus the final state.

    The transient is integrated away first; the last period is projected onto
    cos(Om t) and sin(Om t).  Returning the final state lets the caller use it
    as the initial condition of the next frequency, which is what makes a slow
    sweep follow one branch until it disappears at a fold.
    """
    period = 2.0 * np.pi / om
    t_end = n_periods * period
    t_eval = t_end + np.linspace(0.0, period, 401)
    sol = solve_ivp(rhs, [0.0, t_eval[-1]], y0, args=(om,), t_eval=t_eval,
                    rtol=RTOL, atol=ATOL, method="DOP853")
    t, x = sol.t - t_eval[0], sol.y[0]
    w = np.ones_like(t)
    w[0] = w[-1] = 0.5             # trapezoidal weights over exactly one period
    a = 2.0 / (len(t) - 1) * np.sum(w * x * np.cos(om * t))
    b = 2.0 / (len(t) - 1) * np.sum(w * x * np.sin(om * t))
    return np.hypot(a, b), sol.y[:, -1]


def sweep(f_nl, force, om_values):
    """Slow frequency sweep; each step starts from the previous steady state."""
    rhs = make_rhs(f_nl, force)
    amps, y = [], np.array([0.0, 0.0])
    for om in om_values:
        a, y = steady_amplitude(rhs, om, y)
        amps.append(a)
    return np.array(amps)


# ----------------------------------------------------------------------
# Van der Pol oscillator: equivalent damping and limit cycle
# ----------------------------------------------------------------------
def vdp_limit_cycle_amplitude(mu, t_end=400.0):
    """Amplitude of the van der Pol limit cycle from time integration.

    x'' - mu(1 - x^2) x' + w0^2 x = 0.  After the transient the amplitude is
    read as max|x| over the last few cycles.
    """
    def rhs(t, y):
        x, v = y
        return [v, mu * (1.0 - x ** 2) * v - W0 ** 2 * x]

    sol = solve_ivp(rhs, [0.0, t_end], [0.1, 0.0], rtol=1e-10, atol=1e-12,
                    dense_output=True, method="DOP853")
    t = np.linspace(t_end - 60.0, t_end, 20000)
    return np.max(np.abs(sol.sol(t)[0]))


def vdp_equivalent_damping(a, mu):
    """c_eq(A) = -mu*(1 - A^2/4): fundamental of -mu(1-x^2)x' over A*Omega."""
    return -mu * (1.0 - np.asarray(a, dtype=float) ** 2 / 4.0)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    print("4.3 Equivalent Linearization -- describing functions")
    print(f"    w0 = {W0}, zeta = {ZETA}, alpha = {ALPHA}, F = {F_AMP}\n")

    # --- (1)/(2) characteristics and their describing functions -------
    amps = np.linspace(0.05, 3.0, 200)
    cases = [
        ("saturation", f_saturation, n_saturation_exact, "tab:blue"),
        ("dead zone", f_dead_zone, n_dead_zone_exact, "tab:orange"),
        (r"cubic $\alpha x^3$", f_cubic, n_cubic_exact, "tab:green"),
        ("relay", f_relay, n_relay_exact, "tab:red"),
    ]
    n_num, n_ex = {}, {}
    for name, f_nl, n_exact, _ in cases:
        n_num[name] = describing_function(f_nl, amps)
        n_ex[name] = n_exact(amps)
        err = np.max(np.abs(n_num[name] - n_ex[name]))
        print(f"    N(A) {name:16s}: max |numeric - closed form| = {err:.2e}")

    print(f"    check: N(A=1) cubic = {describing_function(f_cubic, 1.0):.6f}"
          f"  (3/4*alpha = {0.75 * ALPHA:.6f})")

    # --- (3) Duffing frequency response -------------------------------
    om_grid = np.linspace(0.5, 1.8, 400)
    br_duff = df_branches(om_grid, n_cubic_exact, F_AMP)
    om_sweep = np.linspace(0.5, 1.8, N_SWEEP)
    up_duff = sweep(f_cubic, F_AMP, om_sweep)
    down_duff = sweep(f_cubic, F_AMP, om_sweep[::-1])

    i_peak = np.nanargmax(br_duff[:, 2])
    print(f"\n    Duffing: DF peak amplitude {br_duff[i_peak, 2]:.4f} at "
          f"Omega = {om_grid[i_peak]:.4f}; "
          f"time integration peak {np.max(up_duff):.4f} "
          f"at Omega = {om_sweep[np.argmax(up_duff)]:.4f}")

    # --- (4) saturating spring (softening) ----------------------------
    # Here the *whole* restoring force is the saturating characteristic:
    #     x'' + 2 zeta w0 x' + sat(x) = F cos(Om t).
    # The equivalent stiffness is therefore N_sat(A) itself.  The generic
    # routines above always add w0^2, so the difference N_sat(A) - w0^2 is
    # handed over, and for the time integration the extra force
    # sat(x) - w0^2 x is used.  N_sat(A) is tabulated once and interpolated.
    n_sat_tab = tabulate_df(f_saturation, a_max=12.0)

    def n_sat_numeric(a):
        return n_sat_tab(a) - W0 ** 2

    def f_sat_extra(x):
        return f_saturation(x) - W0 ** 2 * x

    # the softening peak sits below w0, so a shifted frequency window is used
    om_grid_sat = np.linspace(0.4, 1.3, 400)
    om_sweep_sat = np.linspace(0.4, 1.3, N_SWEEP)
    br_sat = df_branches(om_grid_sat, n_sat_numeric, F_SAT, a_max=10.0)
    up_sat = sweep(f_sat_extra, F_SAT, om_sweep_sat)
    down_sat = sweep(f_sat_extra, F_SAT, om_sweep_sat[::-1])
    j_peak = np.nanargmax(br_sat[:, 2])
    print(f"    saturating spring: DF peak {br_sat[j_peak, 2]:.4f} at "
          f"Omega = {om_grid_sat[j_peak]:.4f}; "
          f"time integration peak {np.max(down_sat):.4f} "
          f"at Omega = {om_sweep_sat[::-1][np.argmax(down_sat)]:.4f}")

    # --- (5) van der Pol limit cycle ----------------------------------
    a_lc = np.array([vdp_limit_cycle_amplitude(mu) for mu in MU_LIST])
    print("\n    van der Pol limit cycle (DF prediction A = 2):")
    for mu, a in zip(MU_LIST, a_lc):
        print(f"        mu = {mu:4.2f}:  A_numeric = {a:.4f}  "
              f"(error {100 * (a - 2.0) / 2.0:+.1f} %)")

    # --- (6) accuracy of the DF as a function of nonlinearity ---------
    # At a fixed frequency below resonance the amplitude equation has a single
    # root, so the describing-function amplitude and the time-integration
    # amplitude can be compared directly.  Raising the forcing level raises the
    # response amplitude and with it the nonlinearity measure alpha*A^2/w0^2,
    # which is the quantity the lecture uses to bound the error.
    print(f"\n    accuracy test at Omega = {OM_ERR}:")
    nl_measure, rel_error = {}, {}
    for alpha in ALPHA_LIST:
        def n_a(a, alpha=alpha):
            return 0.75 * alpha * np.asarray(a, dtype=float) ** 2

        def f_a(x, alpha=alpha):
            return alpha * x ** 3

        a_df = np.array([df_amplitudes(OM_ERR, n_a, f, a_max=8.0)[0]
                         for f in F_ERR_LIST])
        a_num = []
        for f in F_ERR_LIST:
            rhs_f = make_rhs(f_a, f)
            a_num.append(steady_amplitude(rhs_f, OM_ERR, [0.0, 0.0])[0])
        a_num = np.array(a_num)
        nl_measure[alpha] = alpha * a_num ** 2 / W0 ** 2
        rel_error[alpha] = 100.0 * np.abs(a_df - a_num) / a_num
        print(f"        alpha = {alpha}: amplitude {a_num[0]:.3f} ... "
              f"{a_num[-1]:.3f}, DF error {rel_error[alpha][0]:.2f} % ... "
              f"{rel_error[alpha][-1]:.2f} %")

    # --- plots ---------------------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.8))

    ax = axes[0, 0]
    xx = np.linspace(-2.0, 2.0, 800)
    for name, f_nl, _, col in cases:
        ax.plot(xx, f_nl(xx), color=col, label=name)
    ax.axhline(0.0, color="0.7", lw=0.6)
    ax.axvline(0.0, color="0.7", lw=0.6)
    ax.set_xlabel("$x$")
    ax.set_ylabel(r"$f_{\rm nl}(x)$")
    ax.set_title("Nonlinear characteristics")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    for name, _, _, col in cases:
        ax.plot(amps, n_num[name], color=col, lw=2.0, label=f"{name} (numeric)")
        ax.plot(amps, n_ex[name], "k--", lw=0.9)
    ax.plot([], [], "k--", lw=0.9, label="closed form (lecture)")
    ax.set_ylim(0.0, 2.5)
    ax.set_xlabel("input amplitude $A$")
    ax.set_ylabel("$N(A)$")
    ax.set_title("Describing functions by Fourier integration")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    ax = axes[0, 2]
    ax.plot(om_grid, br_duff[:, 0], "b-", lw=1.5, label="DF, stable branches")
    ax.plot(om_grid, br_duff[:, 2], "b-", lw=1.5)
    ax.plot(om_grid, br_duff[:, 1], "r--", lw=1.5, label="DF, unstable branch")
    a_bb = np.linspace(0.0, 2.6, 100)
    ax.plot(np.sqrt(W0 ** 2 + 0.75 * ALPHA * a_bb ** 2), a_bb, "k:", lw=1.0,
            label="backbone")
    ax.plot(om_sweep, up_duff, "k^", ms=4, mfc="none",
            label=r"time integration, $\Omega\nearrow$")
    ax.plot(om_sweep[::-1], down_duff, "kv", ms=4, mfc="none",
            label=r"time integration, $\Omega\searrow$")
    ax.set_ylim(0.0, 2.6)
    ax.set_xlabel(r"$\Omega$")
    ax.set_ylabel("amplitude $A$")
    ax.set_title(rf"Duffing, hardening ($\alpha={ALPHA}$, $F={F_AMP}$)")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    ax.plot(om_grid_sat, br_sat[:, 0], "b-", lw=1.5,
            label="DF, stable branches")
    ax.plot(om_grid_sat, br_sat[:, 2], "b-", lw=1.5)
    ax.plot(om_grid_sat, br_sat[:, 1], "r--", lw=1.5,
            label="DF, unstable branch")
    ax.axhline(X_SAT, color="0.6", lw=0.8, ls=":")
    ax.text(0.42, X_SAT * 1.05, r"$x_{\rm sat}$", fontsize=7, color="0.4")
    ax.plot(om_sweep_sat, up_sat, "k^", ms=4, mfc="none",
            label=r"time integration, $\Omega\nearrow$")
    ax.plot(om_sweep_sat[::-1], down_sat, "kv", ms=4, mfc="none",
            label=r"time integration, $\Omega\searrow$")
    ax.set_xlabel(r"$\Omega$")
    ax.set_ylabel("amplitude $A$")
    ax.set_title(rf"Saturating spring, softening ($x_{{\rm sat}}={X_SAT}$, "
                 rf"$F={F_SAT}$)")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    a_vdp = np.linspace(0.0, 3.0, 300)
    for mu in [0.1, 0.5, 2.0]:
        ax.plot(a_vdp, vdp_equivalent_damping(a_vdp, mu),
                label=rf"$c_{{\rm eq}}$, $\mu={mu}$")
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.axvline(2.0, color="k", ls="--", lw=1.0, label="DF limit cycle $A=2$")
    ax.plot(a_lc, np.zeros_like(a_lc), "ro", ms=6,
            label="time integration (per $\\mu$)")
    txt = "time integration:\n" + "\n".join(
        rf"$\mu={mu}$:  $A={a:.4f}$" for mu, a in zip(MU_LIST, a_lc))
    ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=7,
            bbox=dict(fc="white", ec="0.8", alpha=0.9))
    ax.set_xlabel("amplitude $A$")
    ax.set_ylabel(r"equivalent damping $c_{\rm eq}(A)$")
    ax.set_title("Van der Pol: limit cycle from $c_{\\rm eq}(A)=0$")
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(alpha=0.3)

    ax = axes[1, 2]
    for alpha in ALPHA_LIST:
        ax.loglog(nl_measure[alpha], rel_error[alpha], "o-", ms=4,
                  label=rf"$\alpha={alpha}$")
    ax.axhline(5.0, color="k", ls=":", lw=1.0, label="5 % error")
    ax.axvline(0.3, color="0.5", ls="--", lw=1.0,
               label=r"$\alpha A^2 = 0.3\,\omega_0^2$")
    ax.set_xlabel(r"nonlinearity measure $\alpha A^2/\omega_0^2$")
    ax.set_ylabel("relative amplitude error [%]")
    ax.set_title(rf"Accuracy of the DF at $\Omega={OM_ERR}$")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3, which="both")

    fig.suptitle("4.3 Equivalent Linearization -- describing functions and "
                 "their frequency responses")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    plt.show()


if __name__ == "__main__":
    main()
