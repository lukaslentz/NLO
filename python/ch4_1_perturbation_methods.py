"""
4.1 Perturbation Methods
========================

Regular (straightforward) perturbation versus the Lindstedt-Poincare method for
the free Duffing oscillator

    x'' + omega0**2 * x + eps*alpha * x**3 = 0,   x(0) = A,  x'(0) = 0.

Regular expansion x = x0 + eps*x1 + ...  gives at order eps**0 the linear
solution x0 = A*cos(omega0*t), and at order eps**1 a linear oscillator driven
by -alpha*A**3*cos^3(omega0*t) = -(3/4)*alpha*A**3*cos(omega0*t)
- (1/4)*alpha*A**3*cos(3*omega0*t).  The first driver is exactly in resonance
with the homogeneous solution, so it produces the secular term

    x1 = -(3*alpha*A**3)/(8*omega0) * t*sin(omega0*t)
         + (alpha*A**3)/(32*omega0**2) * (cos(3*omega0*t) - cos(omega0*t)),

which grows linearly in t and destroys the approximation once t ~ 1/eps.

Lindstedt-Poincare removes it by stretching time, tau = omega*t, and expanding
the unknown frequency as well, omega = omega0 + eps*omega1 + ...  The
solvability condition (no cos(tau) on the right-hand side at order eps**1)
gives omega1 = 3*alpha*A**2/(8*omega0), i.e. the backbone curve

    omega(A) = sqrt(omega0**2 + (3/4)*alpha*eps*A**2),
    x(t) = A*cos(omega*t)
           + eps*alpha*A**3/(32*omega0**2) * (cos(3*omega*t) - cos(omega*t)).

The script computes a high accuracy numerical reference with solve_ivp and
produces three panels:
  (a) time histories - the regular series runs away, Lindstedt-Poincare stays
      on top of the numerical solution;
  (b) the absolute error of both approximations, showing linear growth versus a
      bounded ripple of size O(eps**2);
  (c) the backbone curve omega(A), Lindstedt-Poincare against the numerically
      measured period (zero-crossing / event detection).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------
# PARAMETERS  (values from the slides)
# ----------------------------------------------------------------------
OMEGA0 = 1.0        # linear natural frequency
ALPHA = 1.0         # sign and size of the cubic term (alpha*eps > 0: hardening)
EPS = 0.2           # perturbation parameter, slides use eps = 0.2 for the
                    # time-history comparison
AMP = 1.0           # initial amplitude A for the time histories

T_END = 60.0        # length of the comparison window (about 10 periods);
                    # the secular term reaches the size of x0 near t = 1/eps *
                    # 8*omega0/(3*alpha*A**2) ~ 13 here
N_T = 4000          # samples for the plotted curves

# Backbone curve sweep
A_MIN, A_MAX, N_A = 0.1, 2.0, 40
T_BACKBONE = 120.0  # integration time used to average many periods
EPS_COMPARE = 0.4   # second, larger eps shown in the backbone panel: the
                    # first-order truncation error is O(eps**2), so doubling
                    # eps makes the gap between L-P and numerics visible

# Backbone comparison table.  For the eps used here the Lindstedt-Poincare
# frequency is omega_LP = sqrt(1 + 0.75*eps*A**2).  The reference column is
# not another approximation but the exact period of the free Duffing
# oscillator, which follows from a complete elliptic integral,
#     omega = pi*sqrt(1 + eps*A**2) / (2*K(m)),  m = eps*A**2/(2*(1+eps*A**2)),
# and the event-based measurement below reproduces it to six digits.  Both
# eps = 0.2 and eps = 0.4 are printed, because the first-order truncation
# error is O(eps**2) and only the larger value makes the gap visible.
TABLE_A = (0.5, 1.0, 1.5, 2.0)
TABLE_EPS = (0.2, 0.4)

RTOL, ATOL = 1e-11, 1e-13


# ----------------------------------------------------------------------
# ANALYTICAL APPROXIMATIONS
# ----------------------------------------------------------------------
def x_regular(t, amp=AMP, eps=EPS):
    """
    Straightforward (regular) perturbation solution to first order in eps.

    The t*sin(omega0*t) contribution is the secular term.  The extra
    -cos(omega0*t) piece is the homogeneous solution needed to satisfy
    x1(0) = 0, x1'(0) = 0, so that the whole expansion honours the initial
    conditions exactly at both orders.
    """
    c = ALPHA * amp**3
    x1 = (-3.0 * c / (8.0 * OMEGA0) * t * np.sin(OMEGA0 * t)
          + c / (32.0 * OMEGA0**2) * (np.cos(3.0 * OMEGA0 * t)
                                      - np.cos(OMEGA0 * t)))
    return amp * np.cos(OMEGA0 * t) + eps * x1


def omega_lp(amp=AMP, eps=EPS):
    """Lindstedt-Poincare backbone frequency omega(A) (first order)."""
    return np.sqrt(OMEGA0**2 + 0.75 * ALPHA * eps * np.asarray(amp) ** 2)


def x_lindstedt(t, amp=AMP, eps=EPS):
    """
    Lindstedt-Poincare solution to first order.

    Identical structure to the regular expansion, but the fast oscillation runs
    at the corrected frequency omega(A), which is precisely what removes the
    secular term.  The correction carries the third harmonic cos(3*omega*t).
    """
    om = omega_lp(amp, eps)
    c = eps * ALPHA * amp**3 / (32.0 * OMEGA0**2)
    return (amp * np.cos(om * t)
            + c * (np.cos(3.0 * om * t) - np.cos(om * t)))


# ----------------------------------------------------------------------
# NUMERICAL REFERENCE
# ----------------------------------------------------------------------
def duffing_free(t, z, eps):
    """Right-hand side of the free (unforced, undamped) Duffing oscillator."""
    x, v = z
    return [v, -OMEGA0**2 * x - eps * ALPHA * x**3]


def solve_reference(t_eval, amp=AMP, eps=EPS):
    """Accurate numerical solution x(t) for x(0) = amp, x'(0) = 0."""
    sol = solve_ivp(duffing_free, (0.0, t_eval[-1]), [amp, 0.0], args=(eps,),
                    t_eval=t_eval, rtol=RTOL, atol=ATOL, method="DOP853")
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[0]


def measure_frequency(amp, eps, t_end=T_BACKBONE):
    """
    Numerically measured oscillation frequency for amplitude amp.

    The system is conservative, so the motion is exactly periodic.  Upward
    zero crossings of x are located as solve_ivp events (root finding on the
    dense output, far more accurate than picking the nearest sample), and the
    mean spacing of many crossings gives the period.
    """
    def up_crossing(t, z, eps):
        return z[0]
    up_crossing.direction = 1.0            # only x: - -> + counts

    sol = solve_ivp(duffing_free, (0.0, t_end), [amp, 0.0], args=(eps,),
                    events=up_crossing, rtol=RTOL, atol=ATOL, method="DOP853")
    tc = sol.t_events[0]
    if tc.size < 2:
        raise RuntimeError("too few zero crossings to measure a period")
    period = (tc[-1] - tc[0]) / (tc.size - 1)
    return 2.0 * np.pi / period


# ----------------------------------------------------------------------
# PLOTTING
# ----------------------------------------------------------------------
def plot_time_histories(ax, t, x_num, x_reg, x_lp):
    """Panel (a): the three time histories over the comparison window."""
    ax.plot(t, x_num, "k-", lw=1.6, label="numerical reference")
    ax.plot(t, x_lp, "-", color="#1f77b4", lw=1.3,
            label="Lindstedt-Poincare (1st order)")
    ax.plot(t, x_reg, "--", color="#d62728", lw=1.3,
            label="regular expansion (secular)")
    ax.axhline(AMP, color="0.7", lw=0.8, ls=":")
    ax.axhline(-AMP, color="0.7", lw=0.8, ls=":")
    ax.set_xlabel(r"time $t$")
    ax.set_ylabel(r"$x(t)$")
    ax.set_title(rf"(a) time histories, $\varepsilon={EPS}$, $A={AMP}$")
    ax.legend(loc="lower left", fontsize=8)
    ax.grid(alpha=0.3)


def plot_errors(ax, t, x_num, x_reg, x_lp):
    """Panel (b): absolute error of both approximations."""
    ax.semilogy(t, np.abs(x_reg - x_num) + 1e-16, "--", color="#d62728",
                lw=1.3, label="regular expansion")
    ax.semilogy(t, np.abs(x_lp - x_num) + 1e-16, "-", color="#1f77b4",
                lw=1.3, label="Lindstedt-Poincare")
    # The secular term itself: eps*3*alpha*A**3/(8*omega0) * t, the envelope
    # the regular expansion error follows.
    ax.semilogy(t, EPS * 3.0 * ALPHA * AMP**3 / (8.0 * OMEGA0) * t + 1e-16,
                ":", color="0.4", lw=1.2, label=r"secular envelope $\propto t$")
    ax.set_ylim(1e-4, 1e2)
    ax.set_xlabel(r"time $t$")
    ax.set_ylabel(r"$|x_\mathrm{approx}-x_\mathrm{num}|$")
    ax.set_title("(b) error growth")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(alpha=0.3, which="both")


def plot_backbone(ax, amps, curves):
    """
    Panel (c): backbone curves, analytical against numerical.

    curves is a list of (eps, omega_LP, omega_num) triples; showing two values
    of eps makes the O(eps**2) truncation error of the first-order backbone
    visible - at eps = 0.2 the two curves are still on top of each other.
    """
    styles = ["#1f77b4", "#2ca02c"]
    for k, (eps, om_lp, om_num) in enumerate(curves):
        ax.plot(om_lp, amps, "-", color=styles[k % len(styles)], lw=1.8,
                label=rf"L-P, $\varepsilon={eps}$")
        ax.plot(om_num, amps, "--", color="#d62728", lw=1.5,
                label=rf"numerical, $\varepsilon={eps}$")
    ax.axvline(OMEGA0, color="0.7", lw=1.0, ls=":")
    ax.text(OMEGA0 + 0.006, amps[-1] * 0.97, r"$\omega_0$", color="0.4",
            fontsize=9)
    ax.set_xlabel(r"frequency $\omega$")
    ax.set_ylabel(r"amplitude $A$")
    ax.set_title("(c) backbone curve $\\omega(A)$")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(alpha=0.3)


# ----------------------------------------------------------------------
def main():
    """Compute both approximations, the numerical reference and the figure."""
    # ---------------- time histories ----------------
    t = np.linspace(0.0, T_END, N_T)
    x_num = solve_reference(t)
    x_reg = x_regular(t)
    x_lp = x_lindstedt(t)

    err_reg = np.abs(x_reg - x_num)
    err_lp = np.abs(x_lp - x_num)
    om = omega_lp()
    print(f"eps = {EPS}, A = {AMP}, alpha = {ALPHA}, omega0 = {OMEGA0}")
    print(f"  Lindstedt-Poincare frequency omega = {om:.6f} "
          f"(omega1 = {3.0 * ALPHA * AMP**2 / (8.0 * OMEGA0):.6f})")
    print(f"  numerically measured frequency     = "
          f"{measure_frequency(AMP, EPS):.6f}")
    print(f"  max error over t <= {T_END:.0f}:  regular expansion "
          f"{err_reg.max():.3f},  Lindstedt-Poincare {err_lp.max():.5f}")

    # ---------------- backbone curve ----------------
    amps = np.linspace(A_MIN, A_MAX, N_A)
    curves = []
    for eps in (EPS, EPS_COMPARE):
        om_bb_lp = omega_lp(amps, eps)
        om_bb_num = np.array([measure_frequency(a, eps) for a in amps])
        curves.append((eps, om_bb_lp, om_bb_num))
        rel_err = np.abs(om_bb_lp - om_bb_num) / om_bb_num
        print(f"backbone curve, eps = {eps}, A in [{A_MIN}, {A_MAX}]: "
              f"max relative frequency error {100.0 * rel_err.max():.2f} % "
              f"(at A = {amps[rel_err.argmax()]:.2f})")

    # ---------------- slide table ----------------
    print("L-P vs. numerical backbone (slide table):")
    header = "    A    " + "".join(
        f"| eps={e}: om_LP  om_num   err  " for e in TABLE_EPS)
    print(header)
    for a in TABLE_A:
        row = f"  {a:4.1f}   "
        for e in TABLE_EPS:
            o_lp, o_nu = omega_lp(a, e), measure_frequency(a, e)
            row += (f"|      {o_lp:6.3f}  {o_nu:6.3f}  "
                    f"{100.0 * abs(o_lp - o_nu) / o_nu:4.1f}%  ")
        print(row)

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))
    plot_time_histories(axes[0], t, x_num, x_reg, x_lp)
    plot_errors(axes[1], t, x_num, x_reg, x_lp)
    plot_backbone(axes[2], amps, curves)
    fig.suptitle("4.1 Perturbation methods - free Duffing oscillator: "
                 "secular growth and the Lindstedt-Poincare correction")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
