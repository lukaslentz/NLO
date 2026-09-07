"""
Chapter 1.2 -- Duffing Oscillator

The Duffing oscillator is the prototypical nonlinear oscillator with a cubic
restoring force.  Governing equation (plain text):

    x'' + 2*zeta*om0*x' + om0**2 * x + alpha * x**3 = f * cos(Omega * t)

with damping ratio zeta, linear natural frequency om0, cubic coefficient alpha
(alpha > 0 hardening, alpha < 0 softening), forcing amplitude f and excitation
frequency Omega.

The script produces three results in one figure:

1. Frequency response.  First-order harmonic balance, x(t) ~ A*cos(Omega*t+phi)
   together with cos^3 = (3/4)cos + (1/4)cos(3.), gives the implicit
   amplitude-frequency relation

       [ (om0^2 + (3/4)*alpha*A^2 - Omega^2)^2 + (2*zeta*om0*Omega)^2 ] * A^2 = f^2

   which is a cubic polynomial in A^2 and is solved exactly for every Omega.
   Up to three real amplitudes coexist; their stability follows from the
   Jacobian of the averaged (slow-flow) equations, so the resonance peak splits
   into two stable branches separated by an unstable middle branch.  On top of
   this the true steady state is obtained by numerical continuation: the
   equation of motion is integrated while Omega is stepped up and then down,
   each step starting from the final state of the previous one.  The two sweeps
   follow different branches and jump at the fold points -- the jump phenomenon
   and its hysteresis loop.  The backbone curve omega^2 = om0^2 + (3/4)*alpha*A^2
   is drawn as well; the response curve leans along it.

2. Potential energy landscape V(x) for a linear, a hardening and a double-well
   configuration.

3. Phase portrait of the unforced, undamped double-well oscillator
   x'' - x + x^3 = 0: two centres at x = +-1, a saddle at the origin and the
   figure-eight separatrix through it.

Only numpy, scipy and matplotlib are used.  Runtime is dominated by the two
numerical frequency sweeps (a few seconds); increasing N_SWEEP or N_PERIODS
makes the script proportionally slower.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =====================================================================
# PARAMETERS -- change these to experiment
# =====================================================================

# --- forced Duffing oscillator (frequency response) -------------------------
ZETA = 0.05        # damping ratio
OM0 = 1.0          # linear natural frequency
ALPHA = 0.5        # cubic coefficient, > 0 -> hardening
F_EXC = 0.3        # forcing amplitude
OM_MIN, OM_MAX = 0.6, 1.8      # excitation frequency range of the sweep

# --- numerical continuation (frequency sweep) -------------------------------
N_SWEEP = 90       # frequency steps per sweep (up and down) -- cost is linear
N_PERIODS = 60     # forcing periods integrated per step ...
N_TRANSIENT = 45   # ... of which the first ones are discarded as transient
                   # (with zeta = 0.05 the transient has decayed by then)
METHOD = "DOP853"  # high-order explicit Runge-Kutta: accurate and fast here
RTOL, ATOL = 1e-8, 1e-10

# --- harmonic balance evaluation grid ---------------------------------------
N_HB = 1200        # frequency points of the analytical response curve

# --- double-well oscillator: x'' - om0d^2 * x + alpha_d * x^3 = 0 ------------
OM0D = 1.0         # magnitude of the (negative) linear stiffness
ALPHA_D = 1.0      # cubic coefficient; equilibria at x = +- om0d/sqrt(alpha_d)


# =====================================================================
# EQUATION OF MOTION
# =====================================================================

def duffing(t, z, om_exc, zeta=ZETA, om0=OM0, alpha=ALPHA, f=F_EXC):
    """Forced Duffing oscillator in state-space form, z = [x, dx/dt].

        dx/dt = v
        dv/dt = -2*zeta*om0*v - om0^2*x - alpha*x^3 + f*cos(om_exc*t)
    """
    x, v = z
    return [v, -2.0 * zeta * om0 * v - om0**2 * x - alpha * x**3
            + f * np.cos(om_exc * t)]


# =====================================================================
# HARMONIC BALANCE: ANALYTICAL FREQUENCY RESPONSE
# =====================================================================

def hb_amplitudes(om_exc, zeta=ZETA, om0=OM0, alpha=ALPHA, f=F_EXC):
    """All real positive harmonic-balance amplitudes at one frequency.

    With u = A^2 the amplitude relation expands into the cubic

        (9/16)*alpha^2 * u^3
      + (3/2)*alpha*(om0^2 - Om^2) * u^2
      + [ (om0^2 - Om^2)^2 + (2*zeta*om0*Om)^2 ] * u
      - f^2 = 0,

    which is solved exactly by np.roots.  Only real, positive roots are
    physical amplitudes.
    """
    d = om0**2 - om_exc**2
    coeffs = [0.5625 * alpha**2,
              1.5 * alpha * d,
              d**2 + (2.0 * zeta * om0 * om_exc)**2,
              -f**2]
    if alpha == 0.0:                      # linear case: the cubic degenerates
        coeffs = coeffs[2:]
    roots = np.roots(coeffs)
    real = roots[np.abs(roots.imag) < 1e-10 * (1.0 + np.abs(roots.real))].real
    return np.sqrt(np.sort(real[real > 0.0]))


def hb_is_stable(amp, om_exc, zeta=ZETA, om0=OM0, alpha=ALPHA):
    """Stability of a harmonic-balance solution from the slow-flow Jacobian.

    Averaging the equation of motion over one forcing period gives amplitude
    and phase equations whose Jacobian at a fixed point has

        trace = -2*zeta*om0  (always negative)
        det   = (zeta*om0)^2 + Lambda_a * Lambda_b,

        Lambda_a = (om0^2 - Om^2 + (3/4)*alpha*A^2) / (2*Om),
        Lambda_b = (om0^2 - Om^2 + (9/4)*alpha*A^2) / (2*Om).

    A negative determinant means a saddle: this is the unstable middle branch
    between the two folds of the bent resonance peak.
    """
    lam_a = (om0**2 - om_exc**2 + 0.75 * alpha * amp**2) / (2.0 * om_exc)
    lam_b = (om0**2 - om_exc**2 + 2.25 * alpha * amp**2) / (2.0 * om_exc)
    return (zeta * om0)**2 + lam_a * lam_b > 0.0


def hb_response_curve(om_grid):
    """Evaluate the harmonic-balance response on a frequency grid.

    Returns four flat arrays: frequency and amplitude of the stable solutions
    and of the unstable ones.
    """
    om_s, a_s, om_u, a_u = [], [], [], []
    for om_exc in om_grid:
        for amp in hb_amplitudes(om_exc):
            if hb_is_stable(amp, om_exc):
                om_s.append(om_exc)
                a_s.append(amp)
            else:
                om_u.append(om_exc)
                a_u.append(amp)
    return (np.array(om_s), np.array(a_s), np.array(om_u), np.array(a_u))


def backbone(amp_grid, om0=OM0, alpha=ALPHA):
    """Backbone curve omega(A) = sqrt(om0^2 + (3/4)*alpha*A^2)."""
    return np.sqrt(om0**2 + 0.75 * alpha * amp_grid**2)


def fold_frequencies(om_grid):
    """Frequency interval in which three harmonic-balance solutions coexist.

    The fold (turning) points of the response curve are bracketed by counting
    the number of real solutions along the frequency grid.
    """
    counts = np.array([len(hb_amplitudes(om)) for om in om_grid])
    multi = om_grid[counts >= 3]
    if multi.size == 0:
        return None, None
    return multi.min(), multi.max()


# =====================================================================
# NUMERICAL CONTINUATION (FREQUENCY SWEEP)
# =====================================================================

def sweep(om_grid):
    """Steady-state amplitude along a slow frequency sweep.

    For every excitation frequency the equation of motion is integrated for
    N_PERIODS forcing periods, starting from the final state of the previous
    frequency (that is what makes it a continuation rather than a set of
    independent runs).  The amplitude is read off the last N_PERIODS -
    N_TRANSIENT periods, after the transient has decayed.  Because the sweep
    carries its state along, it stays on one branch until that branch
    disappears at a fold -- and then jumps to the other one.
    """
    amps = np.empty_like(om_grid)
    z0 = [0.0, 0.0]
    for i, om_exc in enumerate(om_grid):
        period = 2.0 * np.pi / om_exc
        t_end = N_PERIODS * period
        t_eval = np.linspace(N_TRANSIENT * period, t_end, 300)
        sol = solve_ivp(duffing, (0.0, t_end), z0, args=(om_exc,),
                        t_eval=t_eval, rtol=RTOL, atol=ATOL, method=METHOD)
        amps[i] = np.max(np.abs(sol.y[0]))
        z0 = sol.y[:, -1]        # warm start for the next frequency
    return amps


# =====================================================================
# POTENTIAL AND DOUBLE-WELL PHASE PORTRAIT
# =====================================================================

def potential(x, om0_sq, alpha):
    """Potential energy V(x) = om0_sq/2 * x^2 + alpha/4 * x^4.

    A double well requires a negative linear stiffness (om0_sq < 0) together
    with alpha > 0; this is the buckled-beam configuration
    x'' + 2*zeta*om0*x' - om0^2*x + alpha*x^3 = f*cos(Om*t) of the slides.  A
    negative alpha alone would make V unbounded from below rather than
    bistable.
    """
    return 0.5 * om0_sq * x**2 + 0.25 * alpha * x**4


def double_well_energy(x, v, om0d=OM0D, alpha_d=ALPHA_D):
    """Total energy of x'' - om0d^2*x + alpha_d*x^3 = 0.

    E = v^2/2 - om0d^2*x^2/2 + alpha_d*x^4/4 is conserved, so the orbits of the
    phase portrait are exactly the level sets of E.  E = 0 through the saddle
    at the origin is the figure-eight separatrix.
    """
    return 0.5 * v**2 + potential(x, -om0d**2, alpha_d)


# =====================================================================
# PLOTTING
# =====================================================================

def plot_frequency_response(ax, om_hb, hb, om_sweep, up, down):
    """Bent resonance peak with stable/unstable branches and both sweeps."""
    om_s, a_s, om_u, a_u = hb

    ax.plot(om_s, a_s, ".", ms=2.5, color="tab:blue",
            label="harmonic balance, stable")
    ax.plot(om_u, a_u, ".", ms=2.5, color="tab:red",
            label="harmonic balance, unstable")

    a_bb = np.linspace(0.0, 1.05 * a_s.max(), 200)
    om_bb = backbone(a_bb)
    ax.plot(om_bb[om_bb <= OM_MAX], a_bb[om_bb <= OM_MAX], "k:", lw=1.4,
            label=r"backbone $\omega^2=\omega_0^2+\frac{3}{4}\alpha A^2$")

    ax.plot(om_sweep, up, "-", lw=1.3, color="tab:green",
            label="sweep up (jumps down at right fold)")
    ax.plot(om_sweep, down, "--", lw=1.3, color="tab:orange",
            label="sweep down (jumps up at left fold)")

    ax.set_xlabel(r"excitation frequency $\Omega$")
    ax.set_ylabel("response amplitude $A$")
    ax.set_title(f"Frequency response of the Duffing oscillator "
                 f"($\\zeta={ZETA}$, $\\alpha={ALPHA}$, $f={F_EXC}$)")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)


def plot_potentials(ax):
    """Linear, hardening and double-well potential energy."""
    x = np.linspace(-1.9, 1.9, 600)
    ax.plot(x, potential(x, OM0**2, 0.0), "k--", lw=1.4,
            label=r"linear, $\alpha=0$")
    ax.plot(x, potential(x, OM0**2, ALPHA), lw=1.6,
            label=fr"hardening, $\alpha={ALPHA}$")
    ax.plot(x, potential(x, -OM0D**2, ALPHA_D), lw=1.6, color="tab:purple",
            label=r"double well, $V=-\frac{1}{2}x^2+\frac{1}{4}x^4$")

    xs = OM0D / np.sqrt(ALPHA_D)      # position of the two minima
    ax.plot([-xs, xs], [potential(-xs, -OM0D**2, ALPHA_D)] * 2, "o",
            color="tab:purple", ms=6)
    ax.plot(0.0, 0.0, "x", color="crimson", ms=8, mew=2)
    ax.annotate("barrier", xy=(0.0, 0.0), xytext=(0.25, 0.35), fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.8))

    ax.set_xlabel("$x$")
    ax.set_ylabel("potential energy $V(x)$")
    ax.set_title("Potential energy landscape")
    ax.legend(loc="upper center", fontsize=8)
    ax.grid(alpha=0.3)


def plot_double_well_phase_portrait(ax):
    """Level sets of the conserved energy of x'' - x + x^3 = 0."""
    xg, vg = np.meshgrid(np.linspace(-2.0, 2.0, 400),
                         np.linspace(-1.4, 1.4, 400))
    energy = double_well_energy(xg, vg)

    levels = np.concatenate([np.linspace(-0.24, -0.02, 6),
                             np.linspace(0.05, 0.9, 7)])
    ax.contour(xg, vg, energy, levels=levels, colors="tab:blue",
               linewidths=0.8, linestyles="solid")
    # E = 0 is the homoclinic figure-eight separatrix through the saddle.
    ax.contour(xg, vg, energy, levels=[0.0], colors="crimson", linewidths=1.8)

    xs = OM0D / np.sqrt(ALPHA_D)
    ax.plot([-xs, xs], [0.0, 0.0], "o", color="k", ms=6)
    ax.plot(0.0, 0.0, "x", color="crimson", ms=9, mew=2)
    ax.text(0.05, 0.12, "saddle", color="crimson", fontsize=8)
    ax.text(xs - 0.25, -0.28, "centre", fontsize=8)
    ax.plot([], [], color="crimson", lw=1.8, label="separatrix $E=0$")
    ax.plot([], [], color="tab:blue", lw=0.8, label="orbits $E=$ const")

    ax.set_xlabel("$x$")
    ax.set_ylabel(r"$\dot{x}$")
    ax.set_title(r"Double well: $\ddot{x}-x+x^3=0$")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.3)


# =====================================================================
# MAIN
# =====================================================================

def main():
    # ---------------- analytical frequency response ----------------
    om_hb = np.linspace(OM_MIN, OM_MAX, N_HB)
    hb = hb_response_curve(om_hb)
    om_lo, om_hi = fold_frequencies(om_hb)

    print("Duffing oscillator: zeta=%.3f, om0=%.2f, alpha=%.2f, f=%.2f"
          % (ZETA, OM0, ALPHA, F_EXC))
    a_peak = hb[1].max()
    om_peak = hb[0][np.argmax(hb[1])]
    print("  peak amplitude (harmonic balance): A = %.4f at Omega = %.4f"
          % (a_peak, om_peak))
    print("  linear peak for alpha = 0 would be A = %.4f at Omega = %.4f"
          % (F_EXC / (2.0 * ZETA * OM0**2), OM0))
    if om_lo is not None:
        print("  three coexisting solutions for %.4f < Omega < %.4f"
              % (om_lo, om_hi))
    for om_test in (0.9, 1.3, 1.6):
        amps = hb_amplitudes(om_test)
        tag = ["stable" if hb_is_stable(a, om_test) else "UNSTABLE"
               for a in amps]
        print("  Omega = %.2f : A = %s  (%s)"
              % (om_test, np.array2string(amps, precision=4), ", ".join(tag)))

    # ---------------- numerical continuation ----------------
    om_sweep = np.linspace(OM_MIN, OM_MAX, N_SWEEP)
    up = sweep(om_sweep)                       # increasing Omega
    down = sweep(om_sweep[::-1])[::-1]         # decreasing Omega
    jump_up_idx = np.argmax(np.abs(np.diff(up)))
    jump_dn_idx = np.argmax(np.abs(np.diff(down)))
    print("  sweep up   jumps down near Omega = %.3f (dA = %.3f)"
          % (om_sweep[jump_up_idx], up[jump_up_idx + 1] - up[jump_up_idx]))
    print("  sweep down jumps up   near Omega = %.3f (dA = %.3f)"
          % (om_sweep[jump_dn_idx], down[jump_dn_idx] - down[jump_dn_idx + 1]))
    print("  hysteresis: the two sweeps differ by up to %.3f in amplitude"
          % np.max(np.abs(up - down)))

    # ---------------- figure ----------------
    fig = plt.figure(figsize=(11.5, 8.5))
    gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.25)
    plot_frequency_response(fig.add_subplot(gs[0, :]), om_hb, hb,
                            om_sweep, up, down)
    plot_potentials(fig.add_subplot(gs[1, 0]))
    plot_double_well_phase_portrait(fig.add_subplot(gs[1, 1]))
    fig.suptitle("1.2 Duffing Oscillator: jump phenomenon, potential "
                 "landscape and double-well phase portrait", fontsize=13)
    plt.show()


if __name__ == "__main__":
    main()
