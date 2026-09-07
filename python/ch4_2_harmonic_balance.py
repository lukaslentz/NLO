"""
4.2 Harmonic Balance
====================

Multi-term harmonic balance (HB) for the damped, harmonically forced Duffing
oscillator

    x'' + 2*D*x' + x + eps*x**3 = f_hat*cos(Omega*t)          (omega_0 = 1)

The steady-state periodic response is written as a truncated Fourier series

    x(t) = sum_{k=1..N} ( a_k*cos(k*Omega*t) + b_k*sin(k*Omega*t) ),

inserted into the equation of motion, and the residual is projected onto every
retained basis function (Galerkin projection).  The cubic force is *not*
expanded by hand: it is evaluated in the time domain on a sample of the period
and transformed back with an FFT.  This is the Alternating Frequency/Time (AFT)
method, so exactly the same code works for any other nonlinearity.

The script produces four panels:

  (1) Frequency response.  The single-term (N = 1) amplitude equation
      (1 - Om^2 + 3*eps/4*a^2)^2 + (2*D*Om)^2 = (f_hat/a)^2
      is solved exactly for Om(a) and gives the classical bent resonance curve
      with its two folds.  On top of it the N-harmonic branch obtained by
      pseudo-arclength continuation is drawn, together with the amplitude of
      the fundamental measured from brute-force time integration during an
      upward and a downward frequency sweep -- the sweeps show the jump
      phenomenon and the hysteresis loop.
  (2) Waveform on the upper branch at a fixed excitation frequency: HB with
      N = 1, 3 and 7 harmonics against the numerically integrated steady state.
      The N = 1 error is essentially the discarded 3*Omega component.
  (3) Convergence: maximum waveform error of the N-term HB solution with
      respect to the time integration, plotted over N (logarithmic axis).
      The error drops geometrically and stagnates at the integration accuracy.
  (4) Amplitude spectrum of the converged solution: only odd harmonics are
      present, because the cubic force preserves the half-wave symmetry, and
      the amplitudes decay by roughly one order of magnitude per odd harmonic.

Runtime is roughly ten seconds; the time-integration sweeps dominate it.
Increasing N_SWEEP or N_PERIODS_TRANSIENT makes the script noticeably slower.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import root

# ----------------------------------------------------------------------
# Parameters -- change these to explore the model
# ----------------------------------------------------------------------
EPS = 0.3          # cubic stiffness coefficient (>0: hardening spring)
D = 0.05           # damping ratio (normalised, D << 1)
F_HAT = 0.2        # forcing amplitude
OM_MIN, OM_MAX = 0.5, 1.8   # frequency window of interest

N_HB = 7           # number of harmonics in the "converged" HB solution
N_LIST = [1, 2, 3, 4, 5, 6, 7, 9]   # truncation orders used in the study
N_AFT = 256        # time samples per period used by the AFT/FFT step

OM_WAVE = 1.20     # frequency for the waveform comparison (upper branch)

# time integration settings
N_PERIODS_TRANSIENT = 80    # periods discarded before the steady state is read
N_SWEEP = 35                # frequencies per sweep direction (drives the runtime)
RTOL, ATOL = 1e-9, 1e-11    # tolerances of the Runge-Kutta integrator
# the waveform used as the reference of the convergence study is integrated
# with much tighter tolerances, otherwise the integration error would mask the
# truncation error of the high-order HB solutions
RTOL_REF, ATOL_REF = 1e-13, 1e-14
N_PERIODS_REF = 300

# continuation settings
DS = 0.02          # arclength step
MAX_STEPS = 4000   # safety limit for the continuation loop


# ----------------------------------------------------------------------
# Harmonic balance with the AFT (Alternating Frequency/Time) method
# ----------------------------------------------------------------------
def hb_unpack(z):
    """Split the unknown vector z = [a_1, b_1, a_2, b_2, ...] into a, b."""
    return z[0::2], z[1::2]


def hb_waveform(z, om, t):
    """Evaluate the truncated Fourier series at the times t."""
    a, b = hb_unpack(z)
    k = np.arange(1, len(a) + 1)[:, None]
    ph = k * om * np.atleast_1d(t)[None, :]
    return (a[:, None] * np.cos(ph) + b[:, None] * np.sin(ph)).sum(axis=0)


def fourier_coeffs(g, n_harm):
    """Fourier cosine/sine coefficients of one period sampled at N_AFT points.

    g contains the samples of a periodic signal at theta = 2*pi*j/M.  The
    returned c[k], s[k] are the coefficients of cos(k*theta) and sin(k*theta)
    for k = 1..n_harm.  This is the "time -> frequency" half of the AFT step.
    """
    m = g.size
    gk = np.fft.rfft(g)[1:n_harm + 1]
    return 2.0 * gk.real / m, -2.0 * gk.imag / m


def hb_residual(z, om, n_harm):
    """Galerkin residual of the Duffing equation for the Fourier unknowns z.

    For every harmonic k the projections of
        x'' + 2*D*x' + x + eps*x^3 - f_hat*cos(Om*t)
    onto cos(k*Om*t) and sin(k*Om*t) must vanish.  The linear terms are
    differentiated analytically, the cubic force is evaluated in the time
    domain and transformed back (AFT), so no trigonometric expansion is needed.
    """
    a, b = hb_unpack(z)
    k = np.arange(1, n_harm + 1)

    # nonlinear force: frequency -> time -> nonlinearity -> frequency
    theta = 2.0 * np.pi * np.arange(N_AFT) / N_AFT
    ph = k[:, None] * theta[None, :]
    x = (a[:, None] * np.cos(ph) + b[:, None] * np.sin(ph)).sum(axis=0)
    c_nl, s_nl = fourier_coeffs(EPS * x ** 3, n_harm)

    lin = 1.0 - (k * om) ** 2          # inertia + linear stiffness
    res_cos = lin * a + 2.0 * D * om * k * b + c_nl
    res_sin = lin * b - 2.0 * D * om * k * a + s_nl
    res_cos[0] -= F_HAT                # forcing acts on the first cosine only

    out = np.empty_like(z)
    out[0::2] = res_cos
    out[1::2] = res_sin
    return out


def hb_solve(om, n_harm, guess):
    """Solve the HB equations at a fixed frequency, starting from `guess`.

    Convergence is judged by the residual itself, not by the solver flag: once
    the residual is at the round-off level the flag may report a stalled step
    even though the solution is perfectly converged.
    """
    sol = root(hb_residual, guess, args=(om, n_harm), method="hybr", tol=1e-13)
    ok = np.max(np.abs(hb_residual(sol.x, om, n_harm))) < 1e-10
    return sol.x, ok


def linear_guess(om, n_harm):
    """Initial guess from the linear frequency response (eps = 0)."""
    z = np.zeros(2 * n_harm)
    den = (1.0 - om ** 2) ** 2 + (2.0 * D * om) ** 2
    z[0] = F_HAT * (1.0 - om ** 2) / den
    z[1] = F_HAT * 2.0 * D * om / den
    return z


# ----------------------------------------------------------------------
# Pseudo-arclength continuation of the HB branch
# ----------------------------------------------------------------------
def continue_branch(n_harm, om_start, om_stop):
    """Trace the HB solution branch through the folds.

    The unknown vector is extended by the frequency, y = [z, Om].  A predictor
    step along the branch tangent plus a corrector that solves the HB equations
    together with the arclength constraint allows the continuation to pass the
    turning points, where a plain sweep in Om would fail.
    """
    z0, ok = hb_solve(om_start, n_harm, linear_guess(om_start, n_harm))
    if not ok:
        raise RuntimeError("continuation could not be started")
    dom = 1e-3
    z1, _ = hb_solve(om_start + dom, n_harm, z0)

    y_prev = np.append(z0, om_start)
    y_cur = np.append(z1, om_start + dom)
    ys = [y_prev, y_cur]

    for _ in range(MAX_STEPS):
        tang = y_cur - y_prev
        tang /= np.linalg.norm(tang)
        y_pred = y_cur + DS * tang

        def aug(y):
            r = hb_residual(y[:-1], y[-1], n_harm)
            return np.append(r, tang @ (y - y_cur) - DS)

        sol = root(aug, y_pred, method="hybr", tol=1e-12)
        if not sol.success:
            break
        y_prev, y_cur = y_cur, sol.x
        ys.append(y_cur)
        # stop once the branch has left the window on the high-frequency side
        if y_cur[-1] > om_stop and y_cur[-1] > y_prev[-1]:
            break
    return np.array(ys)


def fundamental_amplitude(z):
    """Amplitude sqrt(a_1^2 + b_1^2) of the fundamental harmonic."""
    return np.hypot(z[0], z[1])


def split_stability(om_branch):
    """Mark the segment between the two folds.

    Between the fold points the continuation runs backwards in Omega.  For the
    single-degree-of-freedom Duffing oscillator this middle segment is exactly
    the unstable branch (a rigorous check would use Floquet/Hill analysis of
    the periodic solution, see chapter 5.2).
    """
    going_up = np.diff(om_branch) > 0
    return np.append(going_up, going_up[-1])


# ----------------------------------------------------------------------
# Reference solution: direct time integration
# ----------------------------------------------------------------------
def duffing_rhs(t, y, om):
    """State-space form of the forced Duffing oscillator."""
    x, v = y
    return [v, -2.0 * D * v - x - EPS * x ** 3 + F_HAT * np.cos(om * t)]


def steady_state(om, y0, n_transient=N_PERIODS_TRANSIENT,
                 rtol=RTOL, atol=ATOL, n_samples=512):
    """Integrate away the transient and return one sampled period.

    Returns (t, x, y_end): the times and displacements of one period after the
    transient, plus the final state (used as the initial condition of the next
    frequency in a sweep, which is what produces the hysteresis).  The period
    is sampled at n_samples+1 equidistant points, so the last point repeats the
    first one and an FFT can use the first n_samples values.
    """
    period = 2.0 * np.pi / om
    t0 = n_transient * period
    t_out = t0 + np.linspace(0.0, period, n_samples + 1)
    sol = solve_ivp(duffing_rhs, [0.0, t_out[-1]], y0, args=(om,),
                    t_eval=t_out, rtol=rtol, atol=atol, method="DOP853")
    return sol.t, sol.y[0], sol.y[:, -1]


def measured_fundamental(t, x, om):
    """Amplitude of the cos/sin(Om*t) component of a sampled periodic signal."""
    # trapezoidal projection over exactly one period (t[-1] repeats t[0])
    w = np.ones_like(t)
    w[0] = w[-1] = 0.5
    tt = t - t[0]
    a = 2.0 / (len(t) - 1) * np.sum(w * x * np.cos(om * tt))
    b = 2.0 / (len(t) - 1) * np.sum(w * x * np.sin(om * tt))
    return np.hypot(a, b)


def frequency_sweep(om_values, y0):
    """Slow frequency sweep: each frequency starts from the previous state."""
    amps = []
    y = np.array(y0, dtype=float)
    for om in om_values:
        t, x, y = steady_state(om, y)
        amps.append(measured_fundamental(t, x, om))
    return np.array(amps)


# ----------------------------------------------------------------------
# Single-term amplitude equation, solved exactly for Omega(a)
# ----------------------------------------------------------------------
def single_term_curve(a_max=3.0, n=800):
    """Exact solution Om(a) of the N = 1 amplitude equation.

    With u = Om^2 and K = 1 + 3*eps/4*a^2 the amplitude equation
    (K - u)^2 + 4*D^2*u = (f_hat/a)^2 is a quadratic in u; both roots are
    kept, so the full multi-valued response curve is obtained.
    """
    a = np.linspace(1e-3, a_max, n)
    k_eff = 1.0 + 0.75 * EPS * a ** 2
    p = 4.0 * D ** 2 - 2.0 * k_eff
    q = k_eff ** 2 - (F_HAT / a) ** 2
    disc = p ** 2 - 4.0 * q
    ok = disc >= 0.0
    sq = np.sqrt(np.where(ok, disc, 0.0))
    u_lo, u_hi = (-p - sq) / 2.0, (-p + sq) / 2.0
    branch_lo = np.where(ok & (u_lo > 0), np.sqrt(np.abs(u_lo)), np.nan)
    branch_hi = np.where(ok & (u_hi > 0), np.sqrt(np.abs(u_hi)), np.nan)
    return a, branch_lo, branch_hi


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    print("4.2 Harmonic Balance -- Duffing oscillator")
    print(f"    eps = {EPS}, D = {D}, f_hat = {F_HAT}\n")

    # --- (1) frequency response ---------------------------------------
    a_st, om_lo, om_hi = single_term_curve()

    branch = continue_branch(N_HB, OM_MIN, OM_MAX)
    om_br = branch[:, -1]
    a_br = np.array([fundamental_amplitude(y[:-1]) for y in branch])
    stable = split_stability(om_br)

    # fold points of the multi-harmonic branch
    turn = np.where(np.diff(np.sign(np.diff(om_br))) != 0)[0] + 1
    folds = [(om_br[i], a_br[i]) for i in turn]
    for om_f, a_f in folds:
        print(f"    fold point (N={N_HB}): Omega = {om_f:.4f}, a1 = {a_f:.4f}")

    om_up = np.linspace(OM_MIN, OM_MAX, N_SWEEP)
    amp_up = frequency_sweep(om_up, [0.0, 0.0])
    amp_down = frequency_sweep(om_up[::-1], [0.0, 0.0])
    om_down = om_up[::-1]

    # --- (2)/(3) waveform and convergence at OM_WAVE ------------------
    # start high to land on the large-amplitude branch
    t_w, x_w, _ = steady_state(OM_WAVE, [2.5, 0.0], n_transient=N_PERIODS_REF,
                               rtol=RTOL_REF, atol=ATOL_REF)
    a_num = measured_fundamental(t_w, x_w, OM_WAVE)
    print(f"\n    at Omega = {OM_WAVE}: time integration gives "
          f"a1 = {a_num:.4f}, max|x| = {np.max(np.abs(x_w)):.4f}")

    errors, waveforms = [], {}
    guess = None
    for n in N_LIST:
        g = np.zeros(2 * n)
        if guess is None:
            g[0] = a_num          # cosine-dominated start on the upper branch
        else:
            g[:min(len(guess), 2 * n)] = guess[:min(len(guess), 2 * n)]
        z, ok = hb_solve(OM_WAVE, n, g)
        if not ok:
            raise RuntimeError(f"HB solve failed for N = {n}")
        guess = z
        x_hb = hb_waveform(z, OM_WAVE, t_w)
        err = np.max(np.abs(x_hb - x_w))
        errors.append(err)
        waveforms[n] = x_hb
        print(f"    N = {n:2d}: a1 = {fundamental_amplitude(z):.6f}, "
              f"max waveform error = {err:.3e}")

    z_conv, _ = hb_solve(OM_WAVE, N_HB, guess[:2 * N_HB])
    a_k, b_k = hb_unpack(z_conv)
    harm = np.hypot(a_k, b_k)
    print("    harmonic amplitudes (N = %d): " % N_HB
          + ", ".join(f"{h:.3e}" for h in harm))

    # spectrum of the numerical steady state for comparison
    xf = np.fft.rfft(x_w[:-1]) * 2.0 / (len(x_w) - 1)
    harm_num = np.abs(xf[1:N_HB + 1])

    # --- plots ---------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5))

    ax = axes[0, 0]
    # the N = 1 curve is drawn thick and in the background: it deviates from
    # the multi-harmonic branch only slightly, mainly near the upper fold
    ax.plot(om_lo, a_st, "c-", lw=4.0, alpha=0.55, zorder=1,
            label="single-term HB (N = 1)")
    ax.plot(om_hi, a_st, "c-", lw=4.0, alpha=0.55, zorder=1)
    ax.plot(om_br[stable], a_br[stable], "b.", ms=3,
            label=f"HB N = {N_HB} (stable)")
    ax.plot(om_br[~stable], a_br[~stable], "r.", ms=3,
            label=f"HB N = {N_HB} (unstable)")
    a_bb = np.linspace(0.0, 3.0, 200)
    ax.plot(np.sqrt(1.0 + 0.75 * EPS * a_bb ** 2), a_bb, "k--", lw=1.0,
            label=r"backbone $\Omega^2=1+\frac{3\varepsilon}{4}\hat a^2$")
    ax.plot(om_up, amp_up, "k^", ms=4, mfc="none",
            label=r"time integration, $\Omega\nearrow$")
    ax.plot(om_down, amp_down, "kv", ms=4, mfc="none",
            label=r"time integration, $\Omega\searrow$")
    ax.axvline(OM_WAVE, color="0.7", lw=0.8, ls=":")
    ax.set_xlim(OM_MIN, OM_MAX)
    ax.set_ylim(0.0, 3.0)
    ax.set_xlabel(r"excitation frequency $\Omega$")
    ax.set_ylabel(r"fundamental amplitude $\hat a_1$")
    ax.set_title("Frequency response, jump phenomenon and hysteresis")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    tt = (t_w - t_w[0]) * OM_WAVE / (2.0 * np.pi)
    ax.plot(tt, x_w, "k-", lw=2.0, label="time integration")
    for n, style in zip([1, 3, 7], ["r--", "g-.", "b:"]):
        if n in waveforms:
            ax.plot(tt, waveforms[n], style, lw=1.5, label=f"HB, N = {n}")
    ax.set_xlabel(r"$t/T$  (one period, $T=2\pi/\Omega$)")
    ax.set_ylabel(r"$x(t)$")
    ax.set_title(rf"Steady-state waveform at $\Omega = {OM_WAVE}$")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)
    # the truncation error is far too small to be seen in the waveform itself,
    # so the deviation is shown separately in an inset
    axin = ax.inset_axes((0.12, 0.06, 0.42, 0.30))
    for n, style in zip([1, 3, 5], ["r--", "g-.", "b:"]):
        axin.semilogy(tt, np.abs(waveforms[n] - x_w) + 1e-16, style, lw=1.0,
                      label=f"N = {n}")
    axin.set_ylim(1e-8, 1.0)
    axin.set_title(r"$|x_{\rm HB}-x_{\rm num}|$", fontsize=7)
    axin.tick_params(labelsize=6)
    axin.legend(fontsize=6)
    axin.grid(alpha=0.3)

    ax = axes[1, 0]
    ax.semilogy(N_LIST, errors, "bo-")
    ax.set_xlabel("number of retained harmonics $N$")
    ax.set_ylabel(r"$\max_t\,|x_{\rm HB}-x_{\rm num}|$")
    ax.set_title("Convergence of the truncated Fourier ansatz")
    ax.grid(alpha=0.3, which="both")
    ax.annotate("even $N$ add nothing:\nonly odd harmonics exist",
                xy=(2, errors[1]), xytext=(3.6, errors[0] * 0.25),
                fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.8))
    ax.set_ylim(min(errors) * 0.2, max(errors) * 5.0)

    ax = axes[1, 1]
    ks = np.arange(1, N_HB + 1)
    floor = 1e-16
    ax.bar(ks - 0.18, np.maximum(harm, floor), width=0.36,
           label=f"HB, N = {N_HB}")
    ax.bar(ks + 0.18, np.maximum(harm_num, floor), width=0.36,
           label="time integration")
    ax.set_yscale("log")
    ax.set_ylim(1e-12, 10.0)
    ax.set_xticks(ks)
    ax.set_xlabel(r"harmonic order $k$ (frequency $k\Omega$)")
    ax.set_ylabel(r"amplitude $\sqrt{a_k^2+b_k^2}$")
    ax.set_title(rf"Amplitude spectrum at $\Omega = {OM_WAVE}$")
    ax.text(0.55, 0.06, "even harmonics vanish\n(half-wave symmetry of $x^3$)",
            transform=ax.transAxes, ha="center", fontsize=8, color="0.3",
            bbox=dict(fc="white", ec="0.8", alpha=0.9))
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y", which="both")

    fig.suptitle("4.2 Harmonic Balance -- forced Duffing oscillator "
                 rf"($\varepsilon={EPS}$, $D={D}$, $\hat f={F_HAT}$)")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    plt.show()


if __name__ == "__main__":
    main()
