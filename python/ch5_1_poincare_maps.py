"""
Chapter 5.1 -- Poincare Maps

Stroboscopic Poincare sections of the periodically forced Duffing oscillator.

Governing equation (plain text):

    x'' + 2*delta*x' + alpha*x + eps*x**3 = F*cos(Omega*t)

written as the first-order system  z = (x, x')  with

    x'  = v
    v'  = -2*delta*v - alpha*x - eps*x**3 + F*cos(Omega*t).

The right-hand side is T-periodic in t with T = 2*pi/Omega, so the stroboscopic
Poincare map is

    P(z0) = phi(T; z0),      Sigma = { z | t mod T = 0 },

i.e. "photograph the state once per forcing period".  A period-T response is a
single point in Sigma, a period-kT response gives k points, a chaotic response
gives a fractal cloud (the strange attractor).

What this script computes
-------------------------
1. Four stroboscopic sections of the same oscillator for four forcing
   amplitudes F, showing a period-1 orbit, the first period doubling (2 points),
   the second doubling (4 points) and a chaotic attractor.  The period of every
   case is *detected numerically* from the iterates, not assumed.
2. For the period-1 case: the fixed point z* of P is located with a Newton
   iteration, the Jacobian DP(z*) (= monodromy matrix) is built by centred
   finite differences, and its eigenvalues -- the Floquet multipliers -- are
   printed together with the Liouville check det(DP) = exp(-2*delta*T).
3. The largest Lyapunov exponent of each case, obtained by integrating the
   variational equation alongside the trajectory and renormalising once per
   period.  It is negative for the periodic cases and positive for the chaotic
   one, which is the quantitative criterion of the lecture.

The figure shows the four Poincare sections side by side; the title of each
panel reports F, the detected period and the measured Lyapunov exponent.

Why the twin-well system
------------------------
The oscillator used here is the classical two-well ("Holmes") Duffing
oscillator, alpha = -1, eps = +1, 2*delta = 0.3, Omega = 1.2.  The twin-well
form is not a matter of taste: the *hardening single-well* oscillator
(alpha = +1) has a single potential minimum and no saddle, hence no homoclinic
tangle, and it stays locked on a period-1 response for every forcing amplitude
in this range -- its largest Lyapunov exponent remains pinned at -delta.  Set
ALPHA = +1.0 below and the cascade disappears; that is worth trying once.

Runtime: about 30 s.  The chaotic panel dominates it; N_SAMPLE_CHAOS is the
knob -- doubling it doubles the runtime.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------
# PARAMETERS  (edit here)
# ----------------------------------------------------------------------
DELTA = 0.15          # damping ratio; the equation contains 2*delta*x'
ALPHA = -1.0          # linear stiffness (-1 = two-well potential, +1 = hardening)
EPS = 1.0             # cubic stiffness
OMEGA = 1.2           # forcing frequency  ->  T = 2*pi/OMEGA

# Forcing amplitudes to sample, with the expected response type.
# (F, expected period k, colour, number of stroboscopic points to plot)
CASES = [
    (0.20, 1, "#1f77b4"),
    (0.28, 2, "#2ca02c"),
    (0.29, 4, "#ff7f0e"),
    (0.50, None, "#d62728"),   # None = chaotic, no finite period expected
]

Z0 = (1.0, 0.0)       # initial condition: right-hand well, at rest
N_TRANSIENT = 200     # forcing periods discarded before sampling
N_SAMPLE = 400        # stroboscopic points kept for a periodic case
N_SAMPLE_CHAOS = 3000  # stroboscopic points kept for the chaotic case
N_LYAP = 400          # periods used for the Lyapunov exponent

RTOL, ATOL = 1e-9, 1e-12   # integrator tolerances (tight: chaos amplifies error)
METHOD = "DOP853"          # high-order explicit RK, cheapest at these tolerances

PERIOD_TOL = 1e-5     # tolerance of the numerical period detection in Sigma
PERIOD_MAX = 16       # largest period k the detector looks for
FD_STEP = 1e-6        # step of the centred difference for DP (slide: 1e-6..1e-5)

T_FORCE = 2.0 * np.pi / OMEGA     # forcing period T


# ----------------------------------------------------------------------
# Model
# ----------------------------------------------------------------------
def duffing_rhs(t, z, F):
    """Right-hand side of the forced Duffing oscillator, z = (x, v)."""
    x, v = z
    return [v, -2.0 * DELTA * v - ALPHA * x - EPS * x ** 3 + F * np.cos(OMEGA * t)]


def duffing_jacobian(t, z):
    """Jacobian df/dz of the flow. Only the stiffness entry depends on x,
    because the forcing enters additively and therefore drops out."""
    x = z[0]
    return np.array([[0.0, 1.0],
                     [-ALPHA - 3.0 * EPS * x ** 2, -2.0 * DELTA]])


# ----------------------------------------------------------------------
# Stroboscopic sampling: the Poincare map itself
# ----------------------------------------------------------------------
def stroboscopic_orbit(F, z0=Z0, n_transient=N_TRANSIENT, n_sample=N_SAMPLE):
    """Iterate the Poincare map by integrating the ODE once and reading the
    solution at t = n*T.

    The transient is integrated together with the sampled part; only the times
    n_transient*T ... (n_transient+n_sample-1)*T are requested via t_eval, which
    is the cheap equivalent of solve_ivp(dense_output=True) followed by
    sol.sol(...).  Returns an array of shape (2, n_sample).
    """
    n_end = n_transient + n_sample
    t_eval = (n_transient + np.arange(n_sample)) * T_FORCE
    sol = solve_ivp(duffing_rhs, [0.0, n_end * T_FORCE], list(z0), args=(F,),
                    t_eval=t_eval, method=METHOD, rtol=RTOL, atol=ATOL)
    if not sol.success:
        raise RuntimeError(f"integration failed for F = {F}: {sol.message}")
    return sol.y


def detect_period(points, tol=PERIOD_TOL, k_max=PERIOD_MAX):
    """Smallest k <= k_max with z_{n+k} = z_n for all sampled n, or None.

    A k-cycle of P corresponds to a period-kT orbit of the flow.  The test is
    run on the sampled part only, i.e. after the transient has died out.
    """
    x, v = points
    for k in range(1, k_max + 1):
        err = np.max(np.abs(x[k:] - x[:-k])) + np.max(np.abs(v[k:] - v[:-k]))
        if err < tol:
            return k
    return None


# ----------------------------------------------------------------------
# Fixed point, Jacobian DP and Floquet multipliers
# ----------------------------------------------------------------------
def poincare_map(z, F, k=1):
    """Apply the Poincare map k times: integrate exactly k forcing periods."""
    sol = solve_ivp(duffing_rhs, [0.0, k * T_FORCE], list(z), args=(F,),
                    method=METHOD, rtol=RTOL, atol=ATOL)
    return sol.y[:, -1]


def jacobian_fd(z, F, k=1, h=FD_STEP):
    """Centred finite-difference Jacobian DP^k(z), one column per state
    direction (two extra period integrations each), as on the slide."""
    dp = np.empty((2, 2))
    for j in range(2):
        e = np.zeros(2)
        e[j] = h
        dp[:, j] = (poincare_map(z + e, F, k) - poincare_map(z - e, F, k)) / (2.0 * h)
    return dp


def find_fixed_point(z_guess, F, k=1, tol=1e-11, itmax=25):
    """Newton iteration on P^k(z) - z = 0, using the finite-difference DP^k
    as the corrector matrix (the recipe given on the slide)."""
    z = np.asarray(z_guess, dtype=float)
    for _ in range(itmax):
        residual = poincare_map(z, F, k) - z
        if np.linalg.norm(residual) < tol:
            break
        jac = jacobian_fd(z, F, k) - np.eye(2)
        z = z - np.linalg.solve(jac, residual)
    return z, np.linalg.norm(poincare_map(z, F, k) - z)


# ----------------------------------------------------------------------
# Largest Lyapunov exponent
# ----------------------------------------------------------------------
def largest_lyapunov(F, z0=Z0, n_transient=N_TRANSIENT, n_periods=N_LYAP):
    """Largest Lyapunov exponent from the variational equation.

    The state is augmented with a tangent vector d obeying d' = J(t,z) d.  After
    every forcing period the tangent vector is rescaled to unit length and the
    logarithm of its growth factor accumulated; lambda_1 is the mean growth rate
    per unit time.  Renormalising is essential: otherwise d overflows (chaos) or
    underflows (stable orbit) within a few dozen periods.
    """
    def augmented(t, y, F):
        z = y[:2]
        d = y[2:]
        dz = duffing_rhs(t, z, F)
        dd = duffing_jacobian(t, z) @ d
        return [dz[0], dz[1], dd[0], dd[1]]

    # settle onto the attractor first
    z = poincare_map(np.asarray(z0, dtype=float), F, k=n_transient)

    d = np.array([1.0, 0.0])
    log_sum = 0.0
    for _ in range(n_periods):
        y0 = np.concatenate([z, d])
        sol = solve_ivp(augmented, [0.0, T_FORCE], y0, args=(F,),
                        method=METHOD, rtol=RTOL, atol=ATOL)
        y = sol.y[:, -1]
        z, d = y[:2], y[2:]
        growth = np.linalg.norm(d)
        log_sum += np.log(growth)
        d = d / growth                      # renormalise, direction is kept
    return log_sum / (n_periods * T_FORCE)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    print("Forced Duffing oscillator: x'' + 2*%.3g*x' + (%.3g)*x + %.3g*x**3 "
          "= F*cos(%.3g*t)" % (DELTA, ALPHA, EPS, OMEGA))
    print("forcing period T = %.6f\n" % T_FORCE)

    results = []
    for F, k_expected, colour in CASES:
        n_sample = N_SAMPLE_CHAOS if k_expected is None else N_SAMPLE
        pts = stroboscopic_orbit(F, n_sample=n_sample)
        k_found = detect_period(pts[:, :N_SAMPLE])   # detect on a short window
        lam = largest_lyapunov(F)
        results.append((F, k_expected, k_found, lam, pts, colour))

        label = f"period-{k_found}" if k_found else "aperiodic (chaotic)"
        print(f"F = {F:.3f}: {label:22s} lambda_1 = {lam:+.4f} "
              f"(expected {'chaos' if k_expected is None else 'period-' + str(k_expected)})")
        if k_expected is not None and k_found != k_expected:
            print("   WARNING: detected period differs from the expected one")
        if k_expected is None and (k_found is not None or lam <= 0.0):
            print("   WARNING: this case does not look chaotic")

    # --- Floquet analysis of the period-1 fixed point ------------------
    F1 = CASES[0][0]
    z_star, residual = find_fixed_point(results[0][4][:, 0], F1)
    dp = jacobian_fd(z_star, F1)
    mu = np.linalg.eigvals(dp)
    det_expected = np.exp(-2.0 * DELTA * T_FORCE)
    print("\nFixed point of P at F = %.3f" % F1)
    print("  z* = (%.6f, %.6f),  |P(z*)-z*| = %.2e" % (z_star[0], z_star[1], residual))
    print("  DP(z*) = [[%.6f, %.6f], [%.6f, %.6f]]"
          % (dp[0, 0], dp[0, 1], dp[1, 0], dp[1, 1]))
    print("  Floquet multipliers: mu1 = %s, mu2 = %s"
          % (np.array2string(mu[0], precision=5), np.array2string(mu[1], precision=5)))
    print("  |mu| = %.5f, %.5f  ->  %s"
          % (abs(mu[0]), abs(mu[1]),
             "stable" if max(abs(mu)) < 1.0 else "unstable"))
    print("  Liouville check: det(DP) = %.8f, exp(-2*delta*T) = %.8f"
          % (np.linalg.det(dp), det_expected))

    # --- Figure --------------------------------------------------------
    fig, axes = plt.subplots(1, 4, figsize=(15.0, 4.2), sharex=True, sharey=True,
                             constrained_layout=True)
    # All panels share one window (that of the widest, chaotic case) so that the
    # sequence 1 point -> 2 points -> 4 points -> cloud can be compared directly.
    all_pts = np.hstack([r[4] for r in results])
    pad_x = 0.1 * np.ptp(all_pts[0])
    pad_y = 0.1 * np.ptp(all_pts[1])
    for ax, (F, _, k_found, lam, pts, colour) in zip(axes, results):
        chaotic = k_found is None
        ax.scatter(pts[0], pts[1], s=(0.6 if chaotic else 45.0), c=colour,
                   marker="o", edgecolors="none", zorder=3,
                   label=("strange attractor (%d iterates)" % pts.shape[1]
                          if chaotic else "%d-cycle of P" % k_found))
        kind = "chaotic" if chaotic else "period-%d" % k_found
        ax.set_title("F = %.2f:  %s\n$\\lambda_1 = %+.3f$" % (F, kind, lam))
        ax.set_xlabel(r"$x(nT)$")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    axes[0].set_ylabel(r"$\dot{x}(nT)$")
    axes[0].set_xlim(all_pts[0].min() - pad_x, all_pts[0].max() + pad_x)
    axes[0].set_ylim(all_pts[1].min() - pad_y, all_pts[1].max() + pad_y)

    fig.suptitle(r"Poincare sections of the forced Duffing oscillator "
                 r"$\ddot{x}+2\delta\dot{x}+\alpha x+\varepsilon x^{3}"
                 r"=F\cos\Omega t$   "
                 r"($\delta=%.2f$, $\alpha=%.0f$, $\varepsilon=%.0f$, "
                 r"$\Omega=%.1f$), sampled at $t=nT$"
                 % (DELTA, ALPHA, EPS, OMEGA))
    plt.show()


if __name__ == "__main__":
    main()
