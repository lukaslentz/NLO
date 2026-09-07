"""
Chapter 6.2 -- Chaos and Fractal Basins.

This script demonstrates the three central statements of the chapter on a single
mechanical system, the harmonically forced Duffing oscillator with a twin-well
(negative linear stiffness) restoring force

    x'' + 2*D*x' - x + eps*x**3 = fhat*cos(Omega*t)

written as the first-order system  x' = v,  v' = fhat*cos(Omega*t) - 2*D*v + x - eps*x**3.

Four panels are produced:

(a) The largest Lyapunov exponent lambda_1 as a function of the forcing amplitude
    fhat, computed with the Benettin renormalization algorithm: integrate a
    reference and a perturbed trajectory over one forcing period T, accumulate
    ln(|delta|/|delta_0|), and rescale the separation vector back to |delta_0|.
    lambda_1 changes sign at the onset of chaos; in the regular windows it
    returns to the value -D dictated by the damping (a periodic attractor of a
    system with phase-space contraction rate -2D has exponents 0 and -2D, and
    the period-averaged perturbation decays like exp(-D t) when the marginal
    direction is not exactly hit).
(b) Sensitive dependence on initial conditions in the time domain: two histories
    x(t) whose initial states differ by delta_0 = 1e-8 stay visually identical
    for tens of periods and then decorrelate completely.
(c) The same pair of trajectories on a logarithmic separation axis. The straight
    pre-saturation segment has slope lambda_1; a least-squares fit over that
    window reproduces the Benettin value. After saturation at the attractor
    extent no further information about lambda_1 can be extracted.
(d) The basin of attraction of the two coexisting period-one attractors (one in
    each potential well) at a sub-chaotic forcing amplitude. The boundary
    between the two basins is fractal: no matter how finely the grid is
    resolved, blue and red points remain interleaved, so a finite measurement
    accuracy of the initial state can never decide the long-term outcome.

Deviation from the compressed listing on the slide: the listing uses the
single-well Duffing oscillator (+x instead of -x) at Omega = 0.8. That system is
never chaotic for the amplitudes swept there -- its Lyapunov exponent stays
pinned at -D = -0.08 for every fhat up to 1.0 (verified numerically). The
chapter's phenomena -- positive lambda_1, coexisting attractors, fractal basin
boundaries -- live in the twin-well version, which is the classical Duffing /
Holmes chaos model. The damping D = 0.08 and the cubic coefficient eps = 1 of
the slide are kept, the forcing frequency is raised to Omega = 1.2 so that the
chaotic window falls inside the slide's amplitude sweep 0.10 <= fhat <= 0.45.
"""

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# PARAMETERS  (everything a student may want to change lives here)
# ----------------------------------------------------------------------------
D = 0.08            # viscous damping ratio, equation of motion uses 2*D*x'
EPS = 1.0           # coefficient of the cubic restoring term
OMEGA = 1.2         # forcing frequency
T_FORCE = 2.0 * np.pi / OMEGA          # forcing period

# --- (a) Lyapunov sweep -------------------------------------------------------
F_MIN, F_MAX = 0.10, 0.45   # forcing-amplitude range of the slide listing
N_F = 71                    # number of amplitudes (refine: 141, 281, ... -- the
                            # sweep is vectorised over fhat, so the cost grows
                            # only very weakly with N_F)
N_TRANSIENT = 100           # forcing periods discarded before measuring
N_BENETTIN = 150            # renormalization steps used for lambda_1
                            # (refine: 400 periods gives a smoother curve at
                            # roughly 2.5x the runtime of this panel)
N_SUB = 160                 # RK4 sub-steps per forcing period (refine: 320)
DELTA0 = 1.0e-8             # initial separation of the Benettin pair
Y0 = (1.0, 0.0)             # starting state, inside the right-hand well

# --- (b), (c) sensitive dependence -------------------------------------------
F_CHAOS = 0.35              # amplitude inside the chaotic window
N_SDIC = 40                 # forcing periods shown in the time histories
N_SDIC_TRANSIENT = 60       # periods run first, so both start on the attractor
FIT_WINDOW = (5.0, 60.0)    # time window (in units of t) used for the slope fit
                            # -- must stay below the saturation level

# --- (d) basin of attraction --------------------------------------------------
F_BASIN = 0.15              # amplitude with two coexisting periodic attractors
N_GRID = 120                # grid points per axis (refine: 240 or 480 -- the
                            # cost grows with N_GRID**2, 240 already needs about
                            # two minutes; the fractal structure only gets richer)
X_RANGE = (-2.0, 2.0)       # initial displacement range
V_RANGE = (-2.0, 2.0)       # initial velocity range
N_BASIN_PERIODS = 70        # periods integrated before the outcome is read off
                            # (chaotic transients near the boundary are long;
                            # fewer periods leave many points unclassified)
N_BASIN_SUB = 64            # RK4 sub-steps per period for the basin scan
N_SAMPLE = 8                # Poincare samples used to classify the attractor
SCATTER_TOL = 0.15          # spread of those samples above which a point counts
                            # as "not settled" (still transient / chaotic)


# ----------------------------------------------------------------------------
# VECTORISED FIXED-STEP INTEGRATOR
# ----------------------------------------------------------------------------
# All three experiments integrate many copies of the same non-autonomous ODE at
# once (many forcing amplitudes, or many initial conditions). A fixed-step RK4
# on numpy arrays is used instead of scipy.integrate.solve_ivp for two reasons:
#   * speed -- tens of thousands of trajectories advance in lock-step;
#   * reproducibility -- an adaptive solver picks different step sequences for
#     two neighbouring initial states, which injects a numerical difference of
#     the order of the tolerance. With a separation of only 1e-8 that noise
#     would contaminate the Lyapunov estimate. A fixed step applies exactly the
#     same map to every trajectory, so the only difference is the physical one.

def duffing_rhs(t, x, v, fhat):
    """Right-hand side of the twin-well Duffing oscillator (arrays allowed).

    Returns (dx/dt, dv/dt) for  x'' + 2 D x' - x + eps x^3 = fhat cos(Omega t).
    The linear stiffness is negative: the potential  V(x) = -x^2/2 + eps x^4/4
    has two wells at x = +/- sqrt(1/eps) and an unstable saddle at x = 0.
    """
    return v, fhat * np.cos(OMEGA * t) - 2.0 * D * v + x - EPS * x**3


def rk4_step(t, x, v, dt, fhat):
    """One classical Runge-Kutta step of the forced Duffing oscillator."""
    k1x, k1v = duffing_rhs(t, x, v, fhat)
    k2x, k2v = duffing_rhs(t + 0.5 * dt, x + 0.5 * dt * k1x, v + 0.5 * dt * k1v, fhat)
    k3x, k3v = duffing_rhs(t + 0.5 * dt, x + 0.5 * dt * k2x, v + 0.5 * dt * k2v, fhat)
    k4x, k4v = duffing_rhs(t + dt, x + dt * k3x, v + dt * k3v, fhat)
    x_new = x + dt / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x)
    v_new = v + dt / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v)
    return x_new, v_new


def advance(t, x, v, dt, n_steps, fhat):
    """Advance a whole ensemble by n_steps RK4 steps; returns (t, x, v)."""
    for _ in range(n_steps):
        x, v = rk4_step(t, x, v, dt, fhat)
        t += dt
    return t, x, v


# ----------------------------------------------------------------------------
# (a) LARGEST LYAPUNOV EXPONENT -- BENETTIN RENORMALIZATION
# ----------------------------------------------------------------------------

def lyapunov_sweep(fhat_values):
    """Largest Lyapunov exponent for an array of forcing amplitudes.

    Benettin algorithm, renormalising once per forcing period tau = T:
      1. integrate reference and perturbed state over tau,
      2. measure d = ||delta(tau)|| in the (x, v) phase plane,
      3. accumulate ln(d / delta_0),
      4. pull the perturbed state back to distance delta_0 along the same
         direction and continue.
    After N steps  lambda_1 = (1 / (N*tau)) * sum_k ln(d_k / delta_0).

    Renormalisation is essential: without it the separation would saturate at
    the diameter of the attractor and the average slope would be biased towards
    zero. Rescaling keeps the perturbation infinitesimal while its direction is
    free to align with the most unstable direction of the flow.

    All amplitudes are integrated simultaneously as one numpy ensemble.
    """
    fhat = np.atleast_1d(np.asarray(fhat_values, dtype=float))
    dt = T_FORCE / N_SUB

    # Discard the transient so that the pair sits on the attractor.
    x = np.full(fhat.shape, Y0[0])
    v = np.full(fhat.shape, Y0[1])
    t, x, v = advance(0.0, x, v, dt, N_TRANSIENT * N_SUB, fhat)

    # Perturbed copy, offset purely in displacement.
    xp, vp = x + DELTA0, v.copy()

    log_sum = np.zeros(fhat.shape)
    for _ in range(N_BENETTIN):
        t_ref = t
        t, x, v = advance(t_ref, x, v, dt, N_SUB, fhat)
        _, xp, vp = advance(t_ref, xp, vp, dt, N_SUB, fhat)
        dx, dv = xp - x, vp - v
        dist = np.hypot(dx, dv)
        log_sum += np.log(dist / DELTA0)
        scale = DELTA0 / dist                 # renormalise, keep the direction
        xp, vp = x + dx * scale, v + dv * scale

    return log_sum / (N_BENETTIN * T_FORCE)


# ----------------------------------------------------------------------------
# (b), (c) SENSITIVE DEPENDENCE ON INITIAL CONDITIONS
# ----------------------------------------------------------------------------

def sdic_pair(fhat, n_periods):
    """Two trajectories separated by DELTA0 at t = 0, without renormalisation.

    Returns (t, x_ref, x_pert, separation) sampled at every RK4 step, with the
    time origin placed after a transient so that both start on the attractor.
    This is the raw experiment behind the definition of lambda_1: the log of the
    separation grows linearly until it saturates at the attractor extent.
    """
    dt = T_FORCE / N_SUB
    t, x, v = advance(0.0, np.array(Y0[0]), np.array(Y0[1]), dt,
                      N_SDIC_TRANSIENT * N_SUB, fhat)

    # Restart the clock at a multiple of T, so cos(Omega t) is unchanged.
    t = 0.0
    xp, vp = x + DELTA0, v.copy()

    n_steps = n_periods * N_SUB
    ts = np.empty(n_steps + 1)
    xs = np.empty(n_steps + 1)
    xps = np.empty(n_steps + 1)
    sep = np.empty(n_steps + 1)
    ts[0], xs[0], xps[0] = t, x, xp
    sep[0] = abs(xp - x)
    for k in range(1, n_steps + 1):
        x, v = rk4_step(t, x, v, dt, fhat)
        xp, vp = rk4_step(t, xp, vp, dt, fhat)
        t += dt
        ts[k], xs[k], xps[k] = t, x, xp
        sep[k] = np.hypot(xp - x, vp - v)
    return ts, xs, xps, sep


def fit_divergence_slope(t, separation, window):
    """Least-squares slope of ln|delta| over a pre-saturation time window.

    Only the straight segment may be used: before it the perturbation has not
    yet aligned with the unstable direction, after it the separation saturates
    at the attractor diameter and carries no exponent information any more.
    """
    mask = (t >= window[0]) & (t <= window[1])
    slope, intercept = np.polyfit(t[mask], np.log(separation[mask]), 1)
    return slope, intercept, mask


# ----------------------------------------------------------------------------
# (d) BASIN OF ATTRACTION
# ----------------------------------------------------------------------------

def basin_of_attraction(fhat):
    """Classify a grid of initial conditions by the attractor they end on.

    Every grid point is integrated for N_BASIN_PERIODS forcing periods, then
    sampled stroboscopically (once per period, a Poincare section) N_SAMPLE
    times. For a settled period-one attractor those samples coincide; their
    mean sign says which potential well the motion ended in. A large spread
    means the trajectory has not settled within the integration time -- either a
    long chaotic transient sustained by the chaotic saddle in the boundary, or a
    non-period-one attractor.

    Returns (labels, x_grid, v_grid) with labels in {-1: left well,
    +1: right well, 0: not settled}.
    """
    dt = T_FORCE / N_BASIN_SUB
    x0 = np.linspace(X_RANGE[0], X_RANGE[1], N_GRID)
    v0 = np.linspace(V_RANGE[0], V_RANGE[1], N_GRID)
    X0, V0 = np.meshgrid(x0, v0)
    x, v = X0.ravel().copy(), V0.ravel().copy()

    t, x, v = advance(0.0, x, v, dt, N_BASIN_PERIODS * N_BASIN_SUB, fhat)

    samples = np.empty((N_SAMPLE, x.size))
    for k in range(N_SAMPLE):
        t, x, v = advance(t, x, v, dt, N_BASIN_SUB, fhat)
        samples[k] = x                        # stroboscopic sample, phase fixed

    mean_x = samples.mean(axis=0)
    spread = samples.std(axis=0)
    labels = np.where(mean_x > 0.0, 1.0, -1.0)
    labels = np.where(spread > SCATTER_TOL, 0.0, labels)
    return labels.reshape(N_GRID, N_GRID), x0, v0


def boundary_fraction(labels):
    """Fraction of grid cells that have a differently labelled neighbour.

    This is the finite-resolution measure behind the uncertainty exponent: for a
    smooth boundary it shrinks like the cell size, for a fractal boundary it
    shrinks much more slowly (like eps**alpha with alpha = n - d_b < 1), so a
    large value at a fine grid is the signature of a fractal boundary.
    """
    inner = labels[1:-1, 1:-1]
    differs = ((inner != labels[:-2, 1:-1]) | (inner != labels[2:, 1:-1]) |
               (inner != labels[1:-1, :-2]) | (inner != labels[1:-1, 2:]))
    return differs.mean()


# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------

def main():
    # --- (a) Lyapunov exponent versus forcing amplitude ----------------------
    fhat = np.linspace(F_MIN, F_MAX, N_F)
    lam = lyapunov_sweep(fhat)

    chaotic = lam > 0.0
    if np.any(chaotic):
        f_onset = fhat[np.argmax(chaotic)]
        print(f"chaos onset at fhat ~ {f_onset:.3f} "
              f"(first amplitude with lambda_1 > 0)")
        print(f"max lambda_1 = {lam.max():+.4f} 1/s at fhat = "
              f"{fhat[np.argmax(lam)]:.3f}")
    print(f"min lambda_1 = {lam.min():+.4f} 1/s "
          f"(regular windows, expected about -D = {-D:+.3f})")

    # --- (b), (c) sensitive dependence at one chaotic amplitude --------------
    t_sd, x_ref, x_pert, sep = sdic_pair(F_CHAOS, N_SDIC)
    slope, intercept, fit_mask = fit_divergence_slope(t_sd, sep, FIT_WINDOW)
    lam_benettin = float(lyapunov_sweep([F_CHAOS])[0])
    print(f"fhat = {F_CHAOS:.2f}:  lambda_1 (Benettin)     = {lam_benettin:+.4f}")
    print(f"fhat = {F_CHAOS:.2f}:  lambda_1 (log-fit slope) = {slope:+.4f}")

    # Prediction horizon: the time an initial uncertainty delta_0 needs to grow
    # to the attractor extent A, t_max = ln(A/delta_0) / lambda_1.
    extent = x_ref.max() - x_ref.min()
    t_horizon = np.log(extent / DELTA0) / slope
    print(f"attractor extent A = {extent:.2f}, prediction horizon "
          f"t_max = {t_horizon:.1f} ({t_horizon / T_FORCE:.1f} forcing periods)")

    # --- (d) basin of attraction --------------------------------------------
    labels, x0, v0 = basin_of_attraction(F_BASIN)
    frac_right = np.mean(labels > 0.5)
    frac_left = np.mean(labels < -0.5)
    frac_open = np.mean(labels == 0.0)
    print(f"fhat = {F_BASIN:.2f}: basin fractions  right well {frac_right:.2f}, "
          f"left well {frac_left:.2f}, not settled {frac_open:.2f}")
    print(f"cells adjacent to a basin boundary: {boundary_fraction(labels):.3f} "
          f"at {N_GRID}x{N_GRID} resolution")

    # ------------------------------------------------------------------ plots
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5))
    fig.suptitle("6.2 Chaos and fractal basins -- twin-well Duffing oscillator\n"
                 rf"$\ddot x + 2D\dot x - x + \varepsilon x^3 = "
                 rf"\hat f\cos\Omega t$,  $D={D}$, $\varepsilon={EPS}$, "
                 rf"$\Omega={OMEGA}$", fontsize=12)

    # (a) Lyapunov exponent
    ax = axes[0, 0]
    ax.plot(fhat, lam, "-", color="C0", lw=1.2, label=r"$\lambda_1$ (Benettin)")
    ax.plot(fhat[lam > 0], lam[lam > 0], ".", color="C3", ms=5,
            label=r"chaotic, $\lambda_1>0$")
    ax.axhline(0.0, ls="--", color="k", lw=0.8)
    ax.axhline(-D, ls=":", color="gray", lw=0.8,
               label=rf"$-D={-D:.2f}$ (periodic attractor)")
    ax.axvline(F_CHAOS, ls="-.", color="C2", lw=0.8,
               label=rf"$\hat f={F_CHAOS}$ (panels b, c)")
    ax.axvline(F_BASIN, ls="-.", color="C1", lw=0.8,
               label=rf"$\hat f={F_BASIN}$ (panel d)")
    ax.set_xlabel(r"forcing amplitude $\hat f$")
    ax.set_ylabel(r"largest Lyapunov exponent $\lambda_1$")
    ax.set_title("(a) Benettin renormalization sweep")
    ax.set_ylim(top=lam.max() + 0.14)         # headroom for the legend box
    ax.legend(fontsize=7.5, loc="upper left", framealpha=0.95)
    ax.grid(alpha=0.3)

    # (b) sensitive dependence in the time domain
    ax = axes[0, 1]
    ax.plot(t_sd, x_ref, color="C0", lw=0.8, label=r"$x(t)$ from $x_0$")
    ax.plot(t_sd, x_pert, color="C3", lw=0.8, ls="--",
            label=rf"$x(t)$ from $x_0+\delta_0$, $\delta_0=10^{{-8}}$")
    ax.set_xlabel(r"time $t$")
    ax.set_ylabel(r"displacement $x$")
    ax.set_title(rf"(b) SDIC in the time domain, $\hat f={F_CHAOS}$")
    ax.axvline(np.log(extent / DELTA0) / slope, color="gray", ls=":", lw=0.9)
    ax.annotate("records decorrelate", xy=(t_horizon, 1.55),
                xytext=(t_horizon + 6, 1.9), fontsize=8, color="gray")
    ax.set_ylim(-2.1, 2.6)
    ax.legend(fontsize=8, loc="lower left", ncol=2, framealpha=0.95)
    ax.grid(alpha=0.3)

    # (c) logarithmic divergence and its slope
    ax = axes[1, 0]
    ax.semilogy(t_sd, sep, color="C0", lw=0.9, label=r"$\|\delta(t)\|$")
    ax.semilogy(t_sd[fit_mask], np.exp(intercept + slope * t_sd[fit_mask]),
                color="C3", lw=1.6, ls="--",
                label=rf"fit: $\lambda_1={slope:.3f}$")
    ax.axhline(extent, color="gray", ls=":", lw=0.9,
               label="attractor extent (saturation)")
    ax.axvspan(FIT_WINDOW[0], FIT_WINDOW[1], color="C2", alpha=0.08)
    ax.set_xlabel(r"time $t$")
    ax.set_ylabel(r"separation $\|\delta(t)\|$")
    ax.set_title("(c) Exponential separation, slope = $\\lambda_1$")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3, which="both")

    # (d) basin of attraction
    ax = axes[1, 1]
    ax.imshow(labels, origin="lower", cmap="coolwarm", vmin=-1.0, vmax=1.0,
              extent=[x0[0], x0[-1], v0[0], v0[-1]], aspect="auto",
              interpolation="nearest")
    ax.set_xlabel(r"initial displacement $x_0$")
    ax.set_ylabel(r"initial velocity $\dot x_0$")
    ax.set_title(rf"(d) Basins of attraction, $\hat f={F_BASIN}$ "
                 f"({N_GRID}x{N_GRID} grid)")
    # Legend by proxy patches -- imshow itself produces no legend handles.
    cmap = plt.get_cmap("coolwarm")
    handles = [plt.Line2D([], [], marker="s", ls="", color=cmap(1.0),
                          label=r"$\to$ right well ($x>0$)"),
               plt.Line2D([], [], marker="s", ls="", color=cmap(0.0),
                          label=r"$\to$ left well ($x<0$)"),
               plt.Line2D([], [], marker="s", ls="", color=cmap(0.5),
                          label="not settled (transient chaos)")]
    ax.legend(handles=handles, fontsize=8, loc="upper right", framealpha=0.9)

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    plt.show()


if __name__ == "__main__":
    main()
