"""
Chapter 4.5 -- Method of Multiple Scales
========================================

The method of multiple scales (MMS) splits the dynamics into a fast carrier
oscillation on T0 = t and a slow modulation on T1 = eps*t.  For the weakly
damped, weakly nonlinear, harmonically forced Duffing oscillator near primary
resonance the solvability condition at order eps yields the *modulation*
(slow-flow) equations for the real amplitude a and the relative phase gamma:

    da/dT1       = -mu*a + (F/2)*sin(gamma)
    a*dgamma/dT1 =  sigma*a - (3/8)*a**3 + (F/2)*cos(gamma)

(the slides use the opposite sign of both forcing terms, i.e. gamma shifted by
pi; the steady states and their stability are of course the same)

This script does four things with them:

  (a) It solves the steady-state condition, the cubic frequency-response
      equation  [mu**2 + (sigma - 3*a**2/8)**2] * a**2 = (F/2)**2, for every
      detuning sigma, classifies each steady state by the eigenvalues of the
      linearised slow flow, and marks the two saddle-node (jump) points.
  (b) It integrates the slow flow itself and compares the predicted envelope
      a(T1) with the envelope of a full numerical simulation of the original
      second-order equation -- the check that the perturbation result is
      quantitatively right.
  (c) It shows the amplitude evolution a(T1) from several initial amplitudes,
      each relaxing onto one of the two stable steady states while the middle
      one repels.
  (d) It draws the slow-flow phase portrait in the Cartesian slow variables
      p = a*cos(gamma), q = a*sin(gamma), where the flow is smooth at a = 0,
      with the two stable foci and the saddle.

Governing equation (plain text):

    x'' + 2*eps*mu*x' + x + eps*x**3 = eps*F*cos(Omega*t),
    Omega = 1 + eps*sigma          (primary resonance with detuning sigma)

Runtime: a few seconds only, dominated by the full-equation simulations (their
length scales with 1/eps, because the slow dynamics live on T1 = eps*t).
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# =====================================================================
# PARAMETERS -- change these to explore the system
# =====================================================================

EPS = 0.05      # small parameter; damping, nonlinearity and forcing are O(eps)
MU = 0.25       # scaled damping coefficient (physical damping = 2*eps*mu)
F = 1.0         # scaled forcing amplitude (physical amplitude = eps*F)

# detuning window for the frequency-response curve
SIGMA_MIN, SIGMA_MAX = -1.0, 3.5

# detuning used for the transient / phase-portrait panels.  It is chosen
# inside the multi-valued window so that two stable steady states coexist.
SIGMA_DEMO = 1.1

# initial amplitudes for the envelope comparison and the a(T1) panel
A0_LIST = (0.05, 1.60, 2.20, 3.00)

# slow time horizon; physical time is T1_END/EPS
T1_END = 40.0

# numerics
RTOL, ATOL = 1e-10, 1e-12
MAX_STEP_FULL = 0.25    # cap on the step of the full simulation, so that the
                        # fast carrier (period ~ 2*pi) is always well resolved
REAL_TOL = 1e-9         # tolerance for accepting a root of the cubic as real


# =====================================================================
# Slow flow (modulation equations)
# =====================================================================

def slow_flow_cartesian(_t1, pq, sigma):
    """Modulation equations in the Cartesian slow variables p, q.

    Polar form (solvability condition 2i*A' + 2i*mu*A + 3|A|^2*A =
    (F/2)*exp(i*sigma*T1) with A = (a/2)*exp(i*psi) and gamma = sigma*T1 - psi):

        da/dT1       = -mu*a + (F/2)*sin(gamma)
        a*dgamma/dT1 =  sigma*a - (3/8)*a**3 + (F/2)*cos(gamma)

    Note on conventions: the slides write both forcing terms with a minus sign,
    which is the same system with gamma replaced by gamma + pi.  The physics --
    and in particular the frequency-response equation below, in which F appears
    squared -- is identical; the signs used here are the ones consistent with
    x ~ a*cos(Omega*t - gamma), so that the phase portrait can be compared
    directly with the full simulation.

    With p = a*cos(gamma) and q = a*sin(gamma) this becomes

        dp/dT1 = -mu*p - (sigma - 3*(p**2+q**2)/8)*q
        dq/dT1 = -mu*q + (sigma - 3*(p**2+q**2)/8)*p + F/2,

    which is polynomial and, unlike the polar form, free of the coordinate
    singularity at a = 0.  The time variable is the slow time T1 = eps*t.
    """
    p, q = pq
    detune = sigma - 0.375 * (p * p + q * q)
    return [-MU * p - detune * q,
            -MU * q + detune * p + 0.5 * F]


def slow_flow_jacobian(p, q, sigma):
    """Jacobian of the Cartesian slow flow, used for the stability test."""
    a2 = p * p + q * q
    detune = sigma - 0.375 * a2
    return np.array([[-MU + 0.75 * p * q, -detune + 0.75 * q * q],
                     [detune - 0.75 * p * p, -MU - 0.75 * p * q]])


def steady_state_amplitudes(sigma):
    """Steady-state amplitudes from the cubic frequency-response equation.

    Setting da/dT1 = dgamma/dT1 = 0 and eliminating gamma via
    sin^2 + cos^2 = 1 gives

        [mu**2 + (sigma - 3*a**2/8)**2] * a**2 = (F/2)**2,

    a cubic in z = a**2:

        (9/64)*z**3 - (3*sigma/4)*z**2 + (mu**2 + sigma**2)*z - F**2/4 = 0.

    Returns the positive real amplitudes in ascending order (one or three).
    """
    coeffs = [9.0 / 64.0, -0.75 * sigma, MU**2 + sigma**2, -0.25 * F**2]
    roots = np.roots(coeffs)
    real_roots = roots[np.abs(roots.imag) < REAL_TOL].real
    return np.sqrt(np.sort(real_roots[real_roots > 0.0]))


def steady_state_pq(amplitude, sigma):
    """Fixed point (p, q) of the slow flow belonging to a steady amplitude.

    From the steady-state equations,
        sin(gamma) = 2*mu*a/F,
        cos(gamma) = (2*a/F)*(3*a**2/8 - sigma),
    hence p = a*cos(gamma) and q = a*sin(gamma).
    """
    p = (2.0 * amplitude**2 / F) * (0.375 * amplitude**2 - sigma)
    q = 2.0 * MU * amplitude**2 / F
    return p, q


def is_stable(amplitude, sigma):
    """True if both eigenvalues of the linearised slow flow have Re < 0.

    Following the slides: perturb the steady state and linearise the
    modulation equations.  The middle branch between the two saddle-nodes has
    one positive real eigenvalue (a saddle) and is never observed.
    """
    p, q = steady_state_pq(amplitude, sigma)
    eig = np.linalg.eigvals(slow_flow_jacobian(p, q, sigma))
    return bool(np.all(eig.real < 0.0))


def frequency_response(sigma_values):
    """Trace the frequency-response curve, split into stable / unstable points.

    Returns (sigma_stable, a_stable, sigma_unstable, a_unstable).
    """
    s_st, a_st, s_un, a_un = [], [], [], []
    for sigma in sigma_values:
        for amp in steady_state_amplitudes(sigma):
            if is_stable(amp, sigma):
                s_st.append(sigma)
                a_st.append(amp)
            else:
                s_un.append(sigma)
                a_un.append(amp)
    return (np.array(s_st), np.array(a_st),
            np.array(s_un), np.array(a_un))


def saddle_node_points(sigma_values):
    """Locate the folds (saddle-nodes) of the frequency-response curve.

    The saddle-nodes are the detunings at which the cubic acquires a double
    root, i.e. where the number of positive real steady states switches
    between one and three.  They are bracketed on the sigma grid; on the side
    of the bracket that still has three roots the two merging branches are the
    closest pair, so their mean is the fold amplitude.

    Returns a list of (sigma_fold, a_fold) pairs.
    """
    counts = np.array([len(steady_state_amplitudes(s)) for s in sigma_values])
    folds = []
    for i in np.nonzero(np.diff(counts) != 0)[0]:
        j = i if counts[i] == 3 else i + 1
        amps = steady_state_amplitudes(sigma_values[j])
        k = int(np.argmin(np.diff(amps)))       # the pair about to merge
        folds.append((sigma_values[j], 0.5 * (amps[k] + amps[k + 1])))
    return folds


# =====================================================================
# Full equation of motion (reference solution)
# =====================================================================

def full_rhs(t, y, sigma):
    """Right-hand side of the original forced Duffing oscillator.

        x'' + 2*eps*mu*x' + x + eps*x**3 = eps*F*cos(Omega*t),
        Omega = 1 + eps*sigma.
    """
    x, v = y
    omega_exc = 1.0 + EPS * sigma
    accel = (-2.0 * EPS * MU * v - x - EPS * x**3
             + EPS * F * np.cos(omega_exc * t))
    return [v, accel]


def simulate_full(a0, sigma):
    """Simulate the full equation from x(0) = a0, x'(0) = 0.

    Because the carrier frequency is close to 1, the analytic-signal amplitude
    of the response is well approximated by sqrt(x**2 + x'**2); that quantity
    is exactly the envelope a(T1) predicted by the multiple-scales solution
    x ~ a*cos(T0 + beta).

    Returns (T1, envelope, t, x) with the slow time T1 = eps*t.
    """
    t_end = T1_END / EPS
    t_eval = np.linspace(0.0, t_end, 20000)
    sol = solve_ivp(full_rhs, [0.0, t_end], [a0, 0.0], t_eval=t_eval,
                    args=(sigma,), method="DOP853",
                    rtol=RTOL, atol=ATOL, max_step=MAX_STEP_FULL)
    envelope = np.hypot(sol.y[0], sol.y[1])
    return EPS * sol.t, envelope, sol.t, sol.y[0]


def simulate_slow(a0, sigma):
    """Integrate the modulation equations from the matching initial condition.

    x(0) = a0, x'(0) = 0 corresponds to amplitude a0 with relative phase
    gamma = 0, i.e. p = a0, q = 0.

    Returns (T1, a, p, q).
    """
    t1 = np.linspace(0.0, T1_END, 4000)
    sol = solve_ivp(slow_flow_cartesian, [0.0, T1_END], [a0, 0.0],
                    t_eval=t1, args=(sigma,), method="DOP853",
                    rtol=RTOL, atol=ATOL)
    p, q = sol.y
    return sol.t, np.hypot(p, q), p, q


# =====================================================================
# Main
# =====================================================================

def main():
    # ---------- (a) frequency response and stability ------------------
    sigma_grid = np.linspace(SIGMA_MIN, SIGMA_MAX, 3000)
    s_st, a_st, s_un, a_un = frequency_response(sigma_grid)
    folds = saddle_node_points(sigma_grid)

    print("Modulation equations: eps=%.3f, mu=%.3f, F=%.3f" % (EPS, MU, F))
    if len(folds) >= 2:
        for s_fold, a_fold in folds:
            print("  saddle-node: sigma = %.3f, a = %.3f" % (s_fold, a_fold))
        print("  -> three steady states coexist for %.3f < sigma < %.3f"
              % (folds[0][0], folds[-1][0]))
    print("  steady-state amplitudes at sigma = %.2f:" % SIGMA_DEMO)
    demo_amps = steady_state_amplitudes(SIGMA_DEMO)
    for amp in demo_amps:
        print("      a_s = %.4f   %s"
              % (amp, "stable" if is_stable(amp, SIGMA_DEMO) else "unstable"))

    # ---------- (b)+(c) transients, MMS vs. full simulation -----------
    print("Comparing slow-flow envelope with the full simulation ...")
    runs = []
    for a0 in A0_LIST:
        t1_full, env_full, t_full, x_full = simulate_full(a0, SIGMA_DEMO)
        t1_slow, a_slow, _, _ = simulate_slow(a0, SIGMA_DEMO)
        runs.append(dict(a0=a0, t1_full=t1_full, env_full=env_full,
                         t_full=t_full, x_full=x_full,
                         t1_slow=t1_slow, a_slow=a_slow))
        print("    a(0)=%.2f -> full %.4f, MMS %.4f  (relative error %.2f %%)"
              % (a0, env_full[-1], a_slow[-1],
                 100.0 * abs(env_full[-1] - a_slow[-1]) / a_slow[-1]))

    # ---------- figure -----------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.0))
    fig.subplots_adjust(hspace=0.34, wspace=0.26, top=0.89, bottom=0.08,
                        left=0.07, right=0.97)
    ax_fr, ax_env, ax_amp, ax_pp = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    # (a) frequency response with stability
    a_bb = np.linspace(0.0, 2.6, 200)
    ax_fr.plot(0.375 * a_bb**2, a_bb, ":", color="0.6", lw=1.2,
               label=r"backbone $\sigma=\frac{3}{8}a^2$")
    ax_fr.plot(s_st, a_st, ".", ms=2.2, color="tab:blue")
    ax_fr.plot(s_un, a_un, ".", ms=2.2, color="tab:red")
    ax_fr.plot([], [], "-", lw=2.0, color="tab:blue", label="stable")
    ax_fr.plot([], [], "-", lw=2.0, color="tab:red", label="unstable")
    for s_fold, a_fold in folds:
        ax_fr.plot(s_fold, a_fold, "o", ms=7, mfc="none", mew=1.4,
                   color="k", zorder=6)
    ax_fr.plot([], [], "o", ms=7, mfc="none", mew=1.4, color="k",
               label="saddle-node")
    ax_fr.axvline(SIGMA_DEMO, color="0.4", lw=0.8, ls="--")
    ax_fr.set_xlim(SIGMA_MIN, SIGMA_MAX)
    ax_fr.set_ylim(0.0, 2.7)
    ax_fr.set_xlabel(r"detuning $\sigma$   ($\Omega=1+\varepsilon\sigma$)")
    ax_fr.set_ylabel(r"steady-state amplitude $a_s$")
    ax_fr.set_title("(a) Frequency response from the modulation equations")
    ax_fr.legend(loc="upper left", fontsize=8)
    ax_fr.grid(alpha=0.25)

    # (b) envelope: full simulation vs. slow flow, for one initial amplitude
    demo = runs[0]
    ax_env.plot(demo["t1_full"], demo["x_full"], lw=0.4, color="0.75",
                label=r"full simulation $x(t)$")
    ax_env.plot(demo["t1_full"], demo["env_full"], lw=1.3, color="tab:blue",
                label=r"envelope of $x(t)$")
    ax_env.plot(demo["t1_slow"], demo["a_slow"], lw=1.6, ls="--",
                color="tab:orange", label=r"MMS slow flow $a(T_1)$")
    ax_env.set_xlabel(r"slow time $T_1=\varepsilon t$")
    ax_env.set_ylabel(r"$x$,  $a$")
    ax_env.set_title(r"(b) Ring-up at $\sigma=%.2f$, $a(0)=%.2f$"
                     % (SIGMA_DEMO, demo["a0"]))
    ax_env.legend(loc="lower right", fontsize=8)
    ax_env.grid(alpha=0.25)

    # (c) amplitude evolution for several initial amplitudes
    for k, run in enumerate(runs):
        color = "C%d" % k
        ax_amp.plot(run["t1_full"], run["env_full"], lw=0.9, color=color,
                    alpha=0.55)
        ax_amp.plot(run["t1_slow"], run["a_slow"], lw=1.6, ls="--",
                    color=color, label=r"$a(0)=%.2f$" % run["a0"])
    for amp in demo_amps:
        stable = is_stable(amp, SIGMA_DEMO)
        ax_amp.axhline(amp, color="k", lw=0.8,
                       ls="-" if stable else ":", alpha=0.6)
    ax_amp.set_xlabel(r"slow time $T_1=\varepsilon t$")
    ax_amp.set_ylabel(r"amplitude $a$")
    ax_amp.set_title("(c) Relaxation onto the stable steady states\n"
                     "(thin: full simulation, dashed: slow flow; "
                     "solid/dotted lines: stable/unstable $a_s$)",
                     fontsize=10)
    ax_amp.legend(loc="center right", fontsize=8)
    ax_amp.grid(alpha=0.25)

    # (d) slow-flow phase portrait in p, q
    lim = 2.9
    grid = np.linspace(-lim, lim, 26)
    pg, qg = np.meshgrid(grid, grid)
    dp, dq = slow_flow_cartesian(0.0, (pg, qg), SIGMA_DEMO)
    speed = np.hypot(dp, dq)
    ax_pp.streamplot(pg, qg, dp, dq, color=speed, cmap="Blues",
                     density=1.1, linewidth=0.7, arrowsize=0.8)
    for run in runs:
        _, _, p_tr, q_tr = simulate_slow(run["a0"], SIGMA_DEMO)
        ax_pp.plot(p_tr, q_tr, lw=1.2, color="0.25", alpha=0.8)
    for amp in demo_amps:
        p_fp, q_fp = steady_state_pq(amp, SIGMA_DEMO)
        stable = is_stable(amp, SIGMA_DEMO)
        ax_pp.plot(p_fp, q_fp, "o", ms=7,
                   color="tab:blue" if stable else "tab:red",
                   mec="k", mew=0.6, zorder=5)
    ax_pp.plot([], [], "o", ms=7, color="tab:blue", mec="k", mew=0.6,
               label="stable focus")
    ax_pp.plot([], [], "o", ms=7, color="tab:red", mec="k", mew=0.6,
               label="saddle")
    ax_pp.set_xlim(-lim, lim)
    ax_pp.set_ylim(-lim, lim)
    ax_pp.set_aspect("equal")
    ax_pp.set_xlabel(r"$p=a\cos\gamma$")
    ax_pp.set_ylabel(r"$q=a\sin\gamma$")
    ax_pp.set_title(r"(d) Slow-flow phase portrait at $\sigma=%.2f$"
                    % SIGMA_DEMO)
    ax_pp.legend(loc="upper left", fontsize=8)

    fig.suptitle(r"Method of multiple scales: "
                 r"$\ddot{x}+2\varepsilon\mu\dot{x}+x+\varepsilon x^3"
                 r"=\varepsilon F\cos\Omega t$,  "
                 r"$\varepsilon=%.2f$, $\mu=%.2f$, $F=%.2f$"
                 % (EPS, MU, F), fontsize=12)
    plt.show()


if __name__ == "__main__":
    main()
