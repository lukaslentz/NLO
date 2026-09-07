"""
Chapter 2.1 -- Phase Space Analysis

Phase-plane analysis of two planar autonomous systems, following the slide set
"2.1 Phase Space Analysis".  For each system the script

  * draws the direction field of  x_dot = f(x)  (normalised arrows, so that only
    the direction and not the local speed is visible),
  * integrates a handful of trajectories with an adaptive Runge-Kutta scheme and
    overlays them on the direction field,
  * locates all equilibria f(x*) = 0 numerically, evaluates the Jacobian
    J = df/dx at each of them and classifies the equilibrium from the
    eigenvalues of J (centre / stable-unstable focus / node / saddle).  The
    classification, the eigenvalues and the Poincare index are printed to the
    console.

System A -- van der Pol oscillator (the slide example, mu = 1.5):

    x'' - mu (1 - x^2) x' + x = 0        i.e.   x1' = x2
                                                x2' = mu (1 - x1^2) x2 - x1

  The single equilibrium at the origin is an unstable focus; by the
  Poincare-Bendixson theorem the bounded, equilibrium-free annulus around it
  must contain a limit cycle, and indeed all trajectories -- started inside and
  outside -- converge onto one isolated closed orbit.

System B -- unforced bistable Duffing oscillator (the nullcline slide,
damping 0.9):

    x'' + 0.9 x' - x + x^3 = 0           i.e.   x1' = x2
                                                x2' = x1 - x1^3 - 0.9 x2

  Three equilibria: a saddle at the origin and two stable foci at x1 = +-1.
  The stable manifold of the saddle is the separatrix that divides the plane
  into the two basins of attraction; it is approximated here by integrating
  backwards in time from the saddle along its stable eigenvector.

The resulting figure has two panels, one per system: direction field in grey,
trajectories in colour, equilibria as markers coded by type, and for system B
the two branches of the stable manifold (separatrix) as a dashed black curve.

Runtime: a few seconds.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# ---------------------------------------------------------------------------
# PARAMETERS -- change these to explore
# ---------------------------------------------------------------------------
MU = 1.5                 # van der Pol nonlinearity parameter (slide value)
ZETA_DUFFING = 0.9       # linear damping of the bistable Duffing system

VDP_XLIM = (-3.0, 3.0)   # phase-plane window, system A
VDP_VLIM = (-4.0, 4.0)
VDP_START = [0.3, 1.0, 2.0, 2.8]   # initial x with x_dot = 0 (slide values)
VDP_TEND = 40.0          # integration time (long enough to reach the cycle)

DUF_XLIM = (-2.0, 2.0)   # phase-plane window, system B
DUF_VLIM = (-2.0, 2.0)
DUF_TEND = 25.0

N_ARROWS = 22            # arrows per axis in the direction field
RTOL, ATOL = 1e-9, 1e-11 # ODE tolerances
EPS_JAC = 1e-6           # step of the central difference for the Jacobian
# ---------------------------------------------------------------------------


def van_der_pol(t, y, mu=MU):
    """Van der Pol oscillator x'' - mu (1 - x^2) x' + x = 0 as a first-order system.

    For |x| < 1 the damping term is negative (energy is pumped in), for
    |x| > 1 it is positive (energy is dissipated).  This is the mechanism that
    creates the isolated limit cycle.
    """
    x, v = y
    return np.array([v, mu * (1.0 - x**2) * v - x])


def duffing_bistable(t, y, zeta=ZETA_DUFFING):
    """Unforced bistable Duffing oscillator x'' + zeta x' - x + x^3 = 0.

    The potential V(x) = -x^2/2 + x^4/4 is a double well with minima at
    x = +-1 and a local maximum at x = 0, hence two stable equilibria and one
    saddle.
    """
    x, v = y
    return np.array([v, x - x**3 - zeta * v])


def jacobian(f, x_star):
    """Numerical Jacobian of the autonomous field f at x_star (central differences).

    A central difference is used because it is second-order accurate, so the
    eigenvalues are correct to roughly 1e-10 for the smooth polynomial fields
    considered here -- accurate enough for a stability classification.
    """
    n = len(x_star)
    J = np.zeros((n, n))
    for j in range(n):
        e = np.zeros(n)
        e[j] = EPS_JAC
        J[:, j] = (f(0.0, x_star + e) - f(0.0, x_star - e)) / (2.0 * EPS_JAC)
    return J


def classify(eigenvalues, tol=1e-8):
    """Classify a planar equilibrium from the two eigenvalues of its Jacobian.

    Returns (label, Poincare index).  The index is -1 for a saddle and +1 for
    every other non-degenerate type; a closed orbit must enclose equilibria
    whose indices sum to +1.
    """
    lam1, lam2 = eigenvalues
    re = np.real(eigenvalues)
    im = np.imag(eigenvalues)

    if np.all(np.abs(im) > tol):                 # complex conjugate pair
        if np.all(np.abs(re) < tol):
            return "centre (linearisation inconclusive)", +1
        return ("stable focus" if re[0] < 0 else "unstable focus"), +1

    # real eigenvalues
    if re[0] * re[1] < -tol:
        return "saddle point", -1
    if np.all(re < -tol):
        return "stable node", +1
    if np.all(re > tol):
        return "unstable node", +1
    return "non-hyperbolic (linearisation inconclusive)", 0


def find_equilibria(f, xlim, vlim, n_seed=15, tol=1e-7):
    """Locate all equilibria in a window by Newton iterations from a seed grid.

    fsolve is started from a coarse grid of initial guesses; the roots found are
    deduplicated and those outside the window (or spurious, i.e. |f| not small)
    are discarded.
    """
    seeds = [np.array([x, v])
             for x in np.linspace(*xlim, n_seed)
             for v in np.linspace(*vlim, n_seed)]
    roots = []
    for s in seeds:
        r, _, ier, _ = fsolve(lambda z: f(0.0, z), s, full_output=True)
        if ier != 1 or np.linalg.norm(f(0.0, r)) > tol:
            continue
        if not (xlim[0] <= r[0] <= xlim[1] and vlim[0] <= r[1] <= vlim[1]):
            continue
        if all(np.linalg.norm(r - q) > 1e-5 for q in roots):
            roots.append(r)
    return sorted(roots, key=lambda r: (r[0], r[1]))


def report_equilibria(name, f, xlim, vlim):
    """Find, classify and print all equilibria of f; return a list of records."""
    print(f"\n{name}")
    print("-" * len(name))
    records = []
    for x_star in find_equilibria(f, xlim, vlim):
        J = jacobian(f, x_star)
        lam = np.linalg.eigvals(J)
        # sort by real part so the printed order is reproducible
        lam = lam[np.argsort(np.real(lam))]
        label, index = classify(lam)
        records.append({"x": x_star, "eig": lam, "label": label, "index": index})
        print(f"  x* = ({x_star[0]:+.4f}, {x_star[1]:+.4f})")
        print(f"     trace J = {np.trace(J):+.4f},  det J = {np.linalg.det(J):+.4f}")
        print(f"     lambda  = {lam[0]:+.4f}, {lam[1]:+.4f}")
        print(f"     type    = {label}   (Poincare index {index:+d})")
    print(f"  sum of indices = {sum(r['index'] for r in records):+d}")
    return records


def direction_field(ax, f, xlim, vlim, n=N_ARROWS):
    """Draw normalised arrows of the vector field f on the axes ax.

    Only the direction matters for the qualitative phase portrait, so every
    arrow is scaled to unit length; the small epsilon avoids a division by zero
    exactly at an equilibrium.
    """
    X, V = np.meshgrid(np.linspace(*xlim, n), np.linspace(*vlim, n))
    DX, DV = f(0.0, [X, V])
    norm = np.hypot(DX, DV) + 1e-12
    ax.quiver(X, V, DX / norm, DV / norm, color="0.6",
              alpha=0.7, pivot="mid", width=0.0035)


def integrate(f, y0, t_end, n_out=4000, backward=False):
    """Integrate f from y0 over [0, t_end] (or backwards) with dense output."""
    t_span = (0.0, -t_end) if backward else (0.0, t_end)
    t_eval = np.linspace(*t_span, n_out)
    sol = solve_ivp(f, t_span, y0, t_eval=t_eval, rtol=RTOL, atol=ATOL,
                    method="RK45")
    return sol.y


MARKERS = {"saddle point": ("s", "black"),
           "stable focus": ("o", "tab:green"),
           "unstable focus": ("o", "tab:red"),
           "stable node": ("o", "tab:green"),
           "unstable node": ("o", "tab:red")}


def plot_equilibria(ax, records):
    """Mark equilibria with a type-dependent symbol; unstable points are hollow."""
    seen = set()
    for rec in records:
        marker, colour = MARKERS.get(rec["label"], ("*", "tab:purple"))
        filled = rec["label"].startswith("stable")
        label = rec["label"].capitalize() if rec["label"] not in seen else None
        seen.add(rec["label"])
        ax.plot(rec["x"][0], rec["x"][1], marker, ms=9, color=colour,
                mfc=colour if filled else "white", mew=1.6, zorder=5,
                label=label)


def stable_manifold(f, records, length=1e-3, t_end=12.0):
    """Approximate the stable manifold (separatrix) of a saddle.

    Start two points a tiny distance along the +- stable eigenvector and
    integrate BACKWARDS in time: backward integration turns the attracting
    stable direction into a repelling one, so the trajectory traces out W^s.
    """
    branches = []
    for rec in records:
        if rec["label"] != "saddle point":
            continue
        J = jacobian(f, rec["x"])
        lam, vec = np.linalg.eig(J)
        k = int(np.argmin(np.real(lam)))          # the negative eigenvalue
        w = np.real(vec[:, k])
        w /= np.linalg.norm(w)
        for sign in (+1.0, -1.0):
            branches.append(integrate(f, rec["x"] + sign * length * w,
                                      t_end, backward=True))
    return branches


def panel_van_der_pol(ax):
    """Left panel: direction field, trajectories and limit cycle of van der Pol."""
    direction_field(ax, van_der_pol, VDP_XLIM, VDP_VLIM)

    # Trajectories started inside and outside the cycle all wind onto it.
    for i, x0 in enumerate(VDP_START):
        y = integrate(van_der_pol, [x0, 0.0], VDP_TEND)
        ax.plot(y[0], y[1], lw=0.9, color="steelblue",
                label="Trajectories" if i == 0 else None)

    # The last quarter of a long run is (numerically) on the limit cycle.
    y = integrate(van_der_pol, [0.3, 0.0], VDP_TEND, n_out=8000)
    tail = y[:, -2000:]
    ax.plot(tail[0], tail[1], lw=2.2, color="crimson", label="Limit cycle")

    records = report_equilibria(
        f"System A: van der Pol oscillator (mu = {MU})",
        van_der_pol, VDP_XLIM, VDP_VLIM)
    plot_equilibria(ax, records)

    amp = np.max(np.abs(tail[0]))
    print(f"  limit-cycle amplitude max|x| = {amp:.4f}")

    ax.set_xlim(VDP_XLIM)
    ax.set_ylim(VDP_VLIM)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$\dot{x}$")
    ax.set_title(rf"van der Pol, $\mu = {MU}$: unstable focus + limit cycle")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    return amp


def panel_duffing(ax):
    """Right panel: direction field, trajectories, separatrix of bistable Duffing."""
    direction_field(ax, duffing_bistable, DUF_XLIM, DUF_VLIM)

    records = report_equilibria(
        f"System B: bistable Duffing oscillator (damping {ZETA_DUFFING})",
        duffing_bistable, DUF_XLIM, DUF_VLIM)

    # A ring of initial conditions shows how the two basins interleave.
    first = True
    for angle in np.linspace(0.0, 2.0 * np.pi, 17)[:-1]:
        y0 = 1.85 * np.array([np.cos(angle), np.sin(angle)])
        y = integrate(duffing_bistable, y0, DUF_TEND)
        ax.plot(y[0], y[1], lw=0.8, color="steelblue",
                label="Trajectories" if first else None)
        first = False

    for i, branch in enumerate(stable_manifold(duffing_bistable, records)):
        ax.plot(branch[0], branch[1], "k--", lw=1.6,
                label=r"Stable manifold $W^s$ (separatrix)" if i == 0 else None)

    plot_equilibria(ax, records)

    ax.set_xlim(DUF_XLIM)
    ax.set_ylim(DUF_VLIM)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$\dot{x}$")
    ax.set_title(r"Bistable Duffing: saddle at $0$, foci at $\pm 1$")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)


def main():
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6))
    panel_van_der_pol(axes[0])
    panel_duffing(axes[1])
    fig.suptitle("2.1 Phase space: direction fields, trajectories and "
                 "classification of equilibria")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
