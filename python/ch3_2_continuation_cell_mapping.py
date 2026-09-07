"""
3.2 Continuation & Cell Mapping
===============================

Two global analysis techniques for the forced Duffing oscillator

    x'' + 2*delta*x' + x + eps*x**3 = fhat*cos(Omega*t)

are demonstrated side by side.

(1) Pseudo-arc-length continuation of the frequency response curve.
    Harmonic balance with the ansatz x(t) = a*cos(Omega*t + phi) reduces the
    steady-state problem to one scalar equation for the amplitude a,

        F(a, Omega) = [(1 - Omega**2) + (3/4)*eps*a**2]**2 * a**2
                      + (2*delta*Omega)**2 * a**2 - fhat**2 = 0 .

    Solving this for a at a prescribed Omega fails near the folds, where the
    branch has a vertical tangent in the (Omega, a) plane.  Parameterising the
    branch by its arc length s instead and solving the extended system
    [F = 0 ; tangent-hyperplane condition = 0] with a predictor-corrector
    Newton scheme walks smoothly around both folds and therefore also traces
    the unstable middle branch.  The fold points themselves are then pinned
    down by solving F = 0 and dF/da = 0 simultaneously.

(2) Simple Cell Mapping (SCM, Hsu 1980) at a frequency between the two folds,
    where three solutions coexist.  The plane (x, x') in [-3, 3]^2 is cut into
    cells; the ODE is integrated from every cell centre over whole excitation
    periods T = 2*pi/Omega, and the cell that contains the landing point
    becomes the image cell.  Following the resulting deterministic cell-to-cell
    map until a cell repeats yields the periodic groups (the discrete
    attractors) and, by bookkeeping, the basin of attraction of each of them.

The figure has two panels: left the frequency response with stable branches
solid, the unstable branch dashed and the two fold points marked; right the
SCM basins of attraction with the attractor cells highlighted.  The two SCM
attractors reproduce the amplitudes of the two stable harmonic balance
branches, which is the numerical check that both parts agree.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve

# ----------------------------------------------------------------------
# PARAMETERS  (values taken from the slides)
# ----------------------------------------------------------------------
DELTA = 0.05        # damping ratio delta  (equation has 2*delta*x')
EPS = 0.1           # cubic stiffness coefficient, > 0 -> hardening spring
FHAT = 0.3          # forcing amplitude

# Excitation frequency for the cell mapping panel.  Cell mapping only shows
# something interesting inside the bistable window, which for delta = 0.05,
# eps = 0.1, fhat = 0.3 is bounded by the two folds at Omega = 1.152 and
# Omega = 1.212 computed by the continuation below.  Rather than hard-coding a
# frequency, OMEGA_CM = None picks the midpoint of that computed interval, so
# the choice stays correct if the parameters above are changed.  Put a number
# here to override (values between the two folds make sense; outside them only
# one solution exists and every cell maps to the same attractor).
OMEGA_CM = None

# --- continuation settings ---
OMEGA_START = 0.40  # start the branch here (single solution, small amplitude)
OMEGA_END = 2.00    # stop once the branch has passed this frequency again
DS = 0.02           # pseudo-arc-length step length
DS_MIN, DS_MAX = 2e-4, 0.05
NEWTON_TOL = 1e-11
NEWTON_MAX = 12
MAX_STEPS = 4000

# --- simple cell mapping settings ---
X_RANGE = (-3.0, 3.0)   # cell region in x
V_RANGE = (-3.0, 3.0)   # cell region in x'
N_CELLS = 60            # cells per direction -> N_CELLS**2 cells in total
N_PERIODS = 10          # excitation periods integrated per cell (see below)
# Runtime scales with N_CELLS**2 * N_PERIODS.  60 x 60 = 3600 cells takes a few
# seconds because all cell centres are integrated as ONE large ODE system.
# Refining to N_CELLS = 200 (40000 cells) sharpens the basin boundary but costs
# roughly 30x more time and memory - well outside the 60 s budget here.
#
# Why N_PERIODS > 1: the textbook description of SCM integrates exactly one
# period per cell.  With a coarse grid and light damping (delta = 0.05) the
# contraction during a single period is smaller than one cell width, so dozens
# of cells around an attractor map onto themselves and SCM reports dozens of
# spurious one-cell periodic groups for one and the same attractor.  Using the
# k-th iterate of the stroboscopic map (k = N_PERIODS) is the standard cure:
# the trajectories contract by exp(-delta*k*T) before the image cell is taken,
# so only the cells that really contain an attractor stay self-mapping.  The
# map is still a deterministic cell-to-cell map and the basins are unchanged.
# Set N_PERIODS = 1 to see the fragmented picture of the literal algorithm.

RTOL, ATOL = 1e-9, 1e-11


# ----------------------------------------------------------------------
# HARMONIC BALANCE RESIDUAL AND ITS DERIVATIVES
# ----------------------------------------------------------------------
def hb_residual(a, om):
    """Harmonic balance residual F(a, Omega) of the Duffing oscillator."""
    r = (1.0 - om**2) + 0.75 * EPS * a**2      # effective detuning term
    return (r**2 + (2.0 * DELTA * om) ** 2) * a**2 - FHAT**2


def hb_dF_da(a, om):
    """Partial derivative dF/da (singular exactly at the fold points)."""
    r = (1.0 - om**2) + 0.75 * EPS * a**2
    return 2.0 * a * r**2 + 3.0 * EPS * a**3 * r + 8.0 * DELTA**2 * om**2 * a


def hb_dF_dom(a, om):
    """Partial derivative dF/dOmega."""
    r = (1.0 - om**2) + 0.75 * EPS * a**2
    return -4.0 * om * a**2 * r + 8.0 * DELTA**2 * om * a**2


def hb_solutions(om):
    """
    All harmonic balance amplitudes at a given Omega, with their stability.

    In the squared amplitude u = a**2 the residual is the cubic

        b**2 u**3 + 2*b*c*u**2 + (c**2 + 4*delta**2*Omega**2) u - fhat**2 = 0,
        c = 1 - Omega**2,  b = (3/4)*eps,

    so the (one or three) coexisting responses come straight from its real
    positive roots - no bracketing guesswork needed.
    """
    c, b = 1.0 - om**2, 0.75 * EPS
    roots = np.roots([b**2, 2.0 * b * c, c**2 + 4.0 * DELTA**2 * om**2,
                      -FHAT**2])
    amps = np.sqrt(sorted(r.real for r in roots
                          if abs(r.imag) < 1e-10 and r.real > 0.0))
    return [(a, is_stable(a, om)) for a in amps]


def is_stable(a, om):
    """
    Stability of a harmonic balance solution.

    Averaging the Duffing oscillator gives a slow-flow system whose Jacobian
    has trace -2*delta < 0 and determinant proportional to dG/d(a**2), with
    G(a**2) the left-hand side of the amplitude equation.  Hence the solution
    is stable exactly where dG/d(a**2) > 0, i.e. where dF/da has the same sign
    as a.  This flips sign at the folds, which is why the middle branch between
    the two folds is the unstable one.
    """
    return hb_dF_da(a, om) / max(a, 1e-14) > 0.0


# ----------------------------------------------------------------------
# PSEUDO-ARC-LENGTH CONTINUATION
# ----------------------------------------------------------------------
def tangent(v, prev_t):
    """
    Unit tangent of the solution branch at v = (a, Omega).

    The branch is the null space of the 1x2 Jacobian [dF/da, dF/dOmega], so the
    tangent is simply (dF/dOmega, -dF/da) normalised.  Its sign is fixed by
    requiring a positive projection on the previous tangent, which keeps the
    continuation going forward around a fold instead of turning back.
    """
    t = np.array([hb_dF_dom(v[0], v[1]), -hb_dF_da(v[0], v[1])])
    t /= np.linalg.norm(t)
    if prev_t is not None and np.dot(t, prev_t) < 0.0:
        t = -t
    return t


def corrector(v_pred, t, ds_ref):
    """
    Newton corrector for the extended system

        F(a, Omega)                      = 0
        (v - v_pred) . t                 = 0   (pseudo-arc-length constraint)

    The second equation restricts the correction to the hyperplane orthogonal
    to the predictor direction.  The 2x2 Jacobian stays regular at a fold,
    where the plain 1x1 problem dF/da = 0 breaks down.
    """
    v = v_pred.copy()
    for _ in range(NEWTON_MAX):
        res = np.array([hb_residual(v[0], v[1]), np.dot(v - v_pred, t)])
        if np.linalg.norm(res) < NEWTON_TOL:
            return v, True
        jac = np.array([[hb_dF_da(v[0], v[1]), hb_dF_dom(v[0], v[1])],
                        [t[0], t[1]]])
        try:
            v = v - np.linalg.solve(jac, res)
        except np.linalg.LinAlgError:
            return v, False
        if v[0] <= 0.0 or not np.all(np.isfinite(v)):
            return v, False
    return v, np.linalg.norm(
        [hb_residual(v[0], v[1]), np.dot(v - v_pred, t)]) < 1e-8 * ds_ref


def continue_branch():
    """
    Trace the whole frequency response curve by pseudo-arc-length continuation.

    Returns the array of branch points (a, Omega) and the array of tangents,
    ordered along the branch.
    """
    # Starting point: at Omega_start the response is single valued, so an
    # ordinary bracketed root find is enough to get on the branch.
    a0 = brentq(hb_residual, 1e-8, 10.0, args=(OMEGA_START,), xtol=1e-14)
    v = np.array([a0, OMEGA_START])

    t = tangent(v, None)
    if t[1] < 0.0:                       # start by increasing Omega
        t = -t

    pts, tangs = [v.copy()], [t.copy()]
    ds = DS
    passed_end = False

    for _ in range(MAX_STEPS):
        v_pred = v + ds * t              # predictor: step along the tangent
        v_new, ok = corrector(v_pred, t, ds)

        if not ok:                       # step size adaptation: retry smaller
            ds = max(ds * 0.5, DS_MIN)
            if ds <= DS_MIN:
                break
            continue

        t_new = tangent(v_new, t)
        v, t = v_new, t_new
        pts.append(v.copy())
        tangs.append(t.copy())
        ds = min(ds * 1.15, DS_MAX)      # converged fast -> grow the step

        # Stop after the branch has come back up beyond OMEGA_END.
        if v[1] > OMEGA_END and t[1] > 0.0:
            passed_end = True
        if passed_end and v[1] > OMEGA_END:
            break

    return np.array(pts), np.array(tangs)


def locate_folds(pts, tangs):
    """
    Find the fold (saddle-node) points along the continued branch.

    A fold is signalled by dOmega/ds = 0, i.e. a sign change of the Omega
    component of the tangent between two consecutive branch points.  The
    bracket is then refined by solving the two-equation system
    F = 0, dF/da = 0 with a Newton solver, which is the defining condition of
    a fold (a bisection on dOmega/ds would do as well but converges slower).
    """
    folds = []
    sign = np.sign(tangs[:, 1])
    idx = np.where(sign[:-1] * sign[1:] < 0.0)[0]
    for i in idx:
        guess = 0.5 * (pts[i] + pts[i + 1])
        sol, _, flag, _ = fsolve(
            lambda v: [hb_residual(v[0], v[1]), hb_dF_da(v[0], v[1])],
            guess, full_output=True)
        if flag == 1:
            folds.append(sol)
    return np.array(folds)


# ----------------------------------------------------------------------
# SIMPLE CELL MAPPING
# ----------------------------------------------------------------------
def duffing_many(t, z, om):
    """
    Right-hand side for many Duffing oscillators at once.

    The state vector holds all positions first and all velocities second, so
    that every cell centre of the grid is advanced within a single call to
    solve_ivp.  That is far cheaper than one solve_ivp call per cell.
    """
    n = z.size // 2
    x, v = z[:n], z[n:]
    return np.concatenate((v, -2.0 * DELTA * v - x - EPS * x**3
                           + FHAT * np.cos(om * t)))


def cell_map_images(om):
    """
    Build the deterministic cell-to-cell map P over the rectangular grid.

    Every cell is integrated from its centre for N_PERIODS excitation periods
    of length T = 2*pi/Omega.  Cells whose image leaves the grid are sent to
    the sink cell with index N_CELLS**2 (escaped trajectory).

    Returns (image, xc, vc) with image[k] the index of the image cell of cell k
    and xc, vc the cell centre coordinates in each direction.
    """
    dx = (X_RANGE[1] - X_RANGE[0]) / N_CELLS
    dv = (V_RANGE[1] - V_RANGE[0]) / N_CELLS
    xc = X_RANGE[0] + (np.arange(N_CELLS) + 0.5) * dx
    vc = V_RANGE[0] + (np.arange(N_CELLS) + 0.5) * dv
    XC, VC = np.meshgrid(xc, vc, indexing="ij")     # index (i, j) = (x, v)

    z0 = np.concatenate((XC.ravel(), VC.ravel()))
    period = N_PERIODS * 2.0 * np.pi / om
    sol = solve_ivp(duffing_many, (0.0, period), z0, args=(om,),
                    method="RK45", rtol=1e-7, atol=1e-9, dense_output=False)
    n = z0.size // 2
    x_end, v_end = sol.y[:n, -1], sol.y[n:, -1]

    # Landing point -> cell index; anything outside the grid goes to the sink.
    i = np.floor((x_end - X_RANGE[0]) / dx).astype(int)
    j = np.floor((v_end - V_RANGE[0]) / dv).astype(int)
    inside = (i >= 0) & (i < N_CELLS) & (j >= 0) & (j < N_CELLS)
    sink = N_CELLS * N_CELLS
    image = np.full(sink + 1, sink, dtype=int)   # last entry: sink -> sink
    image[:sink][inside] = i[inside] * N_CELLS + j[inside]
    return image, xc, vc


def scm_groups(image):
    """
    Hsu's unravelling algorithm for the simple cell map.

    Each cell is followed along its images until the chain hits either a cell
    that has already been processed (then the whole chain inherits that group)
    or a cell of the chain itself (then a new periodic group has been found and
    the closed part of the chain is that group).

    Returns
        group    : group number of every cell, sink included (group 0 = sink)
        steps    : number of iterations needed to reach the group
        periodic : boolean flag marking the cells that form a periodic group
    """
    nc = image.size                      # includes the sink cell
    sink = nc - 1
    group = np.zeros(nc, dtype=int)      # 0 = not yet processed
    steps = np.zeros(nc, dtype=int)
    periodic = np.zeros(nc, dtype=bool)
    group[sink] = 0                      # the sink is its own (trivial) group
    processed = np.zeros(nc, dtype=bool)
    processed[sink] = True
    n_groups = 0

    for start in range(nc - 1):
        if processed[start]:
            continue
        chain, seen = [], {}
        cell = start
        while True:
            if processed[cell]:                    # ran into known territory
                g, base = group[cell], steps[cell]
                for k, c in enumerate(chain):
                    group[c] = g
                    steps[c] = base + len(chain) - k
                    processed[c] = True
                break
            if cell in seen:                       # closed a new periodic group
                n_groups += 1
                first = seen[cell]
                for k, c in enumerate(chain):
                    group[c] = n_groups
                    processed[c] = True
                    if k >= first:                 # part of the closed cycle
                        periodic[c] = True
                        steps[c] = 0
                    else:                          # transient leading into it
                        steps[c] = first - k
                break
            seen[cell] = len(chain)
            chain.append(cell)
            cell = image[cell]

    return group, steps, periodic


# ----------------------------------------------------------------------
# PLOTTING
# ----------------------------------------------------------------------
def plot_frequency_response(ax, pts, folds, om_cm):
    """Frequency response with stable/unstable branches and fold markers."""
    om, amp = pts[:, 1], pts[:, 0]
    stab = np.array([is_stable(a, o) for a, o in zip(amp, om)])

    # Split the branch into maximal stable / unstable runs so that the line
    # style changes exactly at the folds.
    breaks = np.where(stab[:-1] != stab[1:])[0] + 1
    first = True
    for seg in np.split(np.arange(len(om)), breaks):
        if seg.size < 2:
            continue
        if stab[seg[0]]:
            ax.plot(om[seg], amp[seg], "-", color="#1f77b4", lw=2.0,
                    label="stable branch" if first else None)
        else:
            ax.plot(om[seg], amp[seg], "--", color="#d62728", lw=2.0,
                    label="unstable branch" if first else None)
            first = False
    if folds.size:
        ax.plot(folds[:, 1], folds[:, 0], "ko", ms=7, zorder=5,
                label="fold points $F_1$, $F_2$")
        for k, f in enumerate(sorted(folds.tolist(), key=lambda p: p[1])):
            ax.annotate(f"$F_{k + 1}$", (f[1], f[0]),
                        textcoords="offset points", xytext=(8, 6))

    ax.axvline(om_cm, color="0.5", ls=":", lw=1.2)
    ax.text(om_cm + 0.015, 0.08, rf"$\Omega={om_cm:.3f}$", color="0.35",
            fontsize=9)
    ax.set_xlabel(r"excitation frequency $\Omega$")
    ax.set_ylabel(r"amplitude $\hat{x}$")
    ax.set_title("Duffing frequency response\n(pseudo-arc-length continuation)")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.3)


def plot_basins(ax, group, periodic, xc, vc, attractors, om_cm):
    """Basins of attraction from the simple cell map."""
    grid = group[:-1].reshape(N_CELLS, N_CELLS).T   # transpose -> (v, x)
    labels = np.unique(grid)
    colors = ["#9ecae1", "#fdae6b", "#a1d99b", "#c994c7", "#bdbdbd"]
    cmap = ListedColormap(colors[:labels.size])
    norm = BoundaryNorm(np.append(labels - 0.5, labels[-1] + 0.5), labels.size)

    dx = xc[1] - xc[0]
    dv = vc[1] - vc[0]
    edges_x = np.append(xc - dx / 2, xc[-1] + dx / 2)
    edges_v = np.append(vc - dv / 2, vc[-1] + dv / 2)
    ax.pcolormesh(edges_x, edges_v, grid, cmap=cmap, norm=norm)

    handles = []
    for att in attractors:
        ax.plot(att["x"], att["v"], "k*", ms=13, zorder=6)
        colour = colors[list(labels).index(att["group"]) % len(colors)]
        handles.append(Patch(facecolor=colour, edgecolor="k",
                             label=rf"basin of attractor $\hat{{x}}\approx$"
                                   rf"{att['amp']:.2f}"))

    ax.set_xlabel(r"$x(0)$")
    ax.set_ylabel(r"$\dot{x}(0)$")
    ax.set_title(f"SCM basins at $\\Omega={om_cm:.3f}$"
                 f"\n({N_CELLS}$\\times${N_CELLS} cells, "
                 f"{N_PERIODS} periods per cell, stars = periodic groups)")
    ax.legend(handles=handles, loc="lower right", fontsize=8, framealpha=0.9)
    ax.set_aspect("equal")


# ----------------------------------------------------------------------
def main():
    """Run the continuation, run the cell mapping and draw both panels."""
    # ---------------- continuation ----------------
    pts, tangs = continue_branch()
    folds = locate_folds(pts, tangs)
    print(f"continuation: {len(pts)} branch points, "
          f"Omega in [{pts[:, 1].min():.3f}, {pts[:, 1].max():.3f}], "
          f"max amplitude {pts[:, 0].max():.4f}")
    for k, f in enumerate(sorted(folds.tolist(), key=lambda p: p[1])):
        print(f"  fold F{k + 1}: Omega = {f[1]:.5f}, amplitude = {f[0]:.5f}")

    # Pick the frequency for the cell mapping: inside the bistable window
    # bracketed by the two folds (see the comment at OMEGA_CM above).
    if OMEGA_CM is None:
        if folds.shape[0] != 2:
            raise RuntimeError("expected exactly two folds on the branch")
        om_cm = float(np.mean(folds[:, 1]))
    else:
        om_cm = float(OMEGA_CM)

    coex = hb_solutions(om_cm)
    print(f"harmonic balance at Omega = {om_cm:.4f}:")
    for a, st in coex:
        print(f"  amplitude {a:.4f}  ({'stable' if st else 'unstable'})")

    # ---------------- cell mapping ----------------
    image, xc, vc = cell_map_images(om_cm)
    group, steps, periodic = scm_groups(image)

    # Report the periodic groups (the discrete attractors) and their basins.
    cells = group[:-1]
    attractors = []
    print(f"simple cell mapping: {N_CELLS**2} cells, "
          f"{cells.size} mapped, sink cell index {image.size - 1}")
    for g in np.unique(cells):
        mask_basin = cells == g
        mask_per = periodic[:-1] & (cells == g)
        share = 100.0 * mask_basin.sum() / cells.size
        if g == 0:
            print(f"  group 0 (escaped to sink): {share:5.1f} % of the cells")
            continue
        idx = np.where(mask_per)[0]
        if idx.size == 0:
            continue
        i, j = idx // N_CELLS, idx % N_CELLS
        x_a, v_a = xc[i].mean(), vc[j].mean()
        # Rough amplitude of the periodic orbit from its stroboscopic point:
        # for x = a*cos(Omega*t + phi) one has a**2 = x**2 + (x'/Omega)**2.
        amp = np.hypot(xc[i], vc[j] / om_cm).mean()
        attractors.append({"group": int(g), "x": x_a, "v": v_a, "amp": amp})
        print(f"  group {g}: period {idx.size} cell(s), "
              f"centre (x, x') = ({x_a:+.3f}, {v_a:+.3f}), "
              f"approx. amplitude {amp:.3f}, basin {share:5.1f} % of the cells")

    # ---------------- figure ----------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    plot_frequency_response(ax1, pts, folds, om_cm)
    # Link the two panels: the SCM attractors sit on the two stable branches.
    ax1.plot([om_cm] * len(attractors), [a["amp"] for a in attractors],
             "k*", ms=13, zorder=6, label="attractors found by SCM")
    ax1.legend(loc="upper left", fontsize=9)
    plot_basins(ax2, group, periodic, xc, vc, attractors, om_cm)
    fig.suptitle("3.2 Continuation and cell mapping - forced Duffing "
                 rf"($\delta={DELTA}$, $\varepsilon={EPS}$, $\hat{{f}}={FHAT}$)")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
