"""
Chapter 2.2 -- Energy and the Pendulum

Energy-based analysis of the undamped nonlinear pendulum, following the slide
set "2.2 Energy and the Pendulum".  The pendulum of mass m and length l with
angle phi obeys

    phi'' + omega0^2 sin(phi) = 0,        omega0 = sqrt(g/l),

which is conservative: the total mechanical energy

    E = T + V = 1/2 m l^2 phi'^2 + m g l (1 - cos phi)

is a first integral, dE/dt = 0.  The phase portrait therefore needs no time
integration at all -- the level curves E(phi, phi') = const ARE the
trajectories.  The critical level through the saddles at phi = +-pi,
E_sep = 2 m g l, is the separatrix (a homoclinic orbit) and divides libration
(E < E_sep, closed orbits) from rotation (E > E_sep, open orbits).

The script produces a four-panel figure:

  (a) Phase portrait built from energy contours in the (phi, phi') plane, with
      the separatrix highlighted, the centres at phi = 2k*pi and the saddles at
      phi = (2k+1)*pi marked, and one librating, one rotating orbit obtained by
      actual time integration to confirm that the contours are trajectories.
  (b) Potential energy V(phi) = m g l (1 - cos phi) over several wells with
      horizontal energy levels; the intersections V(phi) = E are the turning
      points of the corresponding orbit.
  (c) Energy exchange over exactly one period of a librating orbit: kinetic
      energy T(t), potential energy V(t) and their constant sum.  T and V swap
      twice per period; the numerically constant sum verifies the integrator.
  (d) Exact period T = 4 K(k)/omega0 with k = sin(phi_max/2), the complete
      elliptic integral of the first kind, compared with the small-amplitude
      series T0 (1 + phi_max^2/16 + 11 phi_max^4/3072) and with periods
      measured from the time integration.  T -> infinity as phi_max -> pi.

The console output lists the equilibria with their Jacobian eigenvalues, the
separatrix energy, the energy-conservation error of the integrator and the
period-ratio table of the slides.

Runtime: a few seconds.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import ellipk

# ---------------------------------------------------------------------------
# PARAMETERS -- change these to explore
# ---------------------------------------------------------------------------
G = 9.81                 # gravitational acceleration [m/s^2]
L = 1.0                  # pendulum length [m]
M = 1.0                  # bob mass [kg]

PHI_EXCHANGE = np.deg2rad(120.0)   # amplitude of the orbit shown in panel (c)
PHI_LIB_DEMO = np.deg2rad(150.0)   # librating demo orbit in panel (a)
PHIDOT_ROT_DEMO_FACTOR = 1.15      # rotating demo orbit: phi' at phi = 0 as a
                                   # multiple of the separatrix value 2*omega0

N_GRID = 500             # contour grid points per axis (cheap; 500 is plenty)
PHI_WINDOW = 2.5 * np.pi # phase-portrait window in phi [rad]
PHIDOT_WINDOW = 9.0      # phase-portrait window in phi' [rad/s]
N_LEVELS = 14            # number of energy contours below/above the separatrix

# Amplitudes of the period-ratio table on the slides.  The exact values printed
# below reproduce the slide table up to 90 deg (1.002, 1.017, 1.073, 1.180); at
# 150 deg and 170 deg the exact ratios are 1.762 and 2.439, slightly larger than
# the 1.654 / 2.257 quoted there.  Both the elliptic integral and the direct
# time integration agree on the values printed here.
TABLE_DEG = [10.0, 30.0, 60.0, 90.0, 150.0, 170.0]
RTOL, ATOL = 1e-11, 1e-12                            # tight: E must stay const
# ---------------------------------------------------------------------------

OMEGA0 = np.sqrt(G / L)          # natural (small-amplitude) frequency [rad/s]
T0 = 2.0 * np.pi / OMEGA0        # small-amplitude period [s]
E_SEP = 2.0 * M * G * L          # separatrix energy [J]


def pendulum(t, y):
    """Undamped pendulum phi'' = -omega0^2 sin(phi) as a first-order system."""
    phi, phidot = y
    return [phidot, -OMEGA0**2 * np.sin(phi)]


def potential(phi):
    """Potential energy V(phi) = m g l (1 - cos phi) [J], zero at the bottom."""
    return M * G * L * (1.0 - np.cos(phi))


def kinetic(phidot):
    """Kinetic energy T = 1/2 m l^2 phi'^2 [J]."""
    return 0.5 * M * L**2 * phidot**2


def energy(phi, phidot):
    """Total mechanical energy E = T + V [J]; a first integral of the motion."""
    return kinetic(phidot) + potential(phi)


def exact_period(phi_max):
    """Exact libration period T = 4 K(k) / omega0 with modulus k = sin(phi_max/2).

    SciPy's ellipk takes the PARAMETER m = k^2, not the modulus k -- a classic
    source of factor errors.
    """
    k = np.sin(phi_max / 2.0)
    return 4.0 / OMEGA0 * ellipk(k**2)


def series_period(phi_max):
    """Small-amplitude series T0 (1 + phi^2/16 + 11 phi^4/3072) of the slides."""
    return T0 * (1.0 + phi_max**2 / 16.0 + 11.0 * phi_max**4 / 3072.0)


def measured_period(phi_max):
    """Measure the libration period by integrating from rest at phi = phi_max.

    Released at the turning point, the pendulum needs half a period to reach
    the opposite turning point -phi_max.  The event is placed on phi' = 0
    rather than on phi = -phi_max: at the turning point phi touches -phi_max
    tangentially (phi' = 0 there) and a root finder would miss the sign change,
    whereas phi' crosses zero transversally.  Catching the event gives the
    half period to solver accuracy instead of to output-sampling accuracy.
    """
    def turning(t, y):
        return y[1]
    turning.terminal = True
    turning.direction = +1      # phi' returns to zero from below

    t_max = 4.0 * exact_period(phi_max)          # generous upper bound
    sol = solve_ivp(pendulum, (0.0, t_max), [phi_max, 0.0], events=turning,
                    rtol=RTOL, atol=ATOL, method="DOP853")
    return 2.0 * sol.t_events[0][0]


def equilibrium_report():
    """Print the equilibria of the pendulum and their Jacobian eigenvalues.

    J = [[0, 1], [-omega0^2 cos(phi*), 0]]: at phi* = 2k*pi the eigenvalues are
    +-i*omega0 (a centre, and here the nonlinear system really has closed orbits
    because E is conserved), at phi* = (2k+1)*pi they are +-omega0 (a saddle).
    """
    print("Nonlinear pendulum:  phi'' + omega0^2 sin(phi) = 0")
    print(f"  g = {G} m/s^2, l = {L} m, m = {M} kg")
    print(f"  omega0 = {OMEGA0:.4f} rad/s,  T0 = {T0:.4f} s")
    print(f"  separatrix energy E_sep = 2 m g l = {E_SEP:.4f} J\n")
    print("Equilibria (V'(phi) = m g l sin(phi) = 0):")
    for phi_star, name in ((0.0, "phi* = 0"), (np.pi, "phi* = pi")):
        J = np.array([[0.0, 1.0], [-OMEGA0**2 * np.cos(phi_star), 0.0]])
        lam = np.linalg.eigvals(J)
        vpp = M * G * L * np.cos(phi_star)
        kind = "centre (stable)" if vpp > 0 else "saddle (unstable)"
        print(f"  {name:12s} V'' = {vpp:+.4f} J  ->  {kind}, "
              f"lambda = {lam[0]:+.4f}, {lam[1]:+.4f}")
    print()


def panel_phase_portrait(ax):
    """(a) Energy contours = trajectories, separatrix, centres and saddles."""
    phi = np.linspace(-PHI_WINDOW, PHI_WINDOW, N_GRID)
    phidot = np.linspace(-PHIDOT_WINDOW, PHIDOT_WINDOW, N_GRID)
    PHI, PHIDOT = np.meshgrid(phi, phidot)
    E = energy(PHI, PHIDOT)

    ax.contourf(PHI, PHIDOT, E, levels=40, cmap="viridis", alpha=0.75)
    # Libration levels strictly below E_sep, rotation levels strictly above.
    ax.contour(PHI, PHIDOT, E, colors="white", linewidths=0.7,
               levels=np.linspace(0.05 * E_SEP, 0.95 * E_SEP, N_LEVELS))
    ax.contour(PHI, PHIDOT, E, colors="white", linewidths=0.7,
               levels=np.linspace(1.1 * E_SEP, energy(np.pi, PHIDOT_WINDOW),
                                  N_LEVELS))
    ax.contour(PHI, PHIDOT, E, levels=[E_SEP], colors="red", linewidths=2.2)

    # Confirm by integration that a contour really is a trajectory.
    sol = solve_ivp(pendulum, (0.0, exact_period(PHI_LIB_DEMO)),
                    [PHI_LIB_DEMO, 0.0], rtol=RTOL, atol=ATOL,
                    dense_output=True, method="DOP853")
    t = np.linspace(0.0, sol.t[-1], 2000)
    ax.plot(*sol.sol(t)[:2], "k-", lw=1.8,
            label=r"Libration ($E < E_\mathrm{sep}$)")

    phidot_rot = PHIDOT_ROT_DEMO_FACTOR * 2.0 * OMEGA0
    sol = solve_ivp(pendulum, (0.0, 6.0), [-PHI_WINDOW, phidot_rot],
                    rtol=RTOL, atol=ATOL, dense_output=True, method="DOP853")
    t = np.linspace(0.0, sol.t[-1], 2000)
    y = sol.sol(t)
    keep = np.abs(y[0]) <= PHI_WINDOW
    ax.plot(y[0][keep], y[1][keep], "k-.", lw=1.8,
            label=r"Rotation ($E > E_\mathrm{sep}$)")

    for k in (-1, 0, 1):
        ax.plot(2 * k * np.pi, 0.0, "o", color="white", mec="black", ms=8,
                zorder=5, label="Centre" if k == 0 else None)
        ax.plot((2 * k + 1) * np.pi, 0.0, "s", color="red", mec="black", ms=8,
                zorder=5, label="Saddle" if k == 0 else None)

    ax.plot([], [], "r-", lw=2.2,
            label=rf"Separatrix $E={E_SEP:.2f}$ J")
    ax.set_xlim(-PHI_WINDOW, PHI_WINDOW)
    ax.set_ylim(-PHIDOT_WINDOW, PHIDOT_WINDOW)
    ax.set_xticks(np.pi * np.arange(-2, 3))
    ax.set_xticklabels([r"$-2\pi$", r"$-\pi$", "0", r"$\pi$", r"$2\pi$"])
    ax.set_xlabel(r"$\varphi$ [rad]")
    ax.set_ylabel(r"$\dot{\varphi}$ [rad/s]")
    ax.set_title("(a) Phase portrait from energy contours")
    ax.legend(loc="upper right", fontsize=7, framealpha=0.9)


def panel_potential(ax):
    """(b) Potential well V(phi) with energy levels and turning points."""
    phi = np.linspace(-PHI_WINDOW, PHI_WINDOW, 1000)
    ax.plot(phi, potential(phi), color="tab:blue", lw=2.0,
            label=r"$V(\varphi)=mgl(1-\cos\varphi)$")
    ax.axhline(E_SEP, color="red", lw=2.0,
               label=r"$E_\mathrm{sep}=2mgl$ (separatrix)")

    for frac, style in ((0.25, "--"), (0.6, "--"), (1.4, ":")):
        E = frac * E_SEP
        ax.axhline(E, color="0.4", ls=style, lw=1.0)
        if frac < 1.0:   # turning points exist only below the separatrix
            phi_t = np.arccos(1.0 - E / (M * G * L))
            ax.plot([-phi_t, phi_t], [E, E], "o", color="tab:orange", ms=5)
        ax.text(PHI_WINDOW * 0.98, E + 0.4, f"E = {E:.1f} J",
                ha="right", fontsize=7, color="0.3")

    ax.set_xlim(-PHI_WINDOW, PHI_WINDOW)
    ax.set_xticks(np.pi * np.arange(-2, 3))
    ax.set_xticklabels([r"$-2\pi$", r"$-\pi$", "0", r"$\pi$", r"$2\pi$"])
    ax.set_xlabel(r"$\varphi$ [rad]")
    ax.set_ylabel("Energy [J]")
    ax.set_title("(b) Potential energy and energy levels\n"
                 "(orange: turning points $V(\\varphi)=E$)")
    ax.legend(loc="upper left", fontsize=7, framealpha=0.9)


def panel_energy_exchange(ax):
    """(c) Kinetic and potential energy over exactly one libration period."""
    period = exact_period(PHI_EXCHANGE)
    sol = solve_ivp(pendulum, (0.0, period), [PHI_EXCHANGE, 0.0],
                    t_eval=np.linspace(0.0, period, 2000),
                    rtol=RTOL, atol=ATOL, method="DOP853")
    t, (phi, phidot) = sol.t, sol.y
    T, V = kinetic(phidot), potential(phi)
    E = T + V

    ax.plot(t, T, color="tab:red", lw=1.8, label=r"Kinetic $T$")
    ax.plot(t, V, color="tab:blue", lw=1.8, label=r"Potential $V$")
    ax.plot(t, E, color="black", lw=2.0, ls="--", label=r"Total $E=T+V$")
    ax.set_xlim(0.0, period)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel("Time $t$ [s]")
    ax.set_ylabel("Energy [J]")
    ax.set_title(rf"(c) Energy exchange, $\varphi_\mathrm{{max}}="
                 rf"{np.degrees(PHI_EXCHANGE):.0f}^\circ$, $T={period:.3f}$ s")
    ax.legend(loc="center right", fontsize=8, framealpha=0.9)

    drift = np.max(np.abs(E - E[0])) / E[0]
    print(f"Energy exchange panel (phi_max = "
          f"{np.degrees(PHI_EXCHANGE):.0f} deg):")
    print(f"  E = {E[0]:.6f} J,  max|T| = {T.max():.6f} J,  "
          f"max|V| = {V.max():.6f} J")
    print(f"  relative energy drift of the integrator: {drift:.2e}")
    print(f"  exact period T = {period:.6f} s,  T/T0 = {period / T0:.4f}\n")


def panel_period(ax):
    """(d) Backbone curve T(phi_max): elliptic integral, series and measurement."""
    phi_max = np.linspace(0.01, np.pi - 1e-3, 400)
    ax.plot(np.degrees(phi_max), exact_period(phi_max) / T0,
            color="tab:blue", lw=2.0, label=r"Exact $4K(k)/\omega_0$")
    ax.plot(np.degrees(phi_max), series_period(phi_max) / T0,
            color="tab:green", ls="--", lw=1.5,
            label=r"Series $1+\varphi^2/16+11\varphi^4/3072$")
    ax.axhline(1.0, color="gray", ls=":", lw=1.2,
               label=r"Linear $T_0=2\pi/\omega_0$")

    print("Period ratio T/T0 (slide table):")
    print("  phi_max [deg]      k     T/T0 exact   T/T0 measured   series")
    for deg in TABLE_DEG:
        a = np.deg2rad(deg)
        r_exact = exact_period(a) / T0
        r_meas = measured_period(a) / T0
        print(f"  {deg:9.0f}   {np.sin(a / 2):.3f}     {r_exact:8.4f}"
              f"      {r_meas:9.4f}   {series_period(a) / T0:8.4f}")
        ax.plot(deg, r_meas, "o", color="crimson", ms=6, zorder=5,
                label="Measured by integration" if deg == TABLE_DEG[0] else None)
    print()

    ax.set_xlim(0.0, 180.0)
    ax.set_ylim(0.9, 3.2)
    ax.set_xlabel(r"Amplitude $\varphi_\mathrm{max}$ [deg]")
    ax.set_ylabel(r"$T / T_0$")
    ax.set_title(r"(d) Backbone curve: $T\to\infty$ as "
                 r"$\varphi_\mathrm{max}\to\pi$")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)


def main():
    equilibrium_report()
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.0))
    panel_phase_portrait(axes[0, 0])
    panel_potential(axes[0, 1])
    panel_energy_exchange(axes[1, 0])
    panel_period(axes[1, 1])
    fig.suptitle("2.2 Energy and the pendulum: contours, separatrix, "
                 "energy exchange and amplitude-dependent period")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
