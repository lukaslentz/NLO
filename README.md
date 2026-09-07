# NLO — Nonlinear Oscillations

Python examples accompanying the lecture **Nonlinear Oscillations** (Winter Semester 2026/27)
at Umwelt-Campus Birkenfeld, Hochschule Trier.

Prof. Dr.-Ing. Lukas Lentz

Each script belongs to one chapter of the lecture and is referenced from the corresponding
slide. The slides show a compressed version that fits on a frame; the scripts here are the
full, commented and runnable version.

## Running the examples

```bash
pip install -r requirements.txt
python python/ch1_1_nonlinear_systems.py
```

Every script is standalone: no arguments, no data files, no imports from the other scripts.
Each opens one figure window and prints the key numerical results to the console.
All of them finish in well under a minute; where a resolution can be increased at the cost
of runtime (cell-mapping grids, stability charts, basin resolutions), a comment says so.

Requirements: Python 3.9 or newer, `numpy`, `scipy`, `matplotlib`.

## Contents

| Chapter | Script | Topic |
|---|---|---|
| 1.1 | [`ch1_1_nonlinear_systems.py`](python/ch1_1_nonlinear_systems.py) | Van der Pol limit cycle, amplitude-dependent frequency |
| 1.2 | [`ch1_2_duffing_oscillator.py`](python/ch1_2_duffing_oscillator.py) | Duffing frequency response, double-well potential and phase portrait |
| 2.1 | [`ch2_1_phase_space.py`](python/ch2_1_phase_space.py) | Direction fields, trajectories, classification of equilibria |
| 2.2 | [`ch2_2_energy_pendulum.py`](python/ch2_2_energy_pendulum.py) | Pendulum energy contours, separatrix, energy exchange |
| 2.3 | [`ch2_3_large_amplitude.py`](python/ch2_3_large_amplitude.py) | Exact pendulum period, elliptic integrals, backbone curve |
| 3.1 | [`ch3_1_numerical_integration.py`](python/ch3_1_numerical_integration.py) | Euler / symplectic Euler / RK4, convergence order, energy drift |
| 3.2 | [`ch3_2_continuation_cell_mapping.py`](python/ch3_2_continuation_cell_mapping.py) | Pseudo-arc-length continuation, simple cell mapping, basins |
| 4.1 | [`ch4_1_perturbation_methods.py`](python/ch4_1_perturbation_methods.py) | Secular terms, Lindstedt-Poincaré, backbone curve |
| 4.2 | [`ch4_2_harmonic_balance.py`](python/ch4_2_harmonic_balance.py) | Multi-harmonic balance with AFT, convergence in the harmonic order |
| 4.3 | [`ch4_3_equivalent_linearization.py`](python/ch4_3_equivalent_linearization.py) | Describing functions by Fourier integration, equivalent stiffness |
| 4.4 | [`ch4_4_forced_duffing.py`](python/ch4_4_forced_duffing.py) | Frequency sweeps, jump phenomenon, hysteresis, Poincaré sections |
| 4.5 | [`ch4_5_multiple_scales.py`](python/ch4_5_multiple_scales.py) | Slow-flow equations, modulation envelope, steady-state stability |
| 5.1 | [`ch5_1_poincare_maps.py`](python/ch5_1_poincare_maps.py) | Stroboscopic sections, period doubling, strange attractor |
| 5.2 | [`ch5_2_floquet_theory.py`](python/ch5_2_floquet_theory.py) | Monodromy matrix, Floquet multipliers, Ince-Strutt chart |
| 5.3 | [`ch5_3_parametric_self_excited.py`](python/ch5_3_parametric_self_excited.py) | Mathieu stability tongues, van der Pol energy balance |
| 6.1 | [`ch6_1_bifurcations.py`](python/ch6_1_bifurcations.py) | Saddle-node, transcritical, pitchfork, Hopf, period doubling |
| 6.2 | [`ch6_2_chaos_fractal_basins.py`](python/ch6_2_chaos_fractal_basins.py) | Lyapunov exponent, sensitive dependence, fractal basin boundary |

## A note on the models

Where the compressed slide listing and the physics disagreed, the script follows the
physics and says so in its module docstring. The clearest case is chaos in the Duffing
oscillator: the single-well form is never chaotic, so chapters 5.1 and 6.2 use the twin-well
(Holmes) form. Each such deviation is documented at the top of the script concerned.
