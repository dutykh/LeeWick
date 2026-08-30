# Quasinormal Modes of Lee-Wick Black Holes

[![Phys. Rev. D](https://img.shields.io/badge/Phys.%20Rev.%20D-110%20%282024%29%20084032-005a9c)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.110.084032)
[![DOI](https://img.shields.io/badge/DOI-10.1103%2FPhysRevD.110.084032-1a7f37)](https://doi.org/10.1103/PhysRevD.110.084032)
[![arXiv](https://img.shields.io/badge/arXiv-2410.13245-b31b1b)](https://arxiv.org/abs/2410.13245)
[![License: LGPL v2.1](https://img.shields.io/badge/License-LGPL%20v2.1-4c1)](LICENSE)
[![Maple](https://img.shields.io/badge/Maple-18%2B-ff6f00)](https://www.maplesoft.com/products/Maple/)
[![MATLAB](https://img.shields.io/badge/MATLAB-polyeig-0076a8)](https://www.mathworks.com/help/matlab/ref/polyeig.html)
[![QNMs Hall of Fame](https://img.shields.io/badge/QNMs-Hall%20of%20Fame-8a2be2)](https://qnms.denys-dutykh.com/)
[![GitHub](https://img.shields.io/badge/GitHub-dutykh%2FLeeWick-181717?logo=github)](https://github.com/dutykh/LeeWick)
[![Last commit](https://img.shields.io/github/last-commit/dutykh/LeeWick)](https://github.com/dutykh/LeeWick/commits/main)

Computational supplementary material for the study of quasinormal modes (QNMs) of
Lee-Wick black holes by a unified spectral method. The codes assemble and solve the
quadratic eigenvalue problem that governs scalar, electromagnetic and gravitational
perturbations of these geometries, covering non-extremal, extremal and near-extremal
horizon configurations within a single framework.

Repository: <https://github.com/dutykh/LeeWick>

## Reference

D. Batic, D. Dutykh & B. L. Giacchini, *A Unified Spectral Approach for Quasinormal
Modes of Lee-Wick Black Holes*, **Physical Review D 110**, 084032 (2024).

* Published version: <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.110.084032>
* Preprint: <https://arxiv.org/abs/2410.13245>

> In this paper, we undertake a comprehensive examination of quasinormal modes linked
> to Lee-Wick black holes, delving into scalar, electromagnetic, and gravitational
> perturbations using the spectral method. Such black holes can display a rich
> structure of horizons, and our analysis considers all the representative scenarios,
> including extremal and non-extremal situations. In particular, we show that purely
> imaginary quasinormal modes emerge for extremal and near-extremal configurations,
> suggesting a rapid return to equilibrium without oscillation.

## Quasi-Normal Modes Hall of Fame

This work is part of the [Quasi-Normal Modes Hall of Fame](https://qnms.denys-dutykh.com/),
a curated collection of studies on the "celestial music" of compact objects: the
characteristic ringing frequencies of black holes and wormholes in general relativity and
in its modified and quantum-corrected extensions. The collection is maintained by
Dr. Davide Batic and Dr. Denys Dutykh (Mathematics Department, Khalifa University of
Science and Technology, Abu Dhabi, UAE) and gathers, for each published manuscript, the
preprint, the journal reference and, where available, the accompanying computational
material. The present repository is the computational entry associated with the Lee-Wick
paper listed there, and the Hall of Fame is the place to look for the companion codes of
the neighbouring spacetimes treated with the same spectral machinery.

## Authors

| Author | Affiliation |
| :--- | :--- |
| Dr. Davide Batic | Khalifa University of Science and Technology, Abu Dhabi, UAE |
| Dr. Denys Dutykh | Khalifa University of Science and Technology, Abu Dhabi, UAE |
| Dr. Breno L. Giacchini | Charles University, Prague, Czech Republic |

## Physical and numerical setting

The Lee-Wick black hole is described by the lapse function

```
f(x) = 1 - h(x)/x,
h(x) = 1 - e^(-a x) [ (1 + a(q^2 + 1) x / 2) cos(q a x)
                      + (1 - q^2 + a(q^2 + 1) x) sin(q a x) / (2 q) ],
```

where `x` is the dimensionless radial coordinate, `a` (also written `alpha`) is the
dimensionless mass parameter and `q` is the ratio fixing the complex-conjugate pair of
Lee-Wick masses. Depending on `(q, a)`, the geometry may possess two horizons, a single
degenerate horizon, or no horizon at all.

Master perturbation equations for spin `s = 0, 1, 2` are compactified onto `y` in
`[-1, 1]`, with the horizon and spatial infinity mapped to the interval endpoints, and
the QNM boundary conditions built into the ansatz. The field is expanded in Chebyshev
polynomials of the first kind,

```
F(y) = sum_{j=0}^{n-1} a_j T_j(y),
```

and the equation is collocated at the Chebyshev roots
`y_i = cos((2i - 1) pi / (2n))`, `i = 1, ..., n`. Since the frequency `Omega` enters the
master equation quadratically, collocation yields a quadratic eigenvalue problem

```
( M0 + Omega M1 + Omega^2 M2 ) a = 0,
```

whose eigenvalues are the quasinormal frequencies. The matrices are assembled in Maple
with arbitrary-precision arithmetic and exported in MATLAB format; the polynomial
eigenvalue problem is then solved in MATLAB with `polyeig`.

## Repository contents

### Maple worksheets (derivations)

| File | Content |
| :--- | :--- |
| `00-f(x).mw` | Analysis of the metric function `f(x)`: horizon structure, critical parameter values, and the degenerate (extremal) configuration in the `(q, alpha)` plane. |
| `01-NONEXTREME CASE-SUPPLEMENTAL MATERIAL.mw` | Full derivation for a non-degenerate horizon: change of variables, asymptotic behaviour of the solution at both endpoints, and the operator coefficients entering `M0`, `M1`, `M2`. |
| `LW-EXTREME CASE-SUPPLEMENTAL MATERIAL.mw` | Same derivation for the extremal (degenerate horizon) case, where the surface gravity vanishes and the standard tortoise mapping degenerates. Introduces the adapted coordinate `eta(y)` and the near-horizon coefficients `h1` and `a`. |

### Maple procedures (matrix assembly)

All three files define a procedure `MatrixAssembler` that builds `M0`, `M1`, `M2` and
writes them as MATLAB ASCII matrices.

| File | Case | Signature |
| :--- | :--- | :--- |
| `matrixassemblerNonExtreme.mpl` | Non-extremal black hole; the horizon `xh` is obtained from `x = h(x)` for user-supplied `(q, a)`. | `MatrixAssembler(d, n, s, L, q, a, p)` |
| `matrixassemblerTypeA.mpl` | Type A extremal configuration: the critical pair `(x_e, alpha)` is found from `f = 0`, `df/dx = 0`, and the perturbation problem is posed at the outer non-degenerate horizon `xh`. | `MatrixAssembler(d, n, s, L, q, p)` |
| `matrixassemblerTypeB.mpl` | Type B extremal configuration: the problem is posed at the degenerate horizon `x_e` itself, using the adapted coordinate `eta(y)`. | `MatrixAssembler(d, n, s, L, q, p)` |

Arguments:

* `d` : number of digits used in the computation (Maple `Digits`);
* `n` : number of Chebyshev modes, equal to the size of the assembled matrices;
* `s` : spin of the perturbation, `0` (scalar), `1` (electromagnetic) or `2` (gravitational);
* `L` : angular momentum multipole index, with `L >= s`;
* `q` : ratio of the Lee-Wick masses;
* `a` : dimensionless mass parameter (non-extremal case only; in the extremal cases it is determined by the degeneracy conditions);
* `p` : path to the directory containing the output subdirectory `data/`.

In the Type A file the computed critical values are checked against the reference values
`x_e = 0.62987`, `alpha = 5.18635` and `xh = 1.040133223`, and the procedure aborts if the
solver drifts away from them.

## Usage

### 1. Assemble the matrices in Maple

Create the output subdirectory `data/` first, then, from a Maple session,

```maple
read "matrixassemblerNonExtreme.mpl":
MatrixAssembler(30, 60, 2, 2, 0.5, 1.0, "/path/to/LeeWick"):
```

which produces `data/M0_60.mat`, `data/M1_60.mat` and `data/M2_60.mat`. The extremal
cases are invoked in the same way, without the argument `a`:

```maple
read "matrixassemblerTypeB.mpl":
MatrixAssembler(30, 60, 2, 2, 0.5, "/path/to/LeeWick"):
```

A large value of `d` is advisable: the collocation matrices of a spectral discretization
become severely ill-conditioned as `n` grows, and the accuracy of the high overtones
depends directly on the working precision.

### 2. Solve the quadratic eigenvalue problem in MATLAB

```matlab
n  = 60;
M0 = load(sprintf('data/M0_%d.mat', n));
M1 = load(sprintf('data/M1_%d.mat', n));
M2 = load(sprintf('data/M2_%d.mat', n));

[~, Omega] = polyeig(M0, M1, M2);
Omega = Omega(isfinite(Omega));
```

The spectrum returned by `polyeig` contains both the physical quasinormal frequencies and
spurious eigenvalues introduced by the discretization. The physical modes are identified
by repeating the computation for two or more values of `n` and retaining the eigenvalues
that remain stable to the requested number of digits.

## Software requirements

* Maple (TM) 18 or later, for the worksheets and the assembly procedures;
* MATLAB (TM), for the polynomial eigenvalue solver `polyeig`.

## Citation

```bibtex
@article{Batic2024,
  author  = {Batic, Davide and Dutykh, Denys and Giacchini, Breno L.},
  title   = {A Unified Spectral Approach for Quasinormal Modes of {L}ee-{W}ick Black Holes},
  journal = {Physical Review D},
  volume  = {110},
  number  = {8},
  pages   = {084032},
  year    = {2024},
  doi     = {10.1103/PhysRevD.110.084032},
  eprint  = {2410.13245},
  archivePrefix = {arXiv},
  primaryClass  = {gr-qc}
}
```

## License

Distributed under the GNU Lesser General Public License v2.1. See [`LICENSE`](LICENSE).

**NB**: Please acknowledge this work if you use these codes.
