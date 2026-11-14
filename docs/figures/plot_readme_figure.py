#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: the lapse function f(x) = 1 - h(x)/x of Lee-Wick black holes.

The extremal configurations are found as the assemblers find them, by solving f = 0 and
df/dx = 0 together for the double root x_e and the mass parameter alpha, here with SciPy.
From this folder,
    python3 plot_readme_figure.py
writes leewick-lapse.png (Python 3 with NumPy, SciPy and Matplotlib). The style is the stamped copy
manuscript_style.py beside it.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.optimize import brentq, fsolve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402


def h(x, a, q):
    """The function h of equation (11) of the article, as in matrixassemblerNonExtreme.mpl."""
    return 1 - np.exp(-a * x) * ((1 + a * (1 + q**2) * x / 2) * np.cos(a * q * x)
                                 + (1 - q**2 + a * (1 + q**2) * x) * np.sin(a * q * x) / (2 * q))


def f(x, a, q):
    return 1 - h(x, a, q) / x


def dfdx(x, a, q, step=1e-6):
    return (f(x + step, a, q) - f(x - step, a, q)) / (2 * step)


def double_root(q, guess):
    """The pair (x_e, alpha_e) where f and df/dx vanish together."""
    solution, info, flag, _ = fsolve(lambda v: [f(v[0], v[1], q), dfdx(v[0], v[1], q)],
                                     guess, xtol=1e-13, full_output=True)
    if flag != 1:
        raise RuntimeError(f"no double root found near {guess} for q = {q}")
    return solution


hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 260
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig, (left, right) = plt.subplots(1, 2, figsize=(hs.WIDTH, 2.75))

# (a) q = 2: the three extremal values of alpha (Table I of the article).
q = 2.0
x = np.linspace(1e-3, 2.5, 5000)
cases = [((1.77, 0.69), "type B"), ((1.07, 4.05), "type B"), ((0.63, 5.19), "type A")]
for k, (guess, kind) in enumerate(cases):
    xe, alpha = double_root(q, guess)
    style = hs.series(k)
    style["marker"] = None
    left.plot(x, f(x, alpha, q), lw=1.1, label=rf"$\alpha = {alpha:.5f}$ ({kind})", **style)
    left.plot([xe], [0.0], linestyle="none", marker="o", markersize=4.2,
              markerfacecolor="white", markeredgecolor=hs.SERIES[k], markeredgewidth=0.9)
    print(f"q = {q}: double root x_e = {xe:.8f} at alpha = {alpha:.8f} ({kind})")
    if kind == "type A":
        xh = brentq(lambda t: f(t, alpha, q), 0.9, 1.3, xtol=1e-14)
        print(f"          event horizon of the type A black hole: x_h = {xh:.8f}")
left.axhline(0.0, color=hs.SLATE, lw=0.6, linestyle=(0, (1.6, 1.6)), zorder=0)
left.set_xlim(0, 2.5)
left.set_ylim(-0.75, 1.1)
left.set_xlabel(r"$x = r/(2M)$")
left.set_ylabel(r"$f(x)$")
hs.panel(left, r"(a) $q = 2$, the three extremal values of $\alpha$")
hs.legend_outside(left, where="below", ncol=1)

# (b) q = 0.1: below, at and above the extremal value of alpha.
q = 0.1
xe, alpha_e = double_root(q, (0.66, 5.09))
print(f"q = {q}: double root x_e = {xe:.8f} at alpha = {alpha_e:.8f}")
x = np.linspace(1e-3, 3.0, 3000)
for k, (alpha, words) in enumerate([(4.1, "no horizon"), (alpha_e, "extremal"), (6.1, "two horizons")]):
    style = hs.series(k)
    style["marker"] = None
    right.plot(x, f(x, alpha, q), lw=1.1, label=rf"$\alpha = {alpha:.6g}$ ({words})", **style)
right.plot([xe], [0.0], linestyle="none", marker="o", markersize=4.2,
           markerfacecolor="white", markeredgecolor=hs.SERIES[1], markeredgewidth=0.9)
right.axhline(0.0, color=hs.SLATE, lw=0.6, linestyle=(0, (1.6, 1.6)), zorder=0)
right.set_xlim(0, 3.0)
right.set_ylim(-0.3, 1.1)
right.set_xlabel(r"$x = r/(2M)$")
right.set_ylabel(r"$f(x)$")
hs.panel(right, r"(b) $q = 0.1$, around its extremal value of $\alpha$")
hs.legend_outside(right, where="below", ncol=1)

out = HERE / "leewick-lapse.png"
hs.save(fig, out)
# The background is opaque: without its alpha channel, and compressed a second time without
# loss, the file is smaller.
Image.open(out).convert("RGB").save(out, optimize=True)
print("written", out)
