"""Plot successive local maxima of x for the classical Lorenz system."""

from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--positive-only", action="store_true",
                    help="Pair successive maxima after retaining only x > 0 events.")
args = parser.parse_args()

def lorenz(t, state):
    x, y, z = state
    return [10.0 * (y - x), x * (28.0 - z) - y, x * y - (8.0 / 3.0) * z]


def x_maximum(t, state):
    return lorenz(t, state)[0]


x_maximum.direction = -1
x_maximum.terminal = False

solution = solve_ivp(
    lorenz, (0.0, 2000.0), [1.0, 1.0, 1.0],
    method="DOP853", rtol=1e-10, atol=1e-12, max_step=0.02,
    events=x_maximum,
)
assert solution.success, solution.message
times = solution.t_events[0]
states = solution.y_events[0]
keep = times > 100.0
times, states = times[keep], states[keep]
maxima = states[:, 0]
# At x = y, the second derivative of x is 10*x*(27-z).
curvatures = 10.0 * states[:, 0] * (27.0 - states[:, 2])
assert np.all(curvatures < 0.0)
assert np.max(np.abs(states[:, 1] - states[:, 0])) < 1e-7

if args.positive_only:
    keep = maxima > 0.0
    times, maxima = times[keep], maxima[keep]

fig, ax = plt.subplots(figsize=(8.4, 7.2), constrained_layout=True)
ax.scatter(maxima[:-1], maxima[1:], s=5, color="#244cce", alpha=0.7, linewidths=0)
lo, hi = maxima.min() - 1, maxima.max() + 1
ax.plot([lo, hi], [lo, hi], ":", color="0.55", linewidth=1.2)
ax.set(xlim=(lo, hi), ylim=(lo, hi), xlabel=r"$x_n$", ylabel=r"$x_{n+1}$")
heading = "Lorenz: return map of positive local maxima of x" if args.positive_only else "Lorenz: return map of local maxima of x"
ax.set_title(heading + "\n"
             r"$\sigma=10,\ \rho=28,\ \beta=8/3$", fontsize=15)
ax.grid(alpha=0.18)
ax.set_aspect("equal", adjustable="box")
event_label = "successive accepted events with x > 0" if args.positive_only else "maxima of x, not |x|"
fig.supxlabel(f"{len(maxima)-1:,} pairs · transient t ≤ 100 discarded\n{event_label}", fontsize=10)
destination = Path(__file__).resolve().parents[1] / "output" / "figures"
destination.mkdir(parents=True, exist_ok=True)
stem = "lorenz-positive-x-maxima" if args.positive_only else "lorenz-x-maxima"
fig.savefig(destination / f"{stem}-return-map.png", dpi=180)
np.savetxt(destination / f"{stem}.csv", np.column_stack((times, maxima)),
           delimiter=",", header="time,x_maximum", comments="")
print(f"Recorded {len(maxima)} maxima; {np.sum(maxima < 0)} are negative.")
print(destination / f"{stem}-return-map.png")
