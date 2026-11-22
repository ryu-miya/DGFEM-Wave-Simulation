import numpy as np
import matplotlib.pyplot as plt
import os
import sys

ne=sys.argv[1]
input_file = f"data/wave_{ne}.dat"
output_dir = f"snapshot/snapshot_{ne}"
os.makedirs(output_dir, exist_ok=True)

with open(input_file) as f:
    blocks = f.read().strip().split("\n\n\n")

fig, ax = plt.subplots(figsize=(8,6))
line_sim = None
line_ana = None

ax.set_xlabel('x[m]')
ax.set_ylabel('u')
ax.set_xlim(0.0, 1.0)
ax.set_ylim(-1.1, 1.1)
ax.grid(True, alpha=0.3)

for k, block in enumerate(blocks):
    lines = block.strip().split("\n")
    t = float(lines[0].split("=")[1])
    frame = np.loadtxt(lines[1:])
    x_ana, x_sim, u_sim, u_ana = frame[:, 0], frame[:, 0].reshape(-1,2), frame[:, 1].reshape(-1,2), frame[:, 2]
    nan_col = np.full((x_sim.shape[0], 1), np.nan) 
    x_sim = np.hstack([x_sim, nan_col]).ravel() 
    u_sim = np.hstack([u_sim, nan_col]).ravel()

    if line_sim is None:
        line_ana, = ax.plot(x_ana, u_ana, color="r", linestyle="-", linewidth=2, label="Analytical")
        line_sim, = ax.plot(x_sim, u_sim, color="b", marker="o", markersize=3, linewidth=1.5, label="DG-FEM")
        ax.legend(loc="upper right")

    else:
        line_ana.set_data(x_ana, u_ana)
        line_sim.set_data(x_sim, u_sim)

    filename = os.path.join(output_dir, f"frame_{k}.png")
    fig.savefig(filename)