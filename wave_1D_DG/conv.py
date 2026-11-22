import numpy as np
import matplotlib.pyplot as plt

L = 1.0

# read lines like "Elements = 50, error = 8.142781e-01"
ne_list = []
err_list = []
with open("conv.dat") as f:
    for line in f:
        parts = line.replace(",", "").split()
        if len(parts) < 6:          # skip empty or broken lines
            continue
        ne_list.append(int(parts[2]))
        err_list.append(float(parts[5]))

ne = np.array(ne_list)
err = np.array(err_list)
h = L / ne                          # element size

# convergence order between neighboring meshes
for i in range(1, len(ne)):
    p = np.log(err[i-1] / err[i]) / np.log(ne[i] / ne[i-1])   # [FIX] h[i-1]/h[i] -> ne[i]/ne[i-1]
    print(f"Ne = {ne[i-1]:4d} -> {ne[i]:4d} : order = {p:.2f}")

fig, ax = plt.subplots(figsize=(6,5))
ax.loglog(ne, err, "o-", color="tab:blue", label="max error")                      # [FIX] h -> ne
ax.loglog(ne, err[-1] * (ne / ne[-1])**(-2), "k--", label="slope -2 (2nd order)")  # [FIX] reference line
ax.set_xlabel("number of elements Ne")                                             # [FIX]
ax.set_ylabel("error")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
fig.savefig("convergence.png", dpi=100)