import numpy as np
import matplotlib.pyplot as plt
import os

dirs_compare = next(os.walk("."))[1]

ode_lib = {
    "ascent": "Ascent",
    "diffeq": "DifferentialEquations.jl",
    "gsl": "GSL",
    "odeint": "Boost::odeint",
    "ponio": "ponio",
    "scipy": "SciPy"
}


# 3d plot
fig = plt.figure(constrained_layout=True, figsize=(8, 8))
gs = fig.add_gridspec(6, 1)
axlist = []

ax = fig.add_subplot(gs[:, 0], projection='3d')
ax.set_xlim3d([-20, 20])
ax.set_ylim3d([-20, 20])
ax.set_zlim3d(bottom=0, top=50)

ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.xaxis.pane.set_edgecolor('w')
ax.yaxis.pane.set_edgecolor('w')
ax.zaxis.pane.set_edgecolor('w')

for d in dirs_compare:
    data = np.loadtxt(os.path.join(d, "lorenz.txt"))
    t, x, y, z = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    ax.plot(x, y, z, linewidth=0.5, label=ode_lib[d])

plt.legend()
plt.show()
