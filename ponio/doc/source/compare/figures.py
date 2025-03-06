import numpy as np
import matplotlib.pyplot as plt
import os
import plotly.graph_objects as go

dirs_compare = next(os.walk("."))[1]

ode_lib = {
    "ascent": "Ascent",
    # "diffeq": "DifferentialEquations.jl",
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

# for d in dirs_compare:
for d in ode_lib.keys():
    print(f"extract data from: {d}")
    data = np.loadtxt(os.path.join(d, "lorenz.txt"))
    t, x, y, z = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    ax.plot(x, y, z, linewidth=0.5, label=ode_lib[d])

plt.legend()
plt.savefig(os.path.join("lorenz.png"))

# plotly version
colors = [
    "#FEA47F", "#25CCF7", "#EAB543", "#55E6C1", "#B33771", "#3B3B98", "#FD7272", "#9AECDB", "#D6A2E8"
]

fig = go.Figure()
for i, d in enumerate(ode_lib.keys()):
    print(f"extract data from: {d}")
    data = np.loadtxt(os.path.join(d, "lorenz.txt"))
    fig.add_trace(go.Scatter3d(
        x=data[:, 1], y=data[:, 2], z=data[:, 3],
        marker=dict(size=1),
        line=dict(
            color=colors[i],
            width=2
        ),
        name=ode_lib[d]
    ))

fig.update_layout(
    autosize=False,
    width=800,
    height=800
)
fig.write_html("lorenz.html")
