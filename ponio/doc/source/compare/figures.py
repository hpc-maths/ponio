import numpy as np
import os

# import matplotlib.pyplot as plt
import plotly.graph_objects as go

dirs_compare = next(os.walk("."))[1]

ode_lib = {
    "ascent": "Ascent",
    # "diffeq": "DifferentialEquations.jl",
    "gsl": "GSL",
    "odeint": "Boost::odeint",
    "petsc": "PETSc",
    "ponio": "ponio",
    "scipy": "SciPy"
}


# matplotlib.pyplot version

# fig = plt.figure(constrained_layout=True, figsize=(8, 8))
# gs = fig.add_gridspec(6, 1)
# axlist = []

# ax = fig.add_subplot(gs[:, 0], projection='3d')
# ax.set_xlim3d([-20, 20])
# ax.set_ylim3d([-20, 20])
# ax.set_zlim3d(bottom=0, top=50)

# ax.xaxis.pane.fill = False
# ax.yaxis.pane.fill = False
# ax.zaxis.pane.fill = False
# ax.xaxis.pane.set_edgecolor('w')
# ax.yaxis.pane.set_edgecolor('w')
# ax.zaxis.pane.set_edgecolor('w')

# # for d in dirs_compare:
# for d in ode_lib.keys():
#     print(f"extract data from: {d}")
#     data = np.loadtxt(os.path.join(d, "lorenz.txt"))
#     t, x, y, z = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
#     ax.plot(x, y, z, linewidth=0.5, label=ode_lib[d])

# plt.legend()
# plt.savefig(os.path.join("lorenz.png"))


# plotly version
colors = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
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
    scene=dict(aspectmode='cube'),
    autosize=False,
    width=600,
    height=600,
    template="plotly_white"
)
fig.write_html("lorenz.html")
