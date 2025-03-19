import numpy as np
import os

import matplotlib.pyplot as plt
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

colors = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
]

# Lorenz attractor
print("lorenz figure")

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
    scene_camera=dict(
        eye=dict(x=1.4, y=-1.1, z=0.6),
    ),
    autosize=False,
    width=600,
    height=600,
    template="plotly_white"
)
fig.write_html("lorenz.html")


# transport
print("transport figure")

fig, axs = plt.subplots(len(ode_lib), 1, sharex='col', figsize=(8.3, 14))

for i, d in enumerate(ode_lib.keys()):
    print(f"extract data from: {d}")
    data = np.loadtxt(os.path.join(d, "transport.txt"))
    # axs[i].imshow(np.transpose(data[:, 1:]))
    # axs[i].set_ylabel("time")

    # axs[i].margins(0.05)

    axs[i].text(1.0, 0.2, ode_lib[d], horizontalalignment="right")
    x = data[:, 0]
    y = data[:, 1:]
    n_sol = y.shape[1]

    for n in range(10):
        yn = y[:, int(n/10*n_sol)]
        axs[i].plot(x, yn, label=f"iteration {int(n/10*n_sol)}")

axs[-1].set_xlabel("$x$")
axs[0].legend(ncols=5, loc=(0., 1.05))
# axs[-1].legend(ncols=5, loc=(0., -0.45))

fig.tight_layout()
plt.savefig(os.path.join("transport.png"))
