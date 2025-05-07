import numpy as np
import os

from scipy.integrate import solve_ivp

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

a = 1.0
n_x = 500

dx = 1./n_x
dt = dx / a

fig, axs = plt.subplots(len(ode_lib), 1, sharex='col', figsize=(7.85, 14))

for i, d in enumerate(ode_lib.keys()):
    print(f"extract data from: {d}")
    data = np.loadtxt(os.path.join(d, "transport.txt"))
    # axs[i].imshow(np.transpose(data[:, 1:]))
    # axs[i].set_ylabel("time")

    # axs[i].margins(0.05)

    x = data[:, 0]
    y = data[:, 1:]
    n_sol = y.shape[1]-1

    axs[i].set_xlim(x[0], x[-1])
    axs[i].text(x[5], 0.9*np.max(y), ode_lib[d], horizontalalignment="left")

    for n in range(10):
        idx = int(n_sol*n/9)
        yn = y[:, idx]
        axs[i].plot(x, yn, label=f"$t^n = {idx * dt}$")

axs[-1].set_xlabel("$x$")
axs[0].legend(ncols=5, loc=(0., 1.05))
# axs[-1].legend(ncols=5, loc=(0., -0.45))

fig.tight_layout()
plt.savefig(os.path.join("transport.png"))

# Arenstorf orbit
print("arenstorf orbit")


def arenstorf_system(t, y):
    mu = 0.012277471

    y1 = y[0]
    y2 = y[1]
    y3 = y[2]
    y4 = y[3]

    r1 = np.sqrt((y1 + mu) * (y1 + mu) + y2 * y2)
    r2 = np.sqrt((y1 - 1. + mu) * (y1 - 1. + mu) + y2 * y2)

    dy1 = y3
    dy2 = y4
    dy3 = y1 + 2 * y4 - (1 - mu) * (y1 + mu) / \
        (r1 * r1 * r1) - mu * (y1 - 1 + mu) / (r2 * r2 * r2)
    dy4 = y2 - 2 * y3 - (1 - mu) * y2 / (r1 * r1 * r1) - \
        mu * y2 / (r2 * r2 * r2)

    return np.asarray([dy1, dy2, dy3, dy4])


y0 = np.array([0.994, 0., 0., -2.00158510637908252240537862224])
t_span = [0., 17.0652165601579625588917206249]
sol = solve_ivp(arenstorf_system, t_span, y0, method="RK45",
                first_step=1e-8, rtol=1e-10, atol=1e-10)

fig, ax = plt.subplots(figsize=(7, 7))
ax.set_aspect('equal', 'box')
ax.grid(visible=True)

ax.plot(sol.y[0, :], sol.y[1, :], "--", color="black")


ode_lib_arenstorf_files = {
    "ascent": [{'method': "RK5(4) 7M", 'file': "arenstorf.txt"}],
    "gsl": [{'method': "RK8(7) 13M", 'file': "arenstorf.txt"}],
    "odeint": [{'method': "RK5(4) 7M", 'file': "arenstorf.txt"}],
    "petsc": [{'method': "RK5(4) 7M", 'file': "arenstorf.txt"}],
    "ponio": [{'method': "RK5(4) 7M", 'file': "arenstorf_rk54_7m.txt"}, {'method': "RK8(7) 13M", 'file': "arenstorf_rk87_13m.txt"}],
    "scipy": [{'method': "RK5(4) 7M", 'file': "arenstorf_rk45.txt"}, {'method': "RK8(7) 13M", 'file': "arenstorf_dop853.txt"}],
}

for i, d in enumerate(ode_lib.keys()):
    print(f"extract data from: {d}")
    for simu in ode_lib_arenstorf_files[d]:
        data = np.loadtxt(os.path.join(d, simu['file']))
        ax.plot(data[:, 1], data[:, 2], "x",
                label=f"{ode_lib[d]} {simu['method']}")

ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
fig.legend()
plt.savefig(os.path.join("arenstorf.png"))
