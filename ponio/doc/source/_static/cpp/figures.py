import numpy as np
import matplotlib.pyplot as plt

### Lorenz ##########################################################


def plot_lorenz(filename: str, output: str):
    data = np.loadtxt(filename)

    fig = plt.figure(constrained_layout=True, figsize=(8, 8))
    gs = fig.add_gridspec(6, 1)

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

    ax.plot(data[:, 1], data[:, 2], data[:, 3])

    plt.savefig(output)


def plot_multi_lorenz(methods: dict[str, str], output: str):
    fig = plt.figure(constrained_layout=True, figsize=(8, 8))
    gs = fig.add_gridspec(6, 1)

    ax_3d = fig.add_subplot(gs[:3, 0], projection='3d')
    ax_x = fig.add_subplot(gs[3, 0])
    ax_y = fig.add_subplot(gs[4, 0])
    ax_z = fig.add_subplot(gs[5, 0])

    for i, (meth, filename) in enumerate(methods.items()):
        data = np.loadtxt(filename)
        t, x, y, z = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
        ax_3d.plot(x, y, z, label=meth, color=f"C{i}", linewidth=0.375)
        ax_x.plot(t, x, color=f"C{i}")
        ax_y.plot(t, y, color=f"C{i}")
        ax_z.plot(t, z, color=f"C{i}")

    ax_3d.set_xlim3d([-20, 20])
    ax_3d.set_ylim3d([-20, 20])
    ax_3d.set_zlim3d(bottom=0, top=50)

    ax_x.set_ylabel("$x$")
    ax_x.set_ylim(-20, 20)

    ax_y.set_ylabel("$y$")
    ax_y.set_ylim(-20, 20)

    ax_z.set_ylabel("$z$")
    ax_z.set_xlabel("$t$")
    ax_z.set_ylim(0, 50)

    ax_list = [ax_3d, ax_x, ax_y, ax_z]
    handles, labels = map(lambda x: sum(x, start=[]), zip(
        *map(lambda ax: ax.get_legend_handles_labels(), ax_list)))
    fig.legend(handles, labels, loc="upper left",
               bbox_to_anchor=(0.01, 0.99), ncol=1)

    plt.savefig(output)


def lorenz_fig():
    # explicit RK
    plot_lorenz("lorenz_rk.txt", "lorenz_rk.png")

    # 4 methods
    methods = {
        'RK(4, 4) 3/8': "lorenz_rk.txt",
        'Lawson': "lorenz_lrk.txt",
        'DIRK': "lorenz_dirk.txt",
        'Strang': "lorenz_split.txt"
    }
    plot_multi_lorenz(methods, "lorenz.png")

### Lotka-Voleterra #################################################


def plot_lv_obs(filename: str, output: str):
    data = np.loadtxt(filename)

    x, y = data[:, 1], data[:, 2]

    # plot x-y
    ax = plt.figure().add_subplot()
    ax.plot(x, y)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    plt.savefig(output)


def plot_lv_V(methods: dict[str, str], output: str, compute_V=list[bool]):
    ax = plt.figure().add_subplot()
    for i, ((meth, filename), cV) in enumerate(zip(methods.items(), compute_V)):
        data = np.loadtxt(filename)
        t, x, y = data[:, 0], data[:, 1], data[:, 2]

        if cV:
            alpha, beta, delta, gamma = 2./3., 4./3., 1., 1.
            V = delta*x - gamma*np.log(x) + beta*y - alpha*np.log(y)
        else:
            V = data[:, 3]

        ax.plot(t, np.abs(V/V[0] - 1.), color=f"C{i}", label=meth)

    ax.legend()
    plt.savefig(output)


def lotka_volterra_fig():
    # each method
    plot_lv_obs("lotka_volterra_fobs.txt", "lotka_volterra_fobs.png")
    plot_lv_obs("lotka_volterra_cobs.txt", "lotka_volterra_cobs.png")
    plot_lv_obs("lotka_volterra_sobs.txt", "lotka_volterra_sobs.png")
    plot_lv_obs("lotka_volterra_uobs.txt", "lotka_volterra_uobs.png")

    methods = {
        "file observer": "lotka_volterra_fobs.txt",
        "cout observer": "lotka_volterra_cobs.txt",
        "stream observer": "lotka_volterra_sobs.txt",
        "user observer": "lotka_volterra_uobs.txt"
    }
    plot_lv_V(methods, "lotka_volterra.png", [True, True, True, False])

### Curtiss-Hirschfelder ############################################


def plot_ch(filename: str, output: str):
    data = np.loadtxt(filename)
    t, y, dt = data[:, 0], data[:, 1], data[:, 2]

    fig, axs = plt.subplots(2)
    axs[0].plot(t, y)
    axs[1].plot(t, dt)

    axs[0].set_ylabel("$y$")
    axs[1].set_ylabel("$\\Delta t$")
    axs[1].set_xlabel('time')

    plt.savefig(output)


def curtiss_hirschfelder_fig():
    plot_ch("curtiss_hirschfelder_solve.txt", "curtiss_hirschfelder_solve.png")
    plot_ch("curtiss_hirschfelder_while.txt", "curtiss_hirschfelder_while.png")
    plot_ch("curtiss_hirschfelder_for.txt", "curtiss_hirschfelder_for.png")


lorenz_fig()
lotka_volterra_fig()
curtiss_hirschfelder_fig()
