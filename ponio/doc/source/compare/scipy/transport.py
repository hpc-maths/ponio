import numpy as np
import scipy
from scipy.integrate import solve_ivp


class RK44(scipy.integrate._ivp.rk.RungeKutta):
    order = 4
    error_estimator_order = 1
    n_stages = 4
    C = np.array([0, 1/2, 1/2, 1])
    A = np.array([
        [0, 0, 0, 0],
        [1/2, 0, 0, 0],
        [0, 1/2, 0, 0],
        [0, 0, 1, 0]
    ])
    B = np.array([1/6, 1/3, 1/3, 1/6])
    E = np.array([0.1, 0.1, 0.1, 0.1, 0.1])
    P = np.array([[1, -4 / 3, 5 / 9, 0],
                  [0, 1, -2/3, 0],
                  [0, 4/3, -8/9, 0],
                  [0, -1, 1, 0]])


# space parameters
n_x = 500
x = np.linspace(0., 1., n_x, endpoint=False)
dx = x[1] - x[0]

# velocity
a = 1.0

# time parameter
t_span = [0., 0.3]
dt = dx / a

# initial condition
y0 = np.array([xi - 0.25 if 0.25 <= xi and xi < 0.5 else -xi +
              0.75 if 0.5 <= xi and xi < 0.75 else 0. for xi in x])


def upwind(t, y):
    dy = np.zeros_like(y)

    dy[0] = - (np.max([a, 0]) * (y[0] - y[-1]) +
               np.min([a, 0]) * (y[1] - y[0])) / dx
    dy[1:-1] = -(np.max([a, 0]) * (y[1:-1] - y[:-2]) +
                 np.min([a, 0]) * (y[2:] - y[1:-1])) / dx
    dy[-1] = - (np.max([a, 0]) * (y[-1] - y[-2]) +
                np.min([a, 0]) * (y[0] - y[-1])) / dx

    return dy


sol = solve_ivp(upwind, t_span, y0, method=RK44,
                first_step=dt, max_step=dt, rtol=1.0, atol=1.0)

# save solution
save_sol = np.zeros(shape=(len(sol.t)+1, n_x))
save_sol[0, :] = x
save_sol[1:, :] = np.transpose(sol.y)

np.savetxt("transport.txt", np.transpose(save_sol), delimiter=" ")
