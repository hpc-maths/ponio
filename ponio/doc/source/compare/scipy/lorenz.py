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


def lorenz_system(t, y):
    sigma = 10.
    rho = 28.
    beta = 8./3.

    return np.asarray([
        sigma * (y[1] - y[0]),
        y[0]*(rho - y[2]) - y[1],
        y[0] * y[1] - beta * y[2]
    ])


y0 = np.array([1., 1., 1.])
t_span = [0., 20.]
dt = 0.01

sol = solve_ivp(lorenz_system, t_span, y0, method=RK44,
                first_step=dt, max_step=dt, rtol=1.0, atol=1.0)

save_sol = np.zeros(shape=(4, len(sol.t)))
save_sol[0, :] = sol.t
save_sol[1:, :] = sol.y

np.savetxt("lorenz.txt", np.transpose(save_sol), delimiter=" ")
