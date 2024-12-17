import numpy as np
from scipy.integrate import solve_ivp


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

sol = solve_ivp(lorenz_system, t_span, y0)

save_sol = np.zeros(shape=(4, len(sol.t)))
save_sol[0, :] = sol.t
save_sol[1:, :] = sol.y

np.savetxt("lorenz.txt", np.transpose(save_sol), delimiter=" ")
