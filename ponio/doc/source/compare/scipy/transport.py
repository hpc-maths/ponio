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


n_x = 100
t_span = [0., 0.3]
dt = 0.01

x = np.linspace(0., 1., n_x, endpoint=False)
dx = x[1] - x[0]

a = 1.0


def transport_system(t, y):
    dy = np.zeros_like(y)

    dy[0] = -a * (y[1] - y[-1]) / (2.*dx)
    dy[1:-1] = -a * (y[2:] - y[:-2]) / (2.*dx)
    dy[-1] = -a * (y[0] - y[-2]) / (2.*dx)

    return dy


# y0 = np.sin(2.*np.pi*x)
sigma = 0.1
# y0 = 1./(sigma*np.sqrt(2.*np.pi))*np.exp(-(x-0.5)**2/sigma**2)
y0 = np.zeros_like(x)
y0 = np.array([xi - 0.25 if 0.25 <= xi and xi < 0.5 else -xi +
              0.75 if 0.5 <= xi and xi < 0.75 else 0. for xi in x])

sol = solve_ivp(transport_system, t_span, y0, method=RK44,
                first_step=dt, max_step=dt, rtol=1.0, atol=1.0)

save_sol = np.zeros(shape=(len(sol.t)+1, n_x))
save_sol[0, :] = x
save_sol[1:, :] = np.transpose(sol.y)

np.savetxt("transport.txt", np.transpose(save_sol), delimiter=" ")
