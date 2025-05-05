import numpy as np
import scipy
from scipy.integrate import solve_ivp


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
dt = 1e-5

sol = solve_ivp(arenstorf_system, t_span, y0, method="DOP853",
                first_step=dt, rtol=1e-5, atol=1e-5)

save_sol = np.zeros(shape=(5, len(sol.t)))
save_sol[0, :] = sol.t
save_sol[1:, :] = sol.y

np.savetxt("arenstorf.txt", np.transpose(save_sol), delimiter=" ")
