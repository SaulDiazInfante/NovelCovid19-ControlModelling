import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
from sixr_dynamics import rhs_sixr_dynamics, rhs_lsixrd_dynamics
from initial_conditions_distribution import initial_distribution_pop_m3

n_pop = 21581000.0
i_0 = 1866
r_0 = 3331
d_0 = 387
#
beta_0 = 0.91
alpha = 0.25
gamma = 1.0 / 5.0
rho = 1.0 / 4
mu = 1.0 / (70.0 * 365.0)
delta = 0.097
p = 0.25
T = 500
z_0 = initial_distribution_pop_m3(n_pop, i_0, r_0, d_0)
# TODO: fix the whole population notation and implements the new initial distribution for each model
parameters_dictionary = {
    "N": n_pop,
    "y_0": z_0,
    "s_0": .5 * n_pop,
    "i_0": .01 * n_pop,
    "x_0": .08 * n_pop,
    "r_0": .005 * n_pop,
    "r_x_0": .005 * n_pop,
    "mu": 1.0 / (70.0 * 365.0),
    "beta_0": 0.91,
    "alpha": 0.25,
    "gamma": 1.0 / 20.0,
    "rho": 1.0/14.0,
    "p": 0.25,
    "T": 180
}
# z_0 = np.array([z_0, s_0, i_0, x_0, r_0, r_x_0])

t_eval = np.linspace(0, T, 100000)
# sol = solve_ivp(rhs_lsixrd_dynamics, [0, T], z_0,
#                dense_output=False,
#                method='RK45',
#                t_eval=t_eval,
#                events=None,
#                vectorized=True,
#                args=(beta_0, alpha, gamma, rho, p, delta)
#                )
sol = odeint(rhs_lockdown, y0, t, args=(b, c))

t = sol.t
y = sol.y.T
#
l = y[:, 1]
s = y[:, 2]
i = y[:, 3]
x = y[:, 4]
r = y[:, 5]
r_x = y[:, 6]
d = y[:, 7]

fig, ((ax_l, ax_s, ax_i, ax_x), (ax_r, ax_r_x, ax_d, ax_0)) = plt.subplots(nrows=2, ncols=4)

ax_l = plt.subplot(241)
ax_l.plot(t, l, label="l")
plt.legend(loc=0)

ax_s = plt.subplot(242)
ax_s.plot(t, s, label="s")
plt.legend(loc=0)

ax_i = plt.subplot(243)
ax_i.plot(t, i, label="i")
plt.legend(loc=0)

ax_x = plt.subplot(244)
ax_x.plot(t, x, label="x")
plt.legend(loc=0)

ax_r = plt.subplot(245)
ax_r.plot(t, r, label="r")
plt.legend(loc=0)

ax_r_x = plt.subplot(246)
ax_r_x.plot(t, r_x, label="r_x")
plt.legend(loc=0)

ax_d = plt.subplot(247)
ax_d.plot(t, d, label="d")
plt.legend(loc=0)
plt.tight_layout()
plt.figure()
plt.plot(t, y)
plt.show()
