from initial_conditions_distribution import *
n_pop = 21581000.0
i_0 = 20000
r_0 = 13500
d_0 = 1500
#
beta_0 = 0.91
alpha = 0.25
gamma = 1.0 / 5.0
rho = 1.0 / 4
mu = 1.0 / (70.0 * 365.0)
p = 0.25
T = 200
initial_partition = initial_distribution_pop_m3(n_pop, i_0, r_0, d_0)
print(initial_partition.sum())
