# Intervals for boundary conditions
import numpy as np
n_whole = 905263.0

p = 0.7
i_s_0 = 25.0 / n_whole
i_a_0 = 75.0 / n_whole
l_0 = p * (n_whole -(i_s_0 + i_a_0)) / n_whole
s_0 = (1.0 - p) * (n_whole -(i_s_0 + i_a_0)) / n_whole


sigma = 1e-6
x_0 = np.array([l_0, s_0, i_s_0, i_a_0], dtype=float)
interval_length = sigma * np.abs(sigma * np.random.randn())
x_0_l =  x_0 - 0.5 * interval_length * x_0
x_0_r =  x_0 + 0.5 * interval_length * x_0
print(x_0_l)
print(x_0_r)
