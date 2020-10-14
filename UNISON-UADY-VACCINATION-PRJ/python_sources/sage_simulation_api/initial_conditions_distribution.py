import numpy as np
import pandas as pd


def initial_distribution_pop_m3(n_whole, i_0, r_0, d_0):
    # observables
    rhs_i_0 = i_0 / n_whole
    rhs_r_0 = r_0 / n_whole
    rhs_d_0 = d_0 / n_whole
    # partially observable
    regularization_factor = np.random.uniform(8.2, 30)
    x_0 = regularization_factor * i_0
    rhs_x_0 = x_0 / n_whole
    observable_pop = i_0 + r_0 + d_0 + x_0
    partial_observable_pop = n_whole - observable_pop
    # TODO: fix the initial random partition in order that the returned
    # distribution sums one

    p = np.random.uniform(.01, .07)
    # l_0 = p * partial_observable_pop
    s_0 = (9999.0 / 10000) * (1.0 - p) * partial_observable_pop
    r_x_0 = (1.0 / 10000) * (1.0 - p) * partial_observable_pop
    e_0 = partial_observable_pop - (s_0 + r_x_0)
    #
    rhs_s_0 = s_0 / n_whole
    rhs_r_0 = (r_0 + r_x_0) / n_whole
    rhs_e_0 = e_0 / n_whole

    dist = np.array([rhs_s_0,
                     rhs_e_0,
                     rhs_i_0,
                     rhs_x_0,
                     rhs_r_0,
                     rhs_d_0,
                     0.0,
                     0.0,
                     0])
    return dist


def initial_conditions_intervals(distribution):
    interval_lengths = np.random.uniform(1e-12, 2e-12, distribution.shape[0])
    lower_bound_intervals = distribution - interval_lengths * distribution
    upper_bound_intervals = distribution + interval_lengths * distribution
    data = np.transpose(np.array([lower_bound_intervals,
                                  distribution,
                                  upper_bound_intervals]))
    df_intervals = pd.DataFrame(data, columns=['lb', 'mean', 'ub'])
    return df_intervals


# Culiacan setting August 17 2020
n_pop = 962871
i_s_0 = 217
ics = 5307
deaths = 869
recovered = np.ceil(.3 * n_pop)
dist = initial_distribution_pop_m3(n_pop, i_s_0, recovered, deaths)
intervals = initial_conditions_intervals(dist)
intervals.to_csv('intervals.csv')
