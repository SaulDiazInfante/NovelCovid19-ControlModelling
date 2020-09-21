import numpy as np
import json
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime


def load_parameters(file_name='exit_lockdown_parameters.json'):
    """ Load the parameters given a JSON file.
    Parameters
    ----------
    file_name: filename with parameters values in JSON format.

    Returns
    -------
    prm: dictionary
    """
    with open(file_name) as json_file:
        prm = json.load(json_file)
    return prm


def data_solution_save(sol, column_names, file_name_prefix='lock_down_data_solution'):
    """ Save scipy.integrate.solve_ivp output as pandas DataFrame
    Parameters
    ----------
    sol: dict_values.
    column_names: Names for the columns of the DataFrame.
    file_name_prefix: filename  of desire DataFrame file.

    Returns
    -------
    df.plk: a pkl pandas DataFrame file.
    """
    df = pd.DataFrame(sol)
    df.columns = column_names
    string_date = str(datetime.date(datetime.now()))
    file_name = file_name_prefix + string_date + ".pkl"
    df.to_pickle(file_name)
    return


def lockdown_population_plot(r_zero, data_file_name='lock_down_data_solution2020-05-29.pkl',
                                         fig_file_name='model_solution.pdf'):
    """
    Plot the output solution given the pandas data frame
    output from scipy.integrate.solve_ivp
    ----------
    file_name: filename with pandas data frame  output from scipy.integrate.solve_ivp

    Returns
    -------
    Shows figures in the default backend and return a pdf figure file
    """
    df = pd.read_pickle(data_file_name)
    fig, ((ax_l, ax_s, ax_e), (ax_i_s, ax_i_a, ax_h), (ax_r, ax_d, ax_cl)) = plt.subplots(nrows=3, ncols=3)
    #
    n_whole = 905263
    t = df['time']
    ax_l = plt.subplot(331)
    ax_l.plot(t, n_whole * df['l'], label="l")
    ax_l.legend(loc=0)
    #
    ax_s = plt.subplot(332)
    ax_s.plot(t, n_whole * df['s'], label="s")
    ax_s.legend(loc=0)
    #
    ax_e = plt.subplot(333)
    ax_e.plot(t, n_whole * df['e'], label="e")
    ax_e.legend(loc=0)
    #
    ax_i_s = plt.subplot(334)
    ax_i_s.plot(t, n_whole * df['i_s'], label="i_s")
    ax_i_s.legend(loc=0)
    #
    ax_i_a = plt.subplot(335)
    ax_i_a.plot(t, n_whole * df['i_a'], label="i_a")
    ax_i_a.legend(loc=0)
    #
    ax_h = plt.subplot(336)
    ax_h.plot(t, n_whole * df['h'], label="h")
    ax_h.legend(loc=0)
    #
    ax_r = plt.subplot(337)
    ax_r.plot(t, n_whole * df['r'], label="r")
    ax_r.legend(loc=0)
    #
    ax_d = plt.subplot(338)
    ax_d.plot(t, n_whole * df['d'], label="d")
    ax_d.legend(loc=0)
    #
    cl = n_whole * (df['l'] + df['s'] + df['e'] + df['i_s'] + df['i_a'] + df['h'] + df['r'] + df['d'])
    ax_cl = plt.subplot(339)
    ax_cl.plot(t, cl, label="cl")
    ax_cl.legend(loc=0)
    #
    plt.tight_layout()
    fig.suptitle("R0: " + str(r_zero))
    plt.savefig(fig_file_name)
    # plt.show()
    return


def reproductive_number(l_0, s_0, beta_s, beta_a, p, alpha_a, alpha_s):
    r_0 = p * (p * beta_s / alpha_s + (1.0 - p) * beta_a / alpha_a) * (l_0 + s_0)
    return r_0


def rhs_lockdown_exit_00(y, t, beta_s, beta_a, epsilon, delta_u, delta_e, p, alpha_a,
                         alpha_h, alpha_s, mu_h, a_1, a_2, a_3):
    """"
        Version designed in May 10-2020, according to the normalized version and
        lockdown process of the above Showing results for Kermack–McKendrick Model.
        The main idea is to split symptomatic
        individuals into Hospitalized, Intensive Care and mild.
        @param a_3:
        @param a_2:
        @param a_1:
        @param alpha_a:
        @param p:
        @param delta_e:
        @param delta_u:
        @param epsilon:
        @param beta_a:
        @param beta_s:
        @param y:
        @param mu_h:
        @param alpha_s:
        @param alpha_h:
        @type t: float dummy
    """

    l, s, e, i_s, i_a, h, r, d = y
    #
    n_aster = l + s + e + i_s + i_a + h + r
    force_infection = (beta_s * i_s + beta_a * i_a) / n_aster
    rhs_l = - l * (epsilon * force_infection + delta_u)
    rhs_s = delta_u * l - force_infection * s
    rhs_e = force_infection * (epsilon * l + s) - delta_e * e
    rhs_i_s = p * delta_e * e - alpha_s * i_s
    rhs_i_a = (1 - p) * delta_e * e - alpha_a * i_a
    rhs_h = a_1 * alpha_s * i_s - (alpha_h + mu_h) * h
    rhs_r = a_2 * alpha_s * i_s + alpha_a * i_a + alpha_h * h
    rhs_d = a_3 * alpha_s * i_s + mu_h * h
    rhs = np.array([rhs_l, rhs_s, rhs_e, rhs_i_s, rhs_i_a, rhs_h, rhs_r, rhs_d])
    return rhs


def rhs_lockdown_exit_01(t, y, beta_s, beta_a, epsilon,
                         delta_u, delta_e, delta_s, p,
                         alpha_a, alpha_h, alpha_s,
                         mu_s, mu_h):
    """"
        Version designed in May 22-2020, according to the normalized version and
        lockdown process of the above Kermack Model.
        @param mu_h:
        @param mu_s:
        @param alpha_s:
        @param alpha_h:
        @param alpha_a:
        @param p:
        @param delta_s:
        @param delta_e:
        @param delta_u:
        @param epsilon:
        @param beta_a:
        @param beta_s:
        @param y:
        @type t: object
    """
    l, s, e, i_s, i_a, h, r, d = y
    #
    n_aster = l + s + e + i_s + i_a + h + r
    force_infection = (beta_s * i_s + beta_a * i_a) / n_aster
    rhs_l = - l * (epsilon * force_infection + delta_u)
    rhs_s = delta_u * l - force_infection * s
    rhs_e = force_infection * (epsilon * l + s) - delta_e * e
    rhs_i_s = p * delta_e * e - (alpha_s + mu_s + delta_s) * i_s
    rhs_i_a = (1 - p) * delta_e * e - alpha_a * i_a
    rhs_h = delta_s * i_s - (alpha_h + mu_h) * h
    rhs_r = alpha_s * i_s + alpha_a * i_a + alpha_h * h
    rhs_d = mu_s * i_s + mu_h * h
    rhs = np.array([rhs_l, rhs_s, rhs_e, rhs_i_s, rhs_i_a, rhs_h, rhs_r, rhs_d])
    return rhs

# TODO: Modify beta parameters to estimate the exit lockdown
