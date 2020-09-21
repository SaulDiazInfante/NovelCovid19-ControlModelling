import numpy as np
import json
import matplotlib.pyplot as plt
import pandas as pd
# from datetime import datetime


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


def data_solution_save(sol, column_names,
                       file_name_prefix='vaccination_data_solution'):
    """ 
    Save scipy.integrate.solve_ivp output as pandas DataFrame
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
    # string_date = str(datetime.date(datetime.now()))
    file_name = file_name_prefix + ".pkl"  # string_date + ".pkl"
    df.to_pickle(file_name)
    return


def vaccination_population_plot(r_zero,
                                data_file_name='vaccination_data_solution.pkl',
                                fig_file_name='model_solution.pdf'):
    """
    Plot the output solution given the pandas data frame
    output from scipy.integrate.solve_ivp
    ----------
    file_name: filename with pandas data frame  output from
        scipy.integrate.solve_ivp

    Returns
    -------
    Shows figures in the default backend and return a pdf figure file
    """
    # plt.ion()
    print(data_file_name)
    df_not = pd.read_pickle("not_vaccination_solution.pkl")
    data_file_name = '/home/saul/Insync/saul.diazinfante@unison.mx' + \
                     '/OneDrive Biz/UNISON/ARTICLES/Covid19/' + \
                     'COVID19-Sonora/python_code/first_soltuions_set/' + \
                     'lambda_v-9.736842e-04.sol.csv'
    print(data_file_name)
    df_constant = pd.read_pickle("vaccination_solution.pkl")
    df_oc = pd.read_csv(data_file_name, header=0)
    plt.figure(0)
    # plt.ion()
    ax_s = plt.subplot2grid((4, 3), (0, 0))
    ax_e = plt.subplot2grid((4, 3), (0, 1))
    ax_i_s = plt.subplot2grid((4, 3), (0, 2))
    #
    ax_i_a = plt.subplot2grid((4, 3), (1, 0))
    ax_r = plt.subplot2grid((4, 3), (1, 1))
    ax_d = plt.subplot2grid((4, 3), (1, 2))
    #
    ax_v = plt.subplot2grid((4, 3), (2, 0))
    ax_vc = plt.subplot2grid((4, 3), (2, 1))
    ax_hz = plt.subplot2grid((4, 3), (2, 2))
    #
    ax_cl = plt.subplot2grid((4, 3), (3, 0))
    ax_uv = plt.subplot2grid((4, 3), (3, 1))
    ax_oc_v = plt.subplot2grid((4, 3), (3, 2))
    #

    prm = load_parameters("model_vaccination_parameters.json")
    n_cdmx_vm = prm['n_whole']
    # n_cul = 905263
    n_whole = 1.0  # without rescaling population
    t = df_constant['time']
    s = df_constant['s']
    e = df_constant['e']
    i_s = df_constant['i_s']
    i_a = df_constant['i_a']
    r = df_constant['r']
    d = df_constant['d']
    v = df_constant['v']
    x_vaccine = n_whole * (df_constant['x_vac_counter'])
    if data_file_name == 'vaccination_treatment_solution.pkl':
        treat = df_constant['treat']
        print('treatment')
        cl = n_whole * (s + e + i_s + i_a + r + d + v + treat)
    else:
        treat = np.zeros(df_constant.values.shape[0])
        cl = n_whole * (s + e + i_s + i_a + r + d + v)
    #
    ax_s.plot(t, n_whole * df_not['s'], 'k--')
    ax_s.plot(t, n_whole * s, alpha=0.7)
    ax_s.plot(df_oc['time'], n_whole * df_oc['s'],
              ls=':',
              lw=2,
              color='orange')
    ax_s.legend(loc=0)
    ax_s.title.set_text('Susceptible')
    #
    ax_e.plot(t, n_whole * df_not['e'], 'k--')
    ax_e.plot(t, n_whole * e, alpha=0.7)
    ax_e.plot(df_oc['time'], n_whole * df_oc['e'],
              ls=':',
              lw=2,
              color='orange')
    ax_e.legend(loc=0)
    ax_e.title.set_text('Exposed')
    #
    ax_i_s.plot(t, n_whole * df_not['i_s'], 'k--')
    ax_i_s.plot(t, n_whole * i_s, label="CC: i_s", alpha=0.7)
    ax_i_s.plot(df_oc['time'], n_whole * df_oc['i_s'],
                ls=':',
                lw=2,
                color='orange')
    ax_i_s.legend(loc=0)
    ax_i_s.title.set_text('Symptomatic')
    #
    ax_i_a.plot(t, n_whole * df_not['i_a'], 'k--')
    ax_i_a.plot(t, n_whole * i_a, alpha=0.7)
    ax_i_a.plot(df_oc['time'], n_whole * df_oc['i_a'],
                ls=':',
                lw=2,
                color='orange')
    ax_i_a.legend(loc=0)
    ax_i_a.title.set_text('Asymptomatic')
    #
    ax_r.plot(t, n_whole * df_not['r'], 'k--')
    ax_r.plot(t, n_whole * r, alpha=0.7)
    ax_r.plot(df_oc['time'], n_whole * df_oc['r'],
              ls=':',
              lw=2,
              color='orange')
    ax_r.legend(loc=0)
    ax_r.title.set_text('Recovered')
    #
    ax_d.plot(t, n_whole * df_not['d'], 'k--')
    ax_d.plot(t, n_whole * d, label="CC: d", alpha=0.7)
    ax_d.plot(df_oc['time'], n_whole * df_oc['d'],
              ls=':',
              lw=2,
              color='orange')
    ax_d.fill_between(t,
                      # n_whole * df_oc['d'],
                      n_whole * d,
                      n_whole * df_not['d'],
                      facecolor='green',
                      interpolate=True,
                      alpha=0.5)
    ax_d.legend(loc=0)
    ax_d.title.set_text('Deaths')
    #
    ax_v.plot(t, n_whole * df_not['v'], 'k--')
    ax_v.plot(t, n_whole * v, alpha=0.7)
    ax_v.plot(df_oc['time'], n_whole * df_oc['v'],
              ls=':',
              lw=2,
              color='orange')
    ax_v.legend(loc=0)
    ax_v.title.set_text('Effective Vaccinated')
    #
    ax_hz.plot(t, n_cdmx_vm * 0.25 * 0.25 * df_not['i_s'], '--k')
    ax_hz.plot(t, n_cdmx_vm * 0.25 * 0.25 * i_s)
    ax_hz.plot(df_oc['time'], n_cdmx_vm * 0.016 * df_oc['i_s'],
               ls=':',
               lw=2,
               color='orange')

# https://www.forbes.com.mx/noticias-capacidad-hospitalaria-zona-metropolitana/
    ax_hz.hlines(8500, -20, t.values[-1],
                 colors='gray',
                 linestyles=':',
                 alpha=0.5,
                 label='ICU')
    #
    ax_hz.hlines(1800, -20, t.values[-1],
                 colors='red',
                 linestyles=':',
                 alpha=0.5,
                 label='ICU')
    ax_hz.title.set_text('Hospitalised')
    #
    ax_cl.plot(t, cl, label="cl")
    ax_cl.legend(loc=0)
    #
    ax_vc.plot(t, x_vaccine)
    ax_vc.hlines(0.8 * n_whole, t.values[0], t.values[-1],
                 colors='gray',
                 linestyle=':',
                 alpha=0.5)
    ax_vc.hlines(0.5 * n_whole, t.values[0], t.values[-1],
                 colors='gray',
                 linestyle='--',
                 alpha=0.5)
    ax_vc.hlines(0.2 * n_whole, t.values[0], t.values[-1],
                 colors='gray',
                 linestyle='-.',
                 alpha=0.5)
    ax_vc.plot(df_oc['time'], n_whole * df_oc['x_vac_counter'],
               ls=':',
               lw=2,
               color='orange')
    ax_vc.legend(loc=0)
    ax_vc.title.set_text('Coverage')
    #
    lambda_v = prm["lambda_v"]
    u_v = 2 * lambda_v * np.ones(t.values.shape[0])
    u_t = 0.0 * t
    ax_uv.plot(t, n_cdmx_vm * u_v, label='$\lambda_v$')
    ax_uv.plot(t, 2 * n_cdmx_vm * u_v, label='$u_{v_{max}} + \lambda_v$')
    ax_uv.legend(loc=0)
    #
    ax_oc_v.plot(df_oc['time'], n_cdmx_vm * df_oc['control'], label='$u_{oc_v}$')
    ax_oc_v.plot(df_oc['time'],
                 n_cdmx_vm * df_oc['control'] +
                 n_cdmx_vm * lambda_v * np.ones(df_oc['time'].values.shape[0]),
                 label='$u_{oc_v}$')
    ax_oc_v.legend(loc=0)
    #
    plt.tight_layout()
    plt.suptitle("R0: " + str(r_zero[0]) + "R_v:" + str(str(r_zero[1])))
    plt.savefig(fig_file_name)
    plt.show()
    return


def reproductive_number(**kwargs):
    beta_s = kwargs['beta_s']
    beta_a = kwargs['beta_a']
    delta_e = kwargs['delta_e']
    alpha_s = kwargs['alpha_s']
    alpha_a = kwargs['alpha_a']
    lambda_v = kwargs['lambda_v']
    epsilon = kwargs['epsilon']
    delta_v = kwargs['delta_v']
    mu = kwargs['mu']
    p = kwargs['p']

    f_1 = (p * delta_e * beta_s) / ((alpha_s + mu) * (mu + delta_e))
    f_2 = ((1.0 - p) * delta_e * beta_a) / ((mu + alpha_a) * (mu + delta_e))
    r_00 = f_1 + f_2
    fac = (mu + delta_v + (1.0 - epsilon) * lambda_v) / \
            (mu + delta_v + lambda_v)
    r_v0 = fac * r_00
    return [r_00, r_v0]


def rhs_vaccination_treatment(t, y, beta_s, beta_a, epsilon,
                              delta_e, delta_v, delta_r,
                              p,
                              alpha_a, alpha_t, alpha_s,
                              mu, theta, mu_t,
                              lambda_v, lambda_t):
    """"
        Version designed in May 22-2020, 
        according to the normalized version and
        lockdown process of the above Kermack Model.
        @param beta_a:
            rate of asymptomatic contact
        @param beta_s:
            rate of asymptomatic contact
        @param epsilon:
            vaccination failure rate
        @param delta_e:
            the inverse denote  latency time
        @param delta_v:
            the inverse denote immunity time
        @param delta_r:
            the inverse denote immunity time
        @param p:
        @param alpha_s:
        @param alpha_a:
        @param alpha_t:
        @type mu:
        @param theta:
        @param mu_t: 
        @param lambda_v:
        @param lambda_t:
        @param y:
        @type t: object
    """
    s, e, i_s, i_a, r, d, v, treat, x_vaccine_counter = y
    #
    a = 1e1
    u_v = 1 / (1 + a * t)
    u_t = 1 / (1 + a * t)
    n_bar = s + e + i_s + i_a + r + v + treat
    force_infection = (beta_s * i_s + beta_a * i_a) / n_bar
    rhs_s = mu * n_bar - force_infection * s - (mu + (lambda_v * u_v)) * s \
            + delta_v * v + delta_r * r
    rhs_e = force_infection * (epsilon * v + s) - (mu + delta_e) * e
    rhs_i_s = p * delta_e * e - (mu + theta + alpha_s + (lambda_t * u_t)) * i_s
    rhs_i_a = (1 - p) * delta_e * e - (mu + alpha_a) * i_a
    rhs_r = alpha_s * i_s + alpha_a * i_a + alpha_t * treat - (mu + delta_r) * r
    rhs_d = theta * i_s + mu_t * treat
    rhs_v = (lambda_v * u_v) * s - epsilon * force_infection * v - \
            (mu + delta_v) * v
    rhs_treat = (lambda_t + u_t) * i_s - (alpha_t + mu_t + mu) * treat
    rhs_x_vaccine_counter = (lambda_v + u_v) * (s + e + i_a + r)
    rhs = np.array([rhs_s, rhs_e,
                    rhs_i_s, rhs_i_a,
                    rhs_r, rhs_d,
                    rhs_v, rhs_treat,
                    rhs_x_vaccine_counter])
    return rhs


def vaccination_control(t, lambda_v):
    u_v_ = 0 * 1.0 * lambda_v
    return u_v_


def rhs_vaccination(t, y, beta_s, beta_a, epsilon,
                    delta_e, delta_v, delta_r,
                    p,
                    alpha_a, alpha_s,
                    mu, theta,
                    lambda_v, control=True):
    """"
        Version designed in May 22-2020,
        according to the normalized version and
        lockdown process of the above Kermack Model.
        @param beta_a:
            rate of asymptomatic contact
        @param beta_s:
            rate of asymptomatic contact
        @param epsilon:
            vaccination failure rate
        @param delta_e:
            the inverse denote  latency time
        @param delta_v:
            the inverse denote immunity time
        @param delta_r:
            the inverse denote immunity time
        @param p:
        @param alpha_s:
        @param alpha_a:
        @type mu:
        @param theta:
        @param lambda_v:
        @param control:
        @param y:
        @type t: object

    Parameters
    ----------
    control
    """
    s, e, i_s, i_a, r, d, v, x_vaccine_counter = y
    #
    if not control:
        lambda_v = 0
    u_v = vaccination_control(t, lambda_v)
    n_bar = s + e + i_s + i_a + r + v
    force_infection = (beta_s * i_s + beta_a * i_a) / n_bar
    rhs_s = mu * n_bar - force_infection * s - (mu + (lambda_v + u_v)) * s \
            + delta_v * v + delta_r * r
    rhs_e = force_infection * ((1.0 - epsilon) * v + s) - (mu + delta_e) * e
    rhs_i_s = p * delta_e * e - (mu + alpha_s) * i_s
    rhs_i_a = (1 - p) * delta_e * e - (mu + alpha_a) * i_a
    rhs_r = theta * alpha_s * i_s + alpha_a * i_a - (mu + delta_r) * r
    rhs_d = (1 - theta) * alpha_s * i_s
    rhs_v = (lambda_v + u_v) * s - (1.0 - epsilon) * force_infection * v \
            - (mu + delta_v) * v
    rhs_x_vaccine_counter = (lambda_v + u_v) * (s + e + i_a + r)
    rhs = np.array([rhs_s, rhs_e,
                    rhs_i_s, rhs_i_a,
                    rhs_r, rhs_d,
                    control * rhs_v,
                    rhs_x_vaccine_counter])
    return rhs
# TODO: Modify plot to stress the area between vaccinated and not vaccinated
# models
