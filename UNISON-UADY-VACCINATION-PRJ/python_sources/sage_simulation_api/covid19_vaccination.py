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


def base_dynamics_plot(data_file_name='base_dynamics.pkl',
                       fig_file_name='base_model_solution.pdf'):
    df_not = pd.read_pickle("base_dynamics.pkl")
    print(data_file_name)
    plt.figure(0)
    # plt.ion()
    ax_l = plt.subplot2grid((3, 3), (0, 0))
    ax_s = plt.subplot2grid((3, 3), (0, 1))
    ax_e = plt.subplot2grid((3, 3), (0, 2))
    #
    ax_i_s = plt.subplot2grid((3, 3), (1, 0))
    ax_i_a = plt.subplot2grid((3, 3), (1, 1))
    ax_h = plt.subplot2grid((3, 3), (1, 2))
    #
    ax_r = plt.subplot2grid((3, 3), (2, 0))
    ax_d = plt.subplot2grid((3, 3), (2, 1))
    ax_cl = plt.subplot2grid((3, 3), (2, 2))
    #
    #
    # prm = load_parameters("vaccination_parameters.json")
    n_cdmx_vm = 26446435.0
   #
    n_whole = 1.0  # without rescaling population
    n_whole = n_cdmx_vm
    t = df_not['time']
    l = df_not['l']
    s = df_not['s']
    e = df_not['e']
    i_s = df_not['i_s']
    i_a = df_not['i_a']
    h = df_not['h']
    r = df_not['r']
    d = df_not['d']
    cl = l + s + e + i_s + i_a + h + r + d
    ax_l.plot(t, n_whole * df_not['l'], 'k--')
    ax_l.legend(loc=0)
    ax_l.title.set_text('Lockdown')
    #
    ax_s.plot(t, n_whole * df_not['s'], 'k--')
    ax_s.legend(loc=0)
    ax_s.title.set_text('Susceptible')
    #
    ax_e.plot(t, n_whole * df_not['e'], 'k--')
    ax_e.legend(loc=0)
    ax_e.title.set_text('Exposed')
    #
    ax_i_s.plot(t, n_whole * df_not['i_s'], 'k--')
    ax_i_s.legend(loc=0)
    ax_i_s.title.set_text('Symptomatic')
    #
    ax_i_a.plot(t, n_whole * df_not['i_a'], 'k--')
    ax_i_a.legend(loc=0)
    ax_i_a.title.set_text('Asymptomatic')
    #
    ax_h.plot(t, n_whole * df_not['h'], 'k--')
    ax_h.legend(loc=0)
    ax_h.title.set_text('Hospitalized')
    #
    ax_r.plot(t, n_whole * df_not['r'], 'k--')
    ax_r.legend(loc=0)
    ax_r.title.set_text('Recovered')
    #
    ax_d.plot(t, n_whole * df_not['d'], 'k--')
    ax_d.legend(loc=0)
    ax_d.title.set_text('Deaths')
    #
    ax_cl.plot(t, cl)
    ax_cl.legend(loc=0)
    ax_cl.title.set_text('CL')
    plt.tight_layout()
    plt.savefig(fig_file_name)
    plt.show()


def vaccination_dynamics_plot(data_file_name='constant_vac_dynamics.pkl',
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
    df_not = pd.read_pickle("base_dynamics.pkl")
    # data_file_name = 'constant_vac_dynamics.pkl'
    # print(data_file_name)
    df_constant = pd.read_pickle("constant_vac_dynamics.pkl")
    # df_oc = pd.read_csv(data_file_name, header=0)
    plt.figure(0)
    # plt.ion()
    ax_l = plt.subplot2grid((4, 3), (0, 0))
    ax_s = plt.subplot2grid((4, 3), (0, 1))
    ax_e = plt.subplot2grid((4, 3), (0, 2))
    #
    ax_i_s = plt.subplot2grid((4, 3), (1, 0))
    ax_i_a = plt.subplot2grid((4, 3), (1, 1))
    ax_h = plt.subplot2grid((4, 3), (1, 2))
   #
    ax_r = plt.subplot2grid((4, 3), (2, 0))
    ax_d = plt.subplot2grid((4, 3), (2, 1))
    ax_v = plt.subplot2grid((4, 3), (2, 2))
   #
    ax_vac = plt.subplot2grid((4, 3), (3, 0))
    ax_cl = plt.subplot2grid((4, 3), (3, 1))
    #
    prm = load_parameters("vaccination_parameters.json")
    n_cdmx_vm = prm['n_pop']
    # n_cul = 905263
    n_whole = 1.0  # without rescaling population
    t = df_constant['time']
    l = df_constant['l']
    s = df_constant['s']
    e = df_constant['e']
    i_s = df_constant['i_s']
    i_a = df_constant['i_a']
    h = df_constant['h']
    r = df_constant['r']
    d = df_constant['d']
    v = df_constant['v']
    x_vaccine = n_whole * (df_constant['x_vac'])
    cl = l + s + e + i_s + i_a + h + r + d + v
    #
    ax_l.plot(t, n_whole * df_not['l'], 'k--')
    ax_l.plot(t, n_whole * l, alpha=0.7)
    ax_l.legend(loc=0)
    ax_l.title.set_text('Lockdown')
    #
    ax_s.plot(t, n_whole * df_not['s'], 'k--')
    ax_s.plot(t, n_whole * s, alpha=0.7)
    ax_s.legend(loc=0)
    ax_s.title.set_text('Susceptible')
    #
    ax_e.plot(t, n_whole * df_not['e'], 'k--')
    ax_e.plot(t, n_whole * e, alpha=0.7)
    ax_e.legend(loc=0)
    ax_e.title.set_text('Exposed')
    #
    ax_i_s.plot(t, n_whole * df_not['i_s'], 'k--')
    ax_i_s.plot(t, n_whole * i_s, label="CC: i_s", alpha=0.7)
    ax_i_s.legend(loc=0)
    ax_i_s.title.set_text('Symptomatic')
    #
    ax_i_a.plot(t, n_whole * df_not['i_a'], 'k--')
    ax_i_a.plot(t, n_whole * i_a, alpha=0.7)
    ax_i_a.legend(loc=0)
    ax_i_a.title.set_text('Asymptomatic')
    #
    ax_h.plot(t, n_whole * df_not['h'], 'k--')
    ax_h.plot(t, n_whole * h, alpha=0.7)
    ax_h.legend(loc=0)
    ax_h.title.set_text('Hospitalized')
    #
    ax_r.plot(t, n_whole * df_not['r'], 'k--')
    ax_r.plot(t, n_whole * r, alpha=0.7)
    ax_r.legend(loc=0)
    ax_r.title.set_text('Recovered')
    #
    ax_d.plot(t, n_whole * df_not['d'], 'k--')
    ax_d.plot(t, n_whole * d, label="CC: d", alpha=0.7)
    ax_d.title.set_text('Deaths')
    #
    ax_v.plot(t, n_whole * v, alpha=0.7)
    ax_v.legend(loc=0)
    ax_v.title.set_text('Vaccinated')
    #
    ax_vac.plot(t, n_whole * x_vaccine)
    ax_vac.legend(loc=0)
    ax_vac.title.set_text('Coverage')
    #
    # https://www.forbes.com.mx/
    # noticias-capacidad-hospitalaria-zona-metropolitana/
    ax_cl.plot(t, cl, label="cl")
    ax_cl.legend(loc=0)
    #
    #
    lambda_v = prm["lambda_v"]
    plt.tight_layout()
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


def rhs_base_dynamics(t, y,
                      beta_s,
                      beta_a,
                      kappa,
                      delta_l,
                      delta_h,
                      delta_r,
                      p,
                      gamma_a,
                      gamma_s,
                      gamma_h,
                      mu,
                      mu_i_s,
                      mu_h,
                      theta,
                      var_epsilon
                      ):
    #
    x_l, x_s, x_e, x_i_s, x_i_a, x_h, x_r, x_d = y
    n_star = x_l + x_s + x_e + x_i_s + x_i_a + x_h + x_r
    infection_force = (beta_s * x_i_s + beta_a * x_i_a) / n_star
    f_l = theta * mu * n_star - \
          (var_epsilon * infection_force + delta_l + mu) * x_l
    f_s = (1.0 - theta) * mu * n_star + delta_l * x_l + delta_r * x_r - \
          (infection_force + mu) * x_s
    f_e = infection_force * (var_epsilon * x_l + x_s) - \
            (kappa + mu) * x_e
    f_i_s = p * kappa * x_e - (gamma_s + delta_h + mu) * x_i_s
    f_i_a = (1.0 - p) * kappa * x_e - (gamma_a + mu) * x_i_a
    f_h = delta_h * x_i_s - (gamma_h + mu_h + mu) * x_h
    f_r = gamma_s * x_i_s + gamma_a * x_i_a + gamma_h * x_h - \
        (delta_r + mu) * x_r
    f_d = mu_h * x_h
    f = np.array([f_l, f_s, f_e, f_i_s, f_i_a, f_h, f_r, f_d])
    return f


def rhs_lockdown_vaccination(t, y, beta_s, beta_a, epsilon, var_epsilon,
                    delta_l, kappa, delta_h, delta_v, delta_r,
                    p, gamma_a, gamma_s, gamma_h,
                    mu, mu_i_s, mu_h, theta, lambda_v):
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
        @param kappa:
            the inverse denote  latency time
        @param delta_l:
        @param delta_v: the inverse denote immunity time
        @param delta_h: hospitalization rate
        @param delta_r: the inverse denotes natural immunity
        @param p:
        @param gamma_s:
        @param gamma_a:
        @param gamma_h:
        @type mu:
        @param mu_h:
        @param theta:
        @param lambda_v:
        @param y:
        @type t: object
    """
    l, s, e, i_s, i_a, h, r, d, v, x_vaccine_counter = y
    #
    n_bar = l + s + e + i_s + i_a + h + r + v
    delta_v = delta_v
    force_infection = (beta_s * i_s + beta_a * i_a) / n_bar
    rhs_l = theta * mu * n_bar - \
            (var_epsilon * force_infection + delta_l + mu) * l
    rhs_s = (1.0 - theta) * mu * n_bar + delta_l * l + delta_v * v + \
            delta_r * r - (force_infection + lambda_v + mu) * s
    rhs_e = (var_epsilon * l + (1-epsilon) * v + s) * force_infection - \
            (kappa + mu) * e
    rhs_i_s = p * kappa * e - (delta_h + gamma_s + mu_i_s + mu) * i_s
    rhs_i_a = (1 - p) * kappa * e - (gamma_a + mu) * i_a
    rhs_h = delta_h * i_s - (gamma_h + mu_h + mu) * h
    rhs_r = gamma_s * i_s + gamma_a * i_a + gamma_h * h - \
            (delta_r + mu) * r
    rhs_d = mu_i_s * i_s + mu_h * h
    rhs_v = lambda_v * s - ((1-epsilon) * force_infection + delta_v + mu) * v
    rhs_x_vaccine_counter = lambda_v * (s + e + i_a + r)
    rhs = np.array([rhs_l, rhs_s, rhs_e,
                    rhs_i_s, rhs_i_a,
                    rhs_h, rhs_r, rhs_d,
                    rhs_v,  rhs_x_vaccine_counter])
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
