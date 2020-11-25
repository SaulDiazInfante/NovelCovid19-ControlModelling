import numpy as np
from covid19_vaccination import *
from scipy.integrate import solve_ivp
"""
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import Slider, Div
from bokeh.io import curdoc
from bokeh.layouts import column, row
"""
import os
import sys
parameters_dir = '/home/saul/sauld@cimat.mx/UNISON/Articles/NovelCovid-19/' +\
       'NovelCovid19-ControlModelling/' + \
       'NovelCovid19-ControlModellingGitHub/' +\
       'UNISON-UADY-VACCINATION-PRJ/python_sources/' +\
       'vaccination_model_parameters/'
entry = 'vaccination_parameters.json'
full_path = os.path.join(parameters_dir, entry)
prm = load_parameters(full_path)
T = prm["T"]
n_whole = prm["n_pop"]
l_zero = prm["L_0"]
s_zero = prm["S_0"]
e_zero = prm["E_0"]
i_s_zero = prm["I_S_0"]
i_a_zero = prm["I_A_0"]
h_zero = prm["H_0"]
r_zero = prm["R_0"]
d_zero = prm["D_0"]
v_zero = prm["V_0"]
#
#
# ------------------------------------------------------------------------------
t_eval = np.linspace(0, T, num=10000, endpoint=True)
# base dynamics without vaccination
z_0 = np.array([l_zero,
                  s_zero,
                  e_zero,
                  i_s_zero,
                  i_a_zero,
                  h_zero,
                  r_zero,
                  d_zero
                  ])
#
kwargs = {"beta_s": prm["beta_s"],
          "beta_a": prm["beta_a"],
          "kappa": prm["kappa"],
          "delta_l": prm["delta_l"],
          "delta_h": prm["delta_h"],
          "delta_r": prm["delta_r"],
          "p": prm["p"],
          "gamma_a": prm["gamma_a"],
          "gamma_s": prm["gamma_s"],
          "gamma_h": prm["gamma_h"],
          "mu": prm["mu"],
          "mu_i_s": prm["mu_i_s"],
          "mu_h": prm["mu_h"],
          "theta": prm["theta"],
          "var_epsilon": prm["var_epsilon"]
          }
# r_00 = reproductive_number(**kwargs)
df_names = ['time', 'l', 's', 'e',
            'i_s', 'i_a', 'h', 'r',
            'd']
#
base_args = (prm["beta_s"],
        prm["beta_a"],
        prm["kappa"],
        prm["delta_l"],
        prm["delta_h"],
        prm["delta_r"],
        prm["p"],
        prm["gamma_a"],
        prm["gamma_s"],
        prm["gamma_h"],
        prm["mu"],
        prm["mu_i_s"],
        prm["mu_h"],
        prm["theta"],
        prm["var_epsilon"]
    )

sol_not_v = solve_ivp(rhs_base_dynamics,
                      [0.0, T],
                      z_0,
                      dense_output=False,
                      method='LSODA',
                      t_eval=t_eval,
                      events=None,
                      vectorized=False,
                      args=base_args
                      )
df_names = ['time',
            'l', 's', 'e',
            'i_s', 'i_a', 'h',
            'r', 'd']
y = np.c_[sol_not_v.t, sol_not_v.y.T]
data_solution_save(y, df_names, file_name_prefix='base_dynamics')
base_dynamics_plot("base_dynamics.pkl")

