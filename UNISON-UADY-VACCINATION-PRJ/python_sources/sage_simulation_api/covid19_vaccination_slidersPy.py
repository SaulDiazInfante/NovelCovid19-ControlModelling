import bokeh.plotting as bp
import pandas as pd
from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import (Div, Slider)
from scipy.integrate import solve_ivp

from covid19_vaccination import *


def update_data(attrname, old, new):
    lambda_v_ = lambda_v.value
    epsilon_ = epsilon.value
    print(epsilon_)
    print(lambda_v_)
    args_slider_ = (prm["beta_s"],
                    prm["beta_a"],
                    epsilon_,
                    prm["delta_e"],
                    prm["delta_v"],
                    prm["delta_r"],
                    prm["p"],
                    prm["alpha_a"],
                    prm["alpha_s"],
                    prm["mu"],
                    prm["mu_s"],
                    lambda_v_,
                    bool(prm["control"])
                    )
    kwargs = {"s_0": prm["s_0"],
              "beta_s": prm["beta_s"],
              "beta_a": prm["beta_a"],
              "delta_e": prm["delta_e"],
              "p": prm["p"],
              "alpha_a": prm["alpha_a"],
              "alpha_s": prm["alpha_s"],
              "mu": prm["mu"],
              "mu_s": prm["mu_s"],
              }
    r_00 = reproductive_number(**kwargs)
    str_r_not = "R_zero = {0:.2f}".format(r_00)
    sol_ = solve_ivp(rhs_vaccination,
                     [0.0, T],
                     z_v_0,
                     dense_output=False,
                     method='LSODA',
                     t_eval=t_eval,
                     events=None,
                     vectorized=False,
                     args=args_slider_
                     )
    s, e, i_s, i_a, r, d, v, x_v = sol_.y
    df['s'] = s
    df['e'] = e
    df['i_s'] = i_s
    df['i_a'] = i_a
    df['r'] = r
    df['d'] = d
    df['v'] = v
    df['x_v'] = x_v
    ax_s.line('time', 's',
              line_color='gray',
              alpha=0.5,
              source=df
              )
    ax_e.line('time', 'e',
              line_color='gray',
              source=df
              )
    ax_i_s.line('time', 'i_s',
                line_color='gray',
                alpha=0.5,
                source=df
                )
    ax_i_a.line('time', 'i_a',
                line_color='gray',
                alpha=0.5,
                source=df
                )
    ax_r.line('time', 'r',
              line_color='gray',
              alpha=0.5,
              source=df
              )
    ax_d.line('time', 'd',
              line_color='gray',
              alpha=0.5,
              source=df
              )
    ax_v.line('time', 'v',
              line_color='gray',
              alpha=0.5,
              source=df
              )
    ax_x_v.line('time', 'x_vac_counter',
                line_color='gray',
                alpha=0.5,
                source=df)
    inputs = column(Div(text="<h3>SEIR Model Input Values</h3>" + str_r_not),
                    lambda_v, epsilon)


prm = load_parameters("model_vaccination_parameters.json")
T = prm["T"]
n_whole = prm["n_whole"]
s_zero = prm["s_0"]
e_zero = prm["e_0"]
i_s_zero = prm["i_s_0"]
i_a_zero = prm["i_a_0"]
r_zero = prm["r_0"]
d_zero = prm["d_0"]
v_zero = prm["v_0"]
treat_zero = prm["treat_0"]
x_vaccination_zero = prm["x_v_0"]
# ------------------------------------------------------------------------------
t_eval = np.linspace(0, T, num=2000, endpoint=True)
#
z_v_0 = np.array([s_zero,
                  e_zero,
                  i_s_zero,
                  i_a_zero,
                  r_zero,
                  d_zero,
                  v_zero,
                  x_vaccination_zero
                  ])
#
lambda_v = Slider(title="lambda_v",
                  value=prm["lambda_v"],
                  start=0.0001,
                  end=0.008,
                  step=0.0001)
epsilon = Slider(title="epsilon",
                 value=prm["epsilon"],
                 start=0.1 * prm["epsilon"],
                 end=1.0,
                 step=.01)
#
args_slider = (prm["beta_s"],
               prm["beta_a"],
               epsilon.value,
               prm["delta_e"], prm["delta_v"], prm["delta_r"],
               prm["p"],
               prm["alpha_a"], prm["alpha_s"],
               prm["mu"], prm["mu_s"],
               lambda_v.value,
               bool(prm["control"])
               )
#
#
kwargs = {"s_0": prm["s_0"],
          "beta_s": prm["beta_s"],
          "beta_a": prm["beta_a"],
          "delta_e": prm["delta_e"],
          "p": prm["p"],
          "alpha_a": prm["alpha_a"],
          "alpha_s": prm["alpha_s"],
          "mu": prm["mu"],
          "mu_s": prm["mu_s"],
          }
# TODO: Implement R0 for vaccination
r_00 = reproductive_number(**kwargs)
# vaccination_population_plot(r_00,
#                           data_file_name='vaccination_treatment_solution.pkl')
# Only vaccination
# TODO: Check parameters
#
args = (prm["beta_s"],
        prm["beta_a"],
        prm["epsilon"],
        prm["delta_e"],
        prm["delta_v"],
        prm["delta_r"],
        prm["p"],
        prm["alpha_a"],
        prm["alpha_s"],
        prm["mu"],
        prm["mu_s"],
        prm["lambda_v"],
        not (bool(prm["control"])))
#
sol_not_v = solve_ivp(rhs_vaccination,
                      [0.0, T],
                      z_v_0,
                      dense_output=False,
                      method='LSODA',
                      t_eval=t_eval,
                      events=None,
                      vectorized=False,
                      args=args
                      )
df_names = ['time', 's', 'e',
            'i_s', 'i_a', 'r',
            'd', 'v',
            'x_vac_counter']
y = np.c_[sol_not_v.t, sol_not_v.y.T]
data_solution_save(y, df_names,
                   file_name_prefix='not_vaccination_solution'
                   )
data_file_name = "not_vaccination_solution.pkl"
df_not = pd.read_pickle(data_file_name)
args = (prm["beta_s"],
        prm["beta_a"],
        prm["epsilon"],
        prm["delta_e"],
        prm["delta_v"],
        prm["delta_r"],
        prm["p"],
        prm["alpha_a"],
        prm["alpha_s"],
        prm["mu"],
        prm["mu_s"],
        prm["lambda_v"],
        bool(prm["control"]))
#
sol_v = solve_ivp(rhs_vaccination,
                  [0.0, T],
                  z_v_0,
                  dense_output=False,
                  method='LSODA',
                  t_eval=t_eval,
                  events=None,
                  vectorized=False,
                  args=args
                  )
y = np.c_[sol_v.t, sol_v.y.T]
data_solution_save(y, df_names,
                   file_name_prefix='vaccination_solution'
                   )
data_file_name = "vaccination_solution.pkl"
df = pd.read_pickle(data_file_name)
str_r_not = "R_zero = {0:.2f}".format(r_00)
#
for w in [lambda_v, epsilon]:
    w.on_change('value', update_data)

# Plot the data on three separate curves for S(t), I(t) and R(t)
ax_s = bp.figure(width=300, height=200, title='Susceptible')
ax_e = bp.figure(width=300, height=200, title='Exposed')
ax_i_s = bp.figure(width=300, height=200, title='Symptomatic')
ax_i_a = bp.figure(width=300, height=200, title='Asymptomatic')
ax_r = bp.figure(width=300, height=200, title='Recovered')
ax_d = bp.figure(width=300, height=200, title='Deaths')
ax_v = bp.figure(width=300, height=200, title='Vaccinated')
ax_x_v = bp.figure(width=300, height=200, title='Vaccinated Coverage')

ax_s.line('time', 's',
          line_color='blue',
          alpha=0.5,
          source=df_not)
ax_s.line('time', 's',
          line_color='gray',
          alpha=0.5,
          source=df)

ax_e.line('time', 'e',
          line_color='orange',
          source=df_not)
ax_e.line('time', 'e',
          line_color='gray',
          source=df)

ax_i_s.line('time', 'i_s',
            line_color='red',
            alpha=0.5,
            source=df_not)
ax_i_s.line('time', 'i_s',
            line_color='gray',
            alpha=0.5,
            source=df)

ax_i_a.line('time', 'i_a',
            line_color='darkgreen',
            alpha=0.5,
            source=df_not)
ax_i_a.line('time', 'i_a',
            line_color='gray',
            alpha=0.5,
            source=df)

ax_r.line('time', 'r',
          line_color='blue',
          alpha=0.5,
          source=df_not)
ax_r.line('time', 'r',
          line_color='gray',
          alpha=0.5,
          source=df)

ax_d.line('time', 'd',
          line_color='blue',
          alpha=0.5,
          source=df_not)
ax_d.line('time', 'd',
          line_color='gray',
          alpha=0.5,
          source=df)

ax_v.line('time', 'v',
          line_color='blue',
          alpha=0.5,
          source=df_not)
ax_v.line('time', 'v',
          line_color='gray',
          alpha=0.5,
          source=df)

ax_x_v.line('time', 'v',
            line_color='orange',
            alpha=0.5,
            source=df_not)
ax_x_v.line('time', 'v',
            line_color='gray',
            alpha=0.5,
            source=df)
# bp.show(fig)
# Set up layouts and add to document
inputs = column(Div(text="<h3>SEIR Model Input Values</h3>" + str_r_not),
                lambda_v, epsilon)
curdoc().add_root(row(inputs, ax_s, ax_e, width=200))
curdoc().add_root(row(ax_i_s, ax_i_a, width=200))
curdoc().add_root(row(ax_r, ax_d, width=200))
curdoc().add_root(row(ax_v, ax_x_v, width=200))
curdoc().title = "Simple SIR Model " + str_r_not
