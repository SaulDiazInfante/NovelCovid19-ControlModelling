from covid19_vaccination import *
from scipy.integrate import solve_ivp
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import Slider, Div
from bokeh.io import curdoc
from bokeh.layouts import column, row

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
kwargs = {"s_0": prm["s_0"],
          "beta_s": prm["beta_s"],
          "beta_a": prm["beta_a"],
          "delta_e": prm["delta_e"],
          "p": prm["p"],
          "alpha_a": prm["alpha_a"],
          "alpha_s": prm["alpha_s"],
          "mu": prm["mu"],
          "theta": prm["theta"],
          "lambda_v": prm["lambda_v"],
          "delta_v": prm["delta_v"],
          "epsilon": prm["epsilon"]
          }
r_00 = reproductive_number(**kwargs)

df_names = ['time', 's', 'e',
            'i_s', 'i_a', 'r',
            'd', 'v', 'treat',
            'x_vac_counter']
# Only vaccination
# TODO: Check parameters
z_0 = np.array([s_zero,
                e_zero,
                i_s_zero,
                i_a_zero,
                r_zero,
                d_zero,
                v_zero,
                x_vaccination_zero
                ])
#
args = (prm["beta_s"], prm["beta_a"], prm["epsilon"],
        prm["delta_e"], prm["delta_v"], prm["delta_r"],
        prm["p"],
        prm["alpha_a"], prm["alpha_s"],
        prm["mu"], prm["theta"],
        prm["lambda_v"], not (bool(prm["control"])))

sol_not_v = solve_ivp(rhs_vaccination,
                      [0.0, T],
                      z_0,
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
args = (prm["beta_s"], prm["beta_a"], prm["epsilon"],
        prm["delta_e"], prm["delta_v"], prm["delta_r"],
        prm["p"],
        prm["alpha_a"], prm["alpha_s"],
        prm["mu"], prm["theta"],
        prm["lambda_v"], bool(prm["control"]))

sol_v = solve_ivp(rhs_vaccination,
                  [0.0, T],
                  z_0,
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

vaccination_population_plot(r_00,
                            data_file_name='vaccination_solution.pkl')
