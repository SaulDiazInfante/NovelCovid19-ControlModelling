from lockdown_exit_dynamics import *
from scipy.integrate import solve_ivp

prm = load_parameters("exit_lockdown_parameters_may_29_2020.json")
T = prm["T"]
n_whole = prm["n_whole"]
l_zero = 0.5  # prm["l_0"]  # * n_whole
s_zero = 0.5 - 2.0 / n_whole  # prm["s_0"]
e_zero = prm["e_0"]
i_s_zero = prm["i_s_0"]
i_a_zero = prm["i_a_0"]
h_zero = prm["h_0"]
r_zero = prm["r_0"]
d_zero = prm["d_0"]


# ---------------------------------------------------------------------------------
t_eval = np.linspace(0, T, num=10000, endpoint=True)
#
z_lockdown_0 = np.array([l_zero, s_zero, e_zero, i_s_zero,
                                           i_a_zero, h_zero, r_zero, d_zero])

#
args = (prm["beta_s"], prm["beta_a"], prm["epsilon"],
        prm["delta_u"], prm["delta_e"], prm["delta_s"], prm["p"],
        prm["alpha_a"], prm["alpha_h"], prm["alpha_s"], prm["mu_s"], prm["mu_h"])
#
#
sol = solve_ivp(rhs_lockdown_exit_01,
                [0.0, T],
                z_lockdown_0,
                dense_output=False,
                method='RK45',
                t_eval=t_eval,
                events=None,
                vectorized=False,
                args=args
                )

kwargs = {"l_0": prm["l_0"],
          "s_0": prm["s_0"],
          "beta_s": prm["beta_s"],
          "beta_a": prm["beta_a"],
          "p": prm["p"],
          "alpha_a": prm["alpha_a"],
          "alpha_s": prm["alpha_s"]
          }
r_zero = reproductive_number(**kwargs)
df_names = ['time', 'l', 's', 'e', 'i_s', 'i_a', 'h', 'r', 'd']
y = np.c_[sol.t, sol.y.T]
data_solution_save(y, df_names)
lockdown_population_plot(r_zero, data_file_name='lock_down_data_solution2020-05-30.pkl')
print("r_zero: ", r_zero)

# TODO: Check parameters
