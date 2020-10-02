import matplotlib.pyplot as plt
from sage_scenarios_simulation import *
from sage_incidence_plots import *
bocop_main_scene_data = './bocop_data/' + 'deltaR_05_Coverage_50_TH_10M_14M/' \
                        + 'delta_V_05_years/' \
                        + 'epsilon-7.000000e-01.sol'
# + 'delta_V_lifelong/'
# + 'deltaR_05_Coverage_50_TH_10M_14M'
# +'/delta_V_05_years/'
# + 'epsilon-7.000000e-01.sol'


scm = CovidNumericalModel(bocop_solution_file=bocop_main_scene_data)
scm.get_run_bocop_parameters()
scm.read_bocop_parameters_values()
scm.get_lsoda_solution(constant_control=False)
scm.get_lsoda_solution(constant_control=True)
scm.get_solution_data()
fig01 = scm.sage_plot_prevalence()
fig02 = scm.sage_plot_optimal_signal()
scm_incidence = CovidNumericalModelIncidence()
scm_incidence.sage_plot_incidence()

