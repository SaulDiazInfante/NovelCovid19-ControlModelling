import matplotlib.pyplot as plt
from sage_scenarios_simulation import *
from sage_incidence_plots import *
path_list = []
f = open("solution_list.txt", "r")

while True:
    path = f.readline()
    if not path:
        break
    path_list.append(path.strip())
    print(path.strip())
f.close()
bocop_main_scene_data = path_list.pop()
scm = CovidNumericalModel(bocop_solution_file=bocop_main_scene_data)
scm.get_run_bocop_parameters()
scm.read_bocop_parameters_values()
scm.get_lsoda_solution(constant_control=False)
scm.get_lsoda_solution(constant_control=True)
scm.get_solution_data()
scm.sage_plot_prevalence()
scm.sage_plot_optimal_signal()
scm_incidence = CovidNumericalModelIncidence()
scm_incidence.sage_plot_incidence()