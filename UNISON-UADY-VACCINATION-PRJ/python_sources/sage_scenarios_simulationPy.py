import matplotlib.pyplot as plt
from sage_scenarios_simulation import *
from sage_incidence_plots import *
path_list = []
f = open("solution_list.txt", "r")
#
while True:
    path = f.readline()
    if not path:
        break
    path_list.append(path.strip())
    # print(path.strip())
f.close()
path_list_0 = [path_list[1]]
#
#
bocop_parameters_json_file = './vaccination_model_parameters/' + \
    'vaccination_parameters.json'
for bocop_main_scene_data in path_list_0:
    print(bocop_main_scene_data)
    scm = CovidNumericalModel(
        bocop_solution_file=bocop_main_scene_data,
        bocop_parameters_json_file= bocop_parameters_json_file)
    scm.__init__(
        bocop_solution_file=bocop_main_scene_data,
        bocop_parameters_json_file= bocop_parameters_json_file)
    scm.get_run_bocop_parameters()
    scm.read_bocop_parameters_values()
    scm.get_lsoda_solution(constant_control=False)
    scm.get_lsoda_solution(constant_control=True)
    scm.get_solution_data()
    scm.load_pkl_data()
    scm.simulation_run_tagger()
    scm.sage_plot_prevalence()
    scm.sage_plot_optimal_signal()
    scm_incidence = CovidNumericalModelIncidence(
        bocop_solution_file=bocop_main_scene_data,
        bocop_parameters_json_file=bocop_parameters_json_file)
    scm_incidence.__init__(
        bocop_solution_file=bocop_main_scene_data,
        bocop_parameters_json_file=bocop_parameters_json_file)
    scm_incidence.simulation_run_tagger()
    scm_incidence.sage_plot_incidence()
