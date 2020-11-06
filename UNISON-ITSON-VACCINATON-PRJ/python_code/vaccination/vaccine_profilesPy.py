import matplotlib.pyplot as plt
from sage_scenarios_simulation import *
import os
import sys
from sage_incidence_plots import *
from vaccine_profiles import *


def get_solutions_file_list(root_dir):
    """
        For the given path,
        get the List of all files in the directory tree
    """
    # create a list of file and sub directories
    # names in the given directory
    list_of_files = os.listdir(root_dir)
    all_files = list()
    # Iterate over all the entries
    for entry in list_of_files:
        # Create full path
        full_path = os.path.join(root_dir, entry)
        # If entry is a directory then get the list of files in this directory
        if os.path.isdir(full_path):
            all_files = all_files + get_solutions_file_list(full_path)
        else:
            all_files.append(full_path)
    return all_files


def main():
    bocop_data_path = \
        '/home/saul/sauld@cimat.mx/UNISON/Articles/NovelCovid-19/' + \
        '/NovelCovid19-ControlModelling/NovelCovid19-ControlModellingGitHub' + \
        '/UNISON-ITSON-VACCINATON-PRJ/python_code/vaccination/bocop_data' + \
        '/Solutions_WHO_Paper_102520/delta_V_365/' + \
        '/coverage_20_6M/delta_R_180/epsilon-5.000000e-01.sol'

    csv_data_path = \
        '/home/saul/sauld@cimat.mx/UNISON/Articles/NovelCovid-19/' + \
        '/NovelCovid19-ControlModelling/NovelCovid19-ControlModellingGitHub' + \
        '/UNISON-ITSON-VACCINATON-PRJ/python_code/vaccination/' + \
        'vaccination_pkl_solutions/csv'
    path_solutions_list = get_solutions_file_list(csv_data_path)
    path_list = list()
    path_list_0 = list()
    for (dir_path, dir_names, filenames) in os.walk(csv_data_path):
        path_list += [os.path.join(dir_path, file) for file in filenames]
    # Print the files
    with open('solution_list.txt', 'w') as f:
        for item in path_list:
            f.write("%s\n" % item)
    # for elem in path_list:
    #    print(elem)
    bocop_parameters_json_file = './vaccination_model_parameters/' + \
                                 'vaccination_parameters.json'
    path_list_0.append(bocop_data_path)
    for bocop_main_scene_data in path_list_0:
        print(bocop_main_scene_data)
        cvp = Covid19VaccineProfiles(
            bocop_solution_file=bocop_main_scene_data,
            bocop_parameters_json_file= bocop_parameters_json_file)
        cvp.__init__(
            bocop_solution_file=bocop_main_scene_data,
            bocop_parameters_json_file= bocop_parameters_json_file)
        cvp.plot_vaccine_efficiency_profile(
                                        coverage=0.5,
                                        vaccine_immunity=730.0,
                                        natural_immunity=0.005555556,
                                        time_horizon=180.0
        )
        cvp.plot_vaccine_efficiency_profile(
            coverage=0.5,
            vaccine_immunity=365.0,
            natural_immunity=0.005555556,
            time_horizon=180.0
        )

if __name__ == "__main__":
    main()
