import decimal
import json
import os
import re
from datetime import *
import bokeh.palettes as bokeh_palettes
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
# from IPython.display import Image
import plotly.io as pio
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
from scipy import integrate
from matplotlib import rcParams, cycler
from pathlib import Path

from sage_scenarios_simulation import *


class Covid19VaccineProfiles(CovidNumericalModel):

    def __init__(self,
                 bocop_solution_file,
                 bocop_parameters_json_file):
        super().__init__(bocop_solution_file,
                         bocop_parameters_json_file)
        self.csv_data_list = list()

    def get_solutions_file_list(self, root_dir):
        """
            For the given path,
            get the List of all files in the directory tree
        """
        # create a list of file and sub directories
        # names in the given directory
        list_of_files = os.listdir(root_dir)
        all_files_list = list()
        # Iterate over all the entries
        for entry in list_of_files:
            # Create full path
            full_path = os.path.join(root_dir, entry)
            # If entry is a directory then get the list of files
            # in this directory
            if os.path.isdir(full_path):
                all_files_list = all_files_list + \
                                 self.get_solutions_file_list(full_path)
            else:
                all_files_list.append(full_path)
        return all_files_list

    def load_data(self):
        bocop_data_path = \
            '/home/saul/sauld@cimat.mx/UNISON/Articles/NovelCovid-19/' + \
            '/NovelCovid19-ControlModelling/' + \
            'NovelCovid19-ControlModellingGitHub' + \
            '/UNISON-ITSON-VACCINATON-PRJ/python_code/' + \
            'vaccination/' + 'vaccination_pkl_solutions/csv/'
        path_solutions_list = self.get_solutions_file_list(bocop_data_path)
        path_list = list()
        for (dir_path, dir_names, filenames) in os.walk(bocop_data_path):
            path_list += [os.path.join(dir_path, file) for file in filenames]
        self.csv_data_list = path_list.sort()

    def plot_vaccine_efficiency_profile(self,
                                        coverage=0.5,
                                        vaccine_immunity=730.0,
                                        natural_immunity=0.005555556,
                                        time_horizon=180.0):
        self.load_data()
        data_path = 'vaccination_pkl_solutions/csv/'
        cov_tag = str(coverage)
        time_horizon_tag = str(time_horizon)
        delta_v_tag = str(vaccine_immunity)
        delta_r_tag = str(natural_immunity)
        epsilon_tag = '*'
        query_tag = 'bocop_run' + \
                    '_coverage_' + cov_tag + \
                    '_tHor_' + time_horizon_tag + \
                    '_delta_r_' + delta_r_tag + \
                    '_delta_v_' + delta_v_tag + \
                    '_epsilon_' + epsilon_tag + \
                    '_solution' + '.csv'
        print("Data for contrast efficiency:")
        print("=============================")
        data_solution_list = list()
        for data_file in Path(data_path).glob(query_tag):
            data_solution_list.append(data_file)
        #
        data_solution_list.sort()
        for data_file in data_solution_list:
            print('\t' + str(data_file))
        plt.style.use('seaborn-paper')
        cmap = plt.cm.inferno
        # TODO: fix color palette
        rcParams['axes.prop_cycle'] = \
        cycler(color=cmap(np.linspace(0, 1, len(data_solution_list))))
        fig, ax = plt.subplots()
        for data_csv_file in data_solution_list:
            df = pd.read_csv(data_csv_file)
            epsilon_tag = re.findall('\d+\.\d+', str(data_csv_file))[-1]
            ax.plot(df['time'], 100000 * df['d'],
                    label= r'$\epsilon = $' + epsilon_tag)
            ax.set_xlabel('Time (days)')
            ax.set_ylabel('Symptomatic cases per 100,000 inhabitants')
        title_str = 'Coverage: ' + cov_tag + \
                    ', T:' + time_horizon_tag + \
                    r', $\delta_R$: ' + delta_r_tag + \
                    r', $1 / \delta_v$: ' + delta_v_tag
        plt.suptitle(title_str, fontsize=12)
        plt.legend(loc=0)
        #fig_name = 'fig_symptomatic_cases' + \
        #            '_coverage_' + cov_tag + \
        #            '_tHor_' + time_horizon_tag + \
        #            '_delta_r_' + delta_r_tag + \
        #            '_delta_v_' + delta_v_tag + \
        #            '_epsilon_' + 'x'\
        #            '_solution' + '.png'
        fig_name = 'fig_deaths_cases' + \
                    '_coverage_' + cov_tag + \
                    '_tHor_' + time_horizon_tag + \
                    '_delta_r_' + delta_r_tag + \
                    '_delta_v_' + delta_v_tag + \
                    '_epsilon_' + 'x'\
                    '_solution' + '.png'
        path = '/home/saul/sauld@cimat.mx/UNISON/Articles/NovelCovid-19/' + \
            '/NovelCovid19-ControlModelling/' + \
            'NovelCovid19-ControlModellingGitHub' + \
            '/UNISON-ITSON-VACCINATON-PRJ/python_code/' + \
            'vaccination/images/png/vaccine_profile/' + fig_name
        fig_size = plt.rcParams["figure.figsize"]
        fig_size[0] = 7.48031
        fig_size[1] = 7.48031 / 1.6180
        plt.rcParams["figure.figsize"] = fig_size
        fig.savefig(path, dpi=1000)
        plt.show()

