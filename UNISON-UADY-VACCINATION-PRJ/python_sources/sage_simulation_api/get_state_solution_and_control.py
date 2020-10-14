import re
import numpy as np
import pandas as pd
from datetime import datetime
import json
from get_key_string_model_parameters import *


def get_solution_data(bocop_parameters_file,
                      col_names,
                      bocop_file_name='bocop_data/'
                                      + 'deltaR_05_Coverage_50_TH_10M_14M/'
                                      + 'delta_V_05_years/'
                                      + 'epsilon-7.000000e-01.sol',
                      out_put_file_prefix_name='bocop_'
                      ):
    state_keys, control_keys, parameter_keys, \
    boundarycond_keys, constraint_keys, constant_keys = \
        read_bocop_parameters_keys(bocop_file_name= bocop_file_name,
                                   parameters_prefix_file_name=
                                   'vaccination_parameters')
    bocop_prm = pd.read_json(bocop_parameters_file, typ='series')
    offset = bocop_prm["discretization.steps"]
    with open(bocop_file_name) as f:
        all_lines = f.readlines()
    i = 0
    t_s = []
    data = []
    control = []
    final_time = 1.0
    pattern_str_time_solution_header = re.compile(r'SOLUTION')
    pattern_str_state_solution_i_header = re.compile(r'State\s\d')
    pattern_str_state_control_i_header = re.compile(r'Control\s\d')
    pattern_str_optimized_time = re.compile(r'Parameters')
    pattern_value = re.compile(r'\d+\.\d+|\d+')
    pattern_sci_float = re.compile('[-+]?[\d]+\.?[\d]*[Ee](?:[-+]?[\d]+)?')
    value = 0.0
    for line in all_lines:
        x_str_time_header = pattern_str_time_solution_header.search(line)
        x_str_state_i_header = pattern_str_state_solution_i_header.search(line)
        x_str_control_i_header = pattern_str_state_control_i_header.search(line)
        x_str_optimized_time_header = pattern_str_optimized_time.search(line)
        if x_str_time_header:
            pointer = i + 13
            n_time_steps = np.int(1 / np.float(all_lines[pointer + 1]))
            print(pointer, all_lines[pointer], 'n_time_steps:', n_time_steps)
            time_str_block = np.array(all_lines[pointer: pointer +
                                                         n_time_steps + 1])
            for j in np.arange(time_str_block.shape[0]):
                x_val_sci = pattern_sci_float.search(time_str_block[j])
                x_val = pattern_value.search(time_str_block[j])
                if x_val_sci:
                    print(j + 1, x_val_sci.string)
                    value = x_val_sci.group()
                elif x_val:
                    print(j + 1, x_val.string)
                    value = x_val.group()
                t_s.append(np.float(value))
        #
        if x_str_state_i_header:
            pointer = i + 1
            data_block = np.array(all_lines[pointer: pointer + offset])
            data_i = []
            for j in np.arange(data_block.shape[0]):
                x_val_sci = pattern_sci_float.search(data_block[j])
                x_val = pattern_value.search(data_block[j])
                if x_val_sci:
                    print(j + 1, x_val_sci.string)
                    value = x_val_sci.group()
                elif x_val:
                    print(j + 1, x_val.string)
                    value = x_val.group()
                data_i.append(np.float(value))
            data_block = np.array(data_i)
            data.append(data_block)
        if x_str_optimized_time_header:
            pointer = i
            final_time = all_lines[pointer + 1]
            x_val_sci = pattern_sci_float.search(final_time)
            x_val = pattern_value.search(final_time)
            if x_val_sci:
                print(j + 1, x_val_sci.string)
                value = x_val_sci.group()
            elif x_val:
                print(j + 1, x_val.string)
                value = x_val.group()
            final_time = np.float(value)
        if x_str_control_i_header:
            pointer = i + 1
            control_block = np.array(all_lines[pointer: pointer + offset])
            control_i = []
            for j in np.arange(control_block.shape[0]):
                x_val_sci = pattern_sci_float.search(control_block[j])
                x_val = pattern_value.search(control_block[j])
                if x_val_sci:
                    print(j + 1, x_val_sci.string)
                    value = x_val_sci.group()
                elif x_val:
                    print(j + 1, x_val.string)
                    value = x_val.group()
                control_i.append(np.float(value))
            control_block = np.array(control_i)
            control.append(control_block)
        i = i + 1
    control = np.array(control).T
    state_control_data = np.array(data).T
    t_s = final_time * np.array(t_s)
    t_s = t_s.reshape([t_s.shape[0], 1])
    # plt.plot(t_s[0:-1:2], state_control_data[:, 1])
    if t_s.shape[0] > state_control_data.shape[0]:
        offset = np.int(t_s.shape[0] / state_control_data.shape[0])
        t_s = t_s[0:-1: offset]
    data = np.concatenate([t_s,
                           state_control_data,
                           control],
                          axis=1)
    time_now = datetime.now()
    dt_string = time_now.strftime("%b-%d-%Y_%H_%M")
    path_pkl = out_put_file_prefix_name + '_' + dt_string + '.pkl'
    df_bocop_data = pd.DataFrame(data=data, columns=col_names)
    file_name = out_put_file_prefix_name + 'solution.pkl'
    df_bocop_data.to_pickle(file_name)
    df_bocop_data.to_pickle(path_pkl)
    return df_bocop_data


col_names_ = ['time', 's', 'e',
              'i_s', 'i_a', 'r', 'd', 'v',
              'x_vac', 'y_inc(t)', 'x_zero',
              'u_V']
base_path = '../' + 'bocop_run_parameters/'
bocop_parameters_file_ = base_path + 'bocop_run_parameters.json'
bocop_main_scene_data = '../bocop_data/'\
                  + 'delta_V_0.5/' \
                  + 'solution_epsilon_085_cV_10.sol'
df = get_solution_data(bocop_parameters_file=bocop_parameters_file_,
                       col_names=col_names_,
                       bocop_file_name=bocop_main_scene_data)
