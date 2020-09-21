import re
import numpy as np
import json
import pandas as pd
from datetime import datetime
from get_key_string_model_parameters import *


def read_bocop_parameters_values(bocop_file_name='solution_deltav_0.sol',
                                 parameters_prefix_file_name=
                                         'vaccination_parameters'):
    state_keys, control_keys, parameter_keys, \
        boundarycond_keys, constraint_keys, constant_keys = \
        read_bocop_parameters_keys(bocop_file_name='solution_deltav_0.sol',
                                   parameters_prefix_file_name=
                                   'vaccination_parameters')
    path = bocop_file_name
    time_now = datetime.now()
    dt_string = time_now.strftime("%b-%d-%Y_%H_%M")
    path_json = parameters_prefix_file_name + '_' + dt_string + '.json'
    constants_values = []
    parameters_values = []
    boundarycond_values = []
    boundary_cond_dim = boundarycond_keys.__len__()
    state_dim = state_keys.__len__()
    constant_dim = constant_keys.__len__()
    with open(path) as f:
        all_lines = f.readlines()
    i = 0
    pattern_str_constants_label = \
        re.compile(r'problem.constants')
    pattern_str_bounds_label = \
        re.compile(r'problem.bounds')
    pattern_str_optimal_time = \
        re.compile(r'# Parameters')
    pattern_value = re.compile(r'\d+\.\d+|\d+')
    pattern_sci_float = re.compile('[-+]?[\d]+\.?[\d]*[Ee](?:[-+]?[\d]+)?')
    value = 0.0
    for line in all_lines:
        str_header_constants = pattern_str_constants_label.search(line)
        str_header_bounds = pattern_str_bounds_label.search(line)
        str_header_optimal_time = pattern_str_optimal_time.search(line)
        if str_header_bounds:
            print(i, str_header_bounds.string)
            pointer = i + 11
            for j in np.arange(pointer, pointer + boundary_cond_dim):
                x_val_sci = pattern_sci_float.search(all_lines[j])
                x_val = pattern_value.search(all_lines[j])
                if x_val_sci:
                    print(j + 1, x_val_sci.string)
                    value = x_val_sci.group()
                elif x_val:
                    print(j + 1, x_val.string)
                    value = x_val.group()
                boundarycond_values.append(np.float(value))
        #
        if str_header_constants:
            print(str_header_constants.string, i)
            pointer = i + 8
            for j in np.arange(pointer, pointer + constant_dim):
                x_val_sci = pattern_sci_float.search(all_lines[j])
                x_val = pattern_value.search(all_lines[j])
                if x_val_sci:
                    print(j + 1, x_val_sci.string)
                    value = x_val_sci.group()
                elif x_val:
                    print(j + 1, x_val.string)
                    value = x_val.group()
                constants_values.append(np.float(value))
        #
        if str_header_optimal_time:
            print(str_header_optimal_time.string, i)
            pointer = i + 1
            x_val_sci = pattern_sci_float.search(all_lines[pointer])
            x_val = pattern_value.search(all_lines[pointer])
            if x_val_sci:
                print(pointer, x_val_sci.string)
                value = x_val_sci.group()
                print(value)
            elif x_val:
                print(pointer, x_val.string)
                value = x_val.group()
            constants_values.append(np.float(value))
            constant_keys.append('T')
        i = i + 1
    data_constants = dict(zip(constant_keys, constants_values))
    data_boundarycond = dict(zip(boundarycond_keys, boundarycond_values))
    data = data_constants
    data.update(data_boundarycond)
    data_json = json.dumps(data)
    print(json.dumps(data, indent=4))
    with open(path_json, 'w') as json_file:
        json.dump(data, json_file, indent=4)
    with open(parameters_prefix_file_name + '.json', 'w') as json_file:
        json.dump(data, json_file, indent=4)
    df_record_parameters = pd.DataFrame(data, index=[1])
    return df_record_parameters


bocop_main_scene_data = 'bocop_data/'\
                  + 'deltaR_05_Coverage_50_TH_10M_14M/delta_V_05_years/' \
                  + 'epsilon-7.000000e-01.sol'

parameters = \
    read_bocop_parameters_values(
        bocop_file_name=bocop_main_scene_data,
        parameters_prefix_file_name='vaccination_parameters')
