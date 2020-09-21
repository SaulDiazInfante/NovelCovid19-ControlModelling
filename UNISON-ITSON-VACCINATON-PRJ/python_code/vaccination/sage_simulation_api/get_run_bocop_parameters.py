import re
import numpy as np
import pandas as pd
import json
from datetime import datetime


def get_run_bocop_parameters(bocop_file_name='solution_deltav_0.sol',
                             prefix_file_name='bocop_run_parameters'
                             ):
    def check_token(x_str_, x_val_, i_):
        if x_str_ and x_val_:
            string_prop_name_ = x_str.group()
            str_prop_val_ = np.int(x_val[0])
            names.append(string_prop_name_)
            values.append(str_prop_val_)
            print(i_, string_prop_name_, str_prop_val_)

    def check_token_str(x_str_, x_val_, i_):
        if x_str_ and x_val_:
            string_prop_name_ = x_str.group()
            str_prop_val_ = x_val_[0]
            names.append(string_prop_name_)
            values.append(str_prop_val_)
            print(i_, string_prop_name_, str_prop_val_)

    path = bocop_file_name
    names = []
    values = []
    with open(path) as f:
        all_lines = f.readlines()
    i = 1
    for line in all_lines:
        # time regular expressions
        pattern_str = re.compile(r'time\.([a-z]*)')
        # pattern_str_none = re.compile(r'none')
        pattern_str_free = re.compile(r'final')
        pattern_int_value = re.compile(r'\d+')

        x_str = pattern_str.search(line)
        x_val = pattern_int_value.findall(line)
        x_str_final = pattern_str_free.findall(line)
        # x_str_none = pattern_str_free.findall(line)
        if x_str and x_val:
            string_prop_name = 'time.' + x_str.group(1)
            str_prop_val = np.int(x_val[0])
            names.append(string_prop_name)
            values.append(str_prop_val)
            print(i, string_prop_name, str_prop_val)

        elif x_str and x_str_final:
            string_prop_name = 'time.' + 'free'
            str_prop_val = x_str_final[0]
            names.append(string_prop_name)
            values.append(str_prop_val)
            print(i, string_prop_name, str_prop_val)
        # Dimension re
        pattern_str_state_dimension = re.compile(r'state.dimension')
        pattern_str_control = re.compile(r'control.dimension')
        pattern_str_algebraic = re.compile(r'algebraic.dimension')
        pattern_str_parameter = re.compile(r'parameter.dimension')
        pattern_str_constant = re.compile(r'constant.dimension')
        pattern_str_boundarycond = re.compile(r'boundarycond.dimension')
        pattern_str_constraint = re.compile(r'constraint.dimension')
        pattern_int_value = re.compile(r'\d+')
        # state
        x_str = pattern_str_state_dimension.search(line)
        x_val = pattern_int_value.findall(line)
        check_token(x_str, x_val, i)
        # control
        x_str = pattern_str_control.search(line)
        x_val = pattern_int_value.findall(line)
        check_token(x_str, x_val, i)
        # algebraic
        x_str = pattern_str_algebraic.search(line)
        x_val = pattern_int_value.findall(line)
        check_token(x_str, x_val, i)
        # parameter
        x_str = pattern_str_parameter.search(line)
        x_val = pattern_int_value.findall(line)
        check_token(x_str, x_val, i)
        # constants
        x_str = pattern_str_constant.search(line)
        x_val = pattern_int_value.findall(line)
        check_token(x_str, x_val, i)
        # boundarycond
        x_str = pattern_str_boundarycond.search(line)
        x_val = pattern_int_value.findall(line)
        check_token(x_str, x_val, i)
        x_str = pattern_str_constraint.search(line)
        x_val = pattern_int_value.findall(line)
        check_token(x_str, x_val, i)
        #
        # Discretization re:
        pattern_str_steps = re.compile(r'discretization.steps')
        pattern_str_method = re.compile(r'discretization.method')
        pattern_value_str_method = \
            re.compile(r'discretization.method string\s+(\S+)')
        #
        x_str = pattern_str_steps.search(line)
        x_val = pattern_int_value.findall(line)
        check_token(x_str, x_val, i)
        #
        x_str = pattern_str_method.search(line)
        x_val = pattern_value_str_method.findall(line)
        check_token_str(x_str, x_val, i)
        #
        # Optimization:
        #
        pattern_str_opt_type = re.compile(r'optimization.type')
        pattern_value_str_opt_type = \
            re.compile(r'optimization.type string\s+(\S+)')
        # pattern_str_opt_batch_type = re.compile(r'batch.type')
        pattern_str_opt_batch_index = re.compile(r'batch.index')
        pattern_str_opt_batch_nrange = re.compile(r'batch.nrange')
        pattern_str_opt_batch_lowerbound = re.compile(r'batch.lowerbound')
        pattern_str_opt_batch_upperbound = re.compile(r'batch.upperbound')
        pattern_str_opt_batch_directory = re.compile(r'batch.directory')
        pattern_value_str_opt_directory = \
            re.compile(r'batch.directory string\s+(\S+)')
        # type
        x_str = pattern_str_opt_type.search(line)
        x_val = pattern_value_str_opt_type.findall(line)
        check_token_str(x_str, x_val, i)
        # index
        x_str = pattern_str_opt_batch_index.search(line)
        x_val = pattern_int_value.findall(line)
        check_token(x_str, x_val, i)
        # nrange
        x_str = pattern_str_opt_batch_nrange.search(line)
        x_val = pattern_int_value.findall(line)
        check_token(x_str, x_val, i)
        # lowerbound
        x_str = pattern_str_opt_batch_lowerbound.search(line)
        x_val = pattern_int_value.findall(line)
        check_token(x_str, x_val, i)
        # upperbound
        x_str = pattern_str_opt_batch_upperbound.search(line)
        x_val = pattern_int_value.findall(line)
        check_token(x_str, x_val, i)
        # directory
        x_str = pattern_str_opt_batch_directory.search(line)
        x_val = pattern_value_str_opt_directory.findall(line)
        check_token_str(x_str, x_val, i)
        # solution file name
        pattern_str_solution_file = re.compile(r'solution.file')
        pattern_value_str_solution_file = \
            re.compile(r'solution.file string\s+(\S+)')

        x_str = pattern_str_solution_file.search(line)
        x_val = pattern_value_str_solution_file.findall(line)
        check_token_str(x_str, x_val, i)
        i = i + 1
    time_now = datetime.now()
    dt_string = time_now.strftime("%b-%d-%Y_%H_%M")
    path_json = prefix_file_name + '_' + dt_string + '.json'
    #
    bocop_data = dict(zip(names, values))
    bocop_data_df = pd.DataFrame(bocop_data, index=[1])
    bocop_data_df.to_json(path_json)
    print(json.dumps(bocop_data, indent=4))
    with open(path_json, 'w') as json_file:
        json.dump(bocop_data, json_file, indent=4)
    with open(prefix_file_name + '.json', 'w') as json_file:
        json.dump(bocop_data, json_file, indent=4)
    return bocop_data_df


df_bcp_rp = get_run_bocop_parameters(bocop_file_name='solution_deltav_0.sol',
                                     prefix_file_name='bocop_run_parameters'
                                     )
