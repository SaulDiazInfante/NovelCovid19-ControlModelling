import re
import numpy as np
import json
import pandas as pd
from datetime import datetime


def read_bocop_parameters_keys(bocop_file_name='solution_deltav_0.sol',
                               parameters_prefix_file_name=\
                                       'vaccination_parameters_'
                               ):
    # TODO: Recover parameters key from bocop_run_parameters.json
    star_pointer_block = 0
    end_pointer_block = 0
    state_keys = []
    control_keys = []
    parameter_keys = []
    boundarycond_keys = []
    constraint_keys = []
    constant_keys = []
    path = bocop_file_name
    time_now = datetime.now()
    dt_string = time_now.strftime("%b-%d-%Y_%H_%M")
    path_json = parameters_prefix_file_name + dt_string + '.json'
    with open(path) as f:
        all_lines = f.readlines()
    i = 0
    pattern_str_names_label = \
        re.compile(r'# # Names :')
    pattern_str_end_block = \
        re.compile(r'problem.bounds')
    pattern_state_token = re.compile('state\.\d')
    pattern_control_token = re.compile('control\.\d')
    pattern_parameter_token = re.compile('parameter\.\d')
    pattern_boundary_token = re.compile('boundarycond\.\d')
    pattern_constraint_token = re.compile('constraint\.\d')
    pattern_constant_token = re.compile('constant\.\d')
    patter_split_word = re.compile('\S+')

    for line in all_lines:
        str_headers_names = pattern_str_names_label.search(line)
        str_footer_names = pattern_str_end_block.search(line)
        if str_headers_names:
            star_pointer_block = i + 1
            print(str_headers_names.string, i)
        if str_footer_names:
            print(str_footer_names.string, i)
            end_pointer_block = i - 6
        i = i + 1
    for j in np.arange(np.int(star_pointer_block), np.int(end_pointer_block)):
        x_str = patter_split_word.findall(all_lines[j])
        x_identifier = x_str[1]
        x_str_state = pattern_state_token.search(x_identifier)
        x_str_control = pattern_control_token.search(x_identifier)
        x_str_parameter = pattern_parameter_token.search(x_identifier)
        x_str_boundarycond = pattern_boundary_token.search(x_identifier)
        x_str_constraint = pattern_constraint_token.search(x_identifier)
        x_str_constant = pattern_constant_token.search(x_identifier)
        if x_str_state:
            str_state = x_str[3]
            #            print('\t', str_state, '\t\tat line:', j)
            state_keys.append(str_state)
        if x_str_control:
            str_control = x_str[3]
            #            print('\t', str_control, '\t\tat line:', j)
            control_keys.append(str_control)
        if x_str_parameter:
            str_parameter = x_str[3]
            #            print('\t', str_parameter, '\t\tat line:', j)
            parameter_keys.append(str_parameter)
        if x_str_boundarycond:
            str_boundarycond = x_str[3]
            #            print('\t', str_boundarycond, '\t\tat line:', j)
            boundarycond_keys.append(str_boundarycond)
        if x_str_constraint:
            str_constraint = x_str[3]
            #            print('\t', str_constraint, '\t\tat line:', j)
            constraint_keys.append(str_constraint)
        if x_str_constant:
            str_constant = x_str[3]
            #            print('\t', str_constant, '\t\tat line:', j)
            constant_keys.append(str_constant)

    return (state_keys, control_keys, parameter_keys, boundarycond_keys,
            constraint_keys, constant_keys)


# state_keys_, control_keys_, parameter_keys_, \
# boundarycond_keys_, constraint_keys_, constant_keys_ = \
#    read_bocop_parameters_keys(bocop_file_name='solution_deltav_0.sol',
#                               parameters_prefix_file_name='vaccination_parameters_')
