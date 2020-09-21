import numpy as np
import re
import json
import pandas as pd
import os
import plotly.express as px
import plotly.graph_objects as go
import decimal
from numerics_covid_ode_models import NumericsCovid19
from datetime import datetime
from plotly.subplots import make_subplots
from scipy import integrate
import bokeh.palettes as bokeh_palettes
# from IPython.display import Image
import plotly.io as pio


class CovidModels(NumericsCovid19):

    def __init__(self,
                 uncontrolled_solution_path='./vaccination_pkl_solutions/' +
                                            'no_vaccination_data_solution.pkl',
                 optimal_controlled_solution='./vaccination_pkl_solutions/' +
                                             'bocop_solution.pkl',
                 constant_controlled_solution='./vaccination_pkl_solutions/' +
                                              'constant_vaccination_' +
                                              'data_solution.pkl',
                 bocop_parameters_json_file='./bocop_run_parameters/' +
                                            'bocop_run_parameters' +
                                            '.json',
                 vaccination_parameters_json_file='./vaccination_model' +
                                                  '_parameters/'
                                                  + 'vaccination_parameters' +
                                                  '.json',
                 bocop_solution_file='./bocop_data/' + 'solution_deltav_0.sol'
                 ):
        super().__init__()
        with open(vaccination_parameters_json_file) as json_file:
            vaccination_prm_ = json.load(json_file)
        self.parameters = vaccination_prm_
        self.r00 = 0.0
        self.r_vac = 0.0
        self.df_uncontrolled_sol = pd.read_pickle(uncontrolled_solution_path)
        self.df_controlled_sol = pd.read_pickle(optimal_controlled_solution)
        self.bocop_solution_file = bocop_solution_file
        self.bocop_parameters_json_file = bocop_parameters_json_file
        self.constant_controlled_solution_file = constant_controlled_solution
        self.uncontrolled_solution_path = uncontrolled_solution_path
        self.optima_controlled_solution_path = optimal_controlled_solution

        # self.df_constant_vaccination_solution = pd.

    def get_run_bocop_parameters(self,
                                 prefix_file_name='bocop_run_parameters'
                                 ):
        bocop_file_name = self.bocop_solution_file

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

    def load_parameters(self,
                        file_name='vaccination_parameters.json'):
        """ Load the parameters given a JSON file.
        Parameters
        ----------
        file_name: filename with parameters values in JSON format.

        Returns
        -------
        prm: dictionary
        """
        with open(file_name) as json_file:
            prm = json.load(json_file)
        self.parameters = prm
        return prm

    def read_bocop_parameters_keys(self,
                                   parameters_prefix_file_name=
                                   'vaccination_parameters_'):
        # TODO: Recover parameters key from bocop_run_parameters.json
        bocop_file_name = self.bocop_solution_file
        star_pointer_block = 0
        end_pointer_block = 0
        state_keys = []
        control_keys = []
        parameter_keys = []
        boundarycond_keys = []
        constraint_keys = []
        constant_keys = []
        path = bocop_file_name
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
        for j in np.arange(np.int(star_pointer_block),
                           np.int(end_pointer_block)):
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
        self.state_keys = state_keys
        self.control_keys = control_keys
        self.parameter_keys = parameter_keys
        self.boundarycond_keys = boundarycond_keys
        self.constraint_keys = constraint_keys
        self.control_keys = constant_keys
        return (state_keys, control_keys, parameter_keys, boundarycond_keys,
                constraint_keys, constant_keys)

    def get_solution_data(self,
                          out_put_file_prefix_name='bocop_'
                          ):
        bocop_file_name = self.bocop_solution_file
        state_keys, control_keys, parameter_keys, boundarycond_keys, \
            constraint_keys, constant_keys = \
            self.read_bocop_parameters_keys(parameters_prefix_file_name=
                                                'vaccination_parameters')

        bocop_parameters_file = self.bocop_parameters_json_file
        bocop_prm = pd.read_json(bocop_parameters_file, typ='series')
        # self.time_state_solution_names
        offset = bocop_prm["discretization.steps"]
        with open(bocop_file_name) as f:
            all_lines = f.readlines()
        i = 0
        j = 0
        t_s = []
        data = []
        control = []
        final_time = 1.0
        pattern_str_time_solution_header = re.compile(r'SOLUTION')
        pattern_str_state_solution_i_header = re.compile(r'State\s\d')
        pattern_str_state_control_i_header = re.compile(r'Control\s\d')
        pattern_str_optimized_time = re.compile(r'Parameters')
        pattern_value = re.compile(r'\d+\.\d+|\d+')
        pattern_sci_float = re.compile(
            '[-+]?[\d]+\.?[\d]*[Ee](?:[-+]?[\d]+)?')
        value = 0.0
        for line in all_lines:
            x_str_time_header = pattern_str_time_solution_header.search(
                line)
            x_str_state_i_header = pattern_str_state_solution_i_header.search(
                line)
            x_str_control_i_header = pattern_str_state_control_i_header.search(
                line)
            x_str_optimized_time_header = pattern_str_optimized_time.search(
                line)
            if x_str_time_header:
                pointer = i + 13
                n_time_steps = np.int(1 / np.float(all_lines[pointer + 1]))
                print(pointer, all_lines[pointer], 'n_time_steps:',
                      n_time_steps)
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
                control_block = np.array(
                    all_lines[pointer: pointer + offset])
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
        data = np.c_[t_s, state_control_data, control]
        time_now = datetime.now()
        dt_string = time_now.strftime("%b-%d-%Y_%H_%M")
        path_pkl = out_put_file_prefix_name + '_' + dt_string + '.pkl'
        col_names = ['time', 's', 'e',
                     'i_s', 'i_a', 'r', 'd', 'v',
                     'x_vac', 'x_zero',
                     'u_V']
        df_bocop_data = pd.DataFrame(data=data, columns=col_names)
        file_name = out_put_file_prefix_name + 'solution.pkl'
        df_bocop_data.to_pickle(file_name)
        df_bocop_data.to_pickle(path_pkl)
        return df_bocop_data

    def reproductive_number(self):
        beta_s = self.parameters['beta_s']
        beta_a = self.parameters['beta_a']
        delta_e = self.parameters['delta_e']
        alpha_s = self.parameters['alpha_s']
        alpha_a = self.parameters['alpha_a']
        lambda_v = self.parameters['lambda_v']
        epsilon = self.parameters['epsilon']
        delta_v = self.parameters['delta_v']
        mu = self.parameters['mu']
        p = self.parameters['p']

        f_1 = (p * delta_e * beta_s) / ((alpha_s + mu) * (mu + delta_e))
        f_2 = ((1.0 - p) * delta_e * beta_a) / ((mu + alpha_a) * (mu + delta_e))
        r_00 = f_1 + f_2

        fac = (mu + delta_v + (1.0 - epsilon) * lambda_v) / \
              (mu + delta_v + lambda_v)
        u_v_0 = self.df_controlled_sol['u_V'][0]
        opt_fac = (mu + delta_v + (1.0 - epsilon) * (lambda_v + u_v_0)) / \
                  (mu + delta_v + (lambda_v + u_v_0))

        r_v_0 = fac * r_00
        opt_r_v_0 = opt_fac * r_00
        return r_00, r_v_0, opt_r_v_0

    def sage_scenarios_plot(self):
        def hex_to_rgb(hex_color: str) -> tuple:
            hex_color = hex_color.lstrip("#")
            if len(hex_color) == 3:
                hex_color = hex_color * 2
            return int(hex_color[0:2], 16), \
                   int(hex_color[2:4], 16), int(hex_color[4:6], 16)

        def constant_vaccination_cost(df_solution):
            """

            """
            a_S = self.parameters['a_S']
            a_D = self.parameters['a_D']
            c_V = self.parameters['c_V']
            u_V = self.parameters['lambda_v']
            t = df_solution['time']
            i_s = df_solution['i_s']
            d = df_solution['d']
            f = a_S * i_s + a_D * d + 0.5 * c_V * u_V ** 2
            cost_t_ = integrate.cumtrapz(f, t, initial=0)
            return cost_t_

        def constant_vaccination_coverage(df_solution):
            """

            """
            lambda_v = self.parameters['lambda_v']
            t = df_solution['time']
            s = df_solution['s']
            e = df_solution['e']
            i_a = df_solution['i_a']
            r = df_solution['r']
            f = lambda_v * (s + e + i_a + r)
            coverage_t_ = integrate.cumtrapz(f, t, initial=0)
            return coverage_t_

        prm = self.parameters
        lambda_v_base = prm["lambda_v"]
        bocop_sol_path = self.optima_controlled_solution_path
        constant_vaccination_path = self.constant_controlled_solution_file
        no_vaccination_sol_path = self.uncontrolled_solution_path
        df_oc = pd.read_pickle(bocop_sol_path)
        df_not_vaccination = pd.read_pickle(no_vaccination_sol_path)
        df_constant_vaccination = pd.read_pickle(constant_vaccination_path)
        border_color_pallet = px.colors.sequential.ice
        fill_color_pallet = bokeh_palettes.all_palettes['Category20'][20]

        fig = make_subplots(
            rows=3, cols=3,
            specs=[
                [{}, {}, {"rowspan": 2}],
                [{}, {}, None],
                [{"colspan": 2}, None, {}]
            ],
            subplot_titles=("Symptomatic",
                            "Hospitalised",
                            "Death",
                            "Asymptomatic",
                            "Coverage",
                            "Vaccination Policy",
                            "Cost"
                            )
        )
        n_cdmx = self.parameters['n_pop']
        fig.add_trace(
            go.Scatter(
                x=df_oc['time'],
                y=n_cdmx * df_oc['i_s'],
                fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[14]), 0.5)}",
                legendgroup='Mitigation',
                line=dict(color=border_color_pallet[1], width=.7),
                name='(OP)',
                showlegend=False
            ),
            row=1, col=1
        )
        fig.add_trace(
            go.Scatter(
                x=df_oc['time'],
                y=0.25 ** 2 * n_cdmx * df_oc['i_s'],
                line=dict(color=border_color_pallet[1], width=.7),
                legendgroup='Optimized resources',
                showlegend=False,
            ),
            row=1, col=2)
        fig.add_trace(
            go.Scatter(
                x=df_oc['time'],
                y=n_cdmx * df_oc['d'],
                # fillcolor=fill_color_pallet[5],
                # fill='tonexty',
                line=dict(color=border_color_pallet[1], width=.7),
                showlegend=False
            ),
            row=1, col=3
        )
        #
        fig.add_trace(
            go.Scatter(
                x=df_oc['time'],
                y=n_cdmx * df_oc['i_a'],
                fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[14]), 0.5)}",
                # fill='tonexty',
                line=dict(color=border_color_pallet[1], width=.7),
                showlegend=False
            ),
            row=2, col=1)
        fig.add_trace(
            go.Scatter(
                x=df_oc['time'],
                y=df_oc['x_vac'],
                line=dict(color=fill_color_pallet[19], width=.7),
                fill='tozeroy',
                fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[19]), 0.5)}",
                legendgroup='Coverage',
                name="(OP)",
                showlegend=True
            ),
            row=2, col=2
        )
        #
        lambda_base_v_t = lambda_v_base * np.ones(df_oc["u_V"].shape[0])
        fig.add_trace(
            go.Scatter(
                x=df_oc['time'],
                y=n_cdmx * (df_oc['u_V'] + lambda_base_v_t),
                fill='tozeroy',
                fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[19]), 0.5)}",
                line=dict(color=fill_color_pallet[0], width=0.7),
                name="(OP)",
                legendgroup='Policy',
                showlegend=False
            ),
            row=3, col=1)
        fig.add_trace(
            go.Scatter(
                x=df_oc['time'],
                y=df_oc['x_zero'],
                line=dict(color=fill_color_pallet[16], width=0.7),
                fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[17]), 0.5)}",
                fill='tozeroy',
                legendgroup='Cost',
                name='(OP)',
                showlegend=True
            ),
            row=3, col=3
        )

        trace_vaccination_base = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * lambda_v_base * np.ones(df_oc["u_V"].shape[0]),
            line=dict(color=fill_color_pallet[0], width=0.7, dash='dot'),
            showlegend=False
        )
        trace_vaccination_max = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * (
                    df_oc["u_V"].max()
                    + lambda_v_base * np.ones(df_oc["u_V"].shape[0])
            ),
            line=dict(color=fill_color_pallet[0], width=0.7, dash='dot'),
            showlegend=False
        )
        ########################################################################
        # constant_vaccination
        ########################################################################
        trace_constant_vac_i_s = go.Scatter(
            x=df_constant_vaccination['time'],
            y=n_cdmx * df_constant_vaccination['i_s'],
            fill='tonexty',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[2]), 0.7)}",
            line=dict(color=fill_color_pallet[2], width=0.7),
            legendgroup='Mitigation',
            name='(CP)'
        )
        trace_constant_vac_i_a = go.Scatter(
            x=df_constant_vaccination['time'],
            y=n_cdmx * df_constant_vaccination['i_a'],
            fill='tonexty',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[2]), 0.7)}",
            line=dict(color=fill_color_pallet[2], width=0.7),
            showlegend=False
        )
        trace_constant_vac_hospitalization = go.Scatter(
            x=df_constant_vaccination['time'],
            y=0.25 ** 2 * n_cdmx * df_constant_vaccination['i_s'],
            fill='tonexty',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[6]), 0.7)}",
            line=dict(color=fill_color_pallet[6], width=0.7),
            legendgroup='Optimized resources',
            name='(CP)',
            showlegend=True
        )
        trace_constant_vac_d = \
            go.Scatter(
                x=df_constant_vaccination['time'],
                y=n_cdmx * df_constant_vaccination['d'],
                fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[4]), 0.7)}",
                fill='tonexty',
                line=dict(color=fill_color_pallet[4], width=0.7),
                legendgroup='Saved lives',
                name='Saved lives (CP)'
            )

        coverage_t = constant_vaccination_coverage(df_constant_vaccination)
        trace_constant_vac_coverage = go.Scatter(
            x=df_constant_vaccination['time'],
            y=coverage_t,
            fill='tozeroy',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[18]), 0.7)}",
            line=dict(color=fill_color_pallet[18], width=0.7),
            legendgroup='Coverage',
            name='(CP)',
            showlegend=True
        )
        cost_t = constant_vaccination_cost(df_constant_vaccination)
        trace_constant_vac_cost = go.Scatter(
            x=df_constant_vaccination['time'],
            y=cost_t,
            fill='tonexty',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[16]), 0.7)}",
            line=dict(color=fill_color_pallet[16], width=0.7),
            legendgroup='Cost',
            name='(CP)',
            showlegend=True
        )
        ########################################################################
        # Solution without vaccination
        ########################################################################
        trace_no_vac_i_s = go.Scatter(
            x=df_not_vaccination['time'],
            y=n_cdmx * df_not_vaccination['i_s'],
            fill='tonexty',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[2]), 0.5)}",
            line=dict(color=border_color_pallet[1], width=0.7, dash='dot'),
            legendgroup='Mitigation',
            name='(OP)'
        )
        fig.add_trace(
            go.Scatter(
                x=df_not_vaccination['time'],
                y=0.25 ** 2 * n_cdmx * df_not_vaccination['i_s'],
                line=dict(color=border_color_pallet[1], width=0.7, dash='dot'),
                fill='tonexty',
                fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[7]), 0.5)}",
                showlegend=True,
                legendgroup='Optimized resources',
                name='OP'
            ),
            row=1, col=2)

        trace_no_vac_i_a = go.Scatter(
            x=df_not_vaccination['time'],
            y=n_cdmx * df_not_vaccination['i_a'],
            fill='tonexty',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[2]), 0.5)}",
            line=dict(color=border_color_pallet[1], width=0.7, dash='dot'),
            legendgroup='Mitigation',
            showlegend=False
        )
        trace_no_vac_d = \
            go.Scatter(
                x=df_not_vaccination['time'],
                y=n_cdmx * df_not_vaccination['d'],
                fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[5]), 0.5)}",
                fill='tonexty',
                line=dict(color=border_color_pallet[1], width=.70, dash='dot'),
                legendgroup='Saved lives',
                name='Saved lives (OP)'
            )

        fig.append_trace(trace_no_vac_i_s, 1, 1)
        fig.append_trace(trace_constant_vac_i_s, 1, 1)
        fig.append_trace(trace_constant_vac_hospitalization, 1, 2)
        fig.append_trace(trace_no_vac_i_a, 2, 1)
        fig.append_trace(trace_constant_vac_i_a, 2, 1)
        fig.append_trace(trace_constant_vac_coverage, 2, 2)
        fig.append_trace(trace_no_vac_d, 1, 3)
        fig.append_trace(trace_constant_vac_d, 1, 3)
        fig.append_trace(trace_vaccination_base, 3, 1)
        fig.append_trace(trace_vaccination_max, 3, 1)
        fig.append_trace(trace_constant_vac_cost, 3, 3)
        #
        r_0, r_v_0, r_opt_v_0 = self.reproductive_number()
        c_r_0 = decimal.Decimal(r_0)
        c_r_v_0 = decimal.Decimal(r_v_0)
        c_r_opt_v_0 = decimal.Decimal(r_opt_v_0)
        fig.update_layout(showlegend=True,
                          title_text="R0: " + str(round(c_r_0, 6))
                                     + '\t\tRv: ' + str(round(c_r_v_0, 6))
                                     + '\t\tOC-Rv: '
                                     + str(round(c_r_opt_v_0, 6)),
                          # paper_bgcolor='rgba(0,0,0,0)',
                          # plot_bgcolor='rgba(0,0,0,0)',
                          # yaxis={
                          #    "tickcolor": "rgba(186,186,186,0.75)",
                          #    "tickwidth": 15,
                          #    "gridcolor": "rgba(186,186,186,0.75)",
                          #    "gridwidth": 1,
                          #    "zerolinecolor": "gray",
                          #    "zerolinewidth": 2,},
                          # xaxis={
                          #  "tickcolor": "gray",
                          #  "tickwidth": 10,
                          #  "gridcolor": "gray",
                          #  "gridwidth": 2,}
                          )
        if not os.path.exists("images"):
            os.mkdir("images")
        # fig.write_image("images/fig1.pdf")
        golden_width = 718  # width in px
        golden_ratio = 1.618
        pio.kaleido.scope.default_format = "eps"
        pio.kaleido.scope.default_width = golden_width
        pio.kaleido.scope.default_height = golden_width / golden_ratio
        pio.kaleido.scope.default_scale = 0.5
        fig.to_image(format="pdf", engine="kaleido")
        fig.write_image("images/fig1.pdf")
        # TODO:  Edit legend respect to groups
        fig.show()
