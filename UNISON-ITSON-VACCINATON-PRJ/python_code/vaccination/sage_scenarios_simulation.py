import decimal
import json
import os
import re
import bokeh.palettes as bokeh_palettes
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from scipy import integrate
from datetime import datetime
from numerics_covid_ode_models import NumericsCovid19


def hex_to_rgb(hex_color: str) -> tuple:
    hex_color = hex_color.lstrip("#")
    if len(hex_color) == 3:
        hex_color = hex_color * 2
    return int(hex_color[0:2], 16), \
           int(hex_color[2:4], 16), int(hex_color[4:6], 16)


class CovidNumericalModel(NumericsCovid19):

    def __init__(self,
                 bocop_solution_file,
                 bocop_parameters_json_file,
                 coverage=0.2,
                 ):
        super().__init__(bocop_solution_file)
        self.r00 = 0.0
        self.r_vac = 0.0
        self.daly_constant_cost = 0.0
        self.daly_optimal_cost = 0.0
        self.coverage = coverage
        self.df_uncontrolled_sol = pd.DataFrame()
        self.df_controlled_sol = pd.DataFrame()
        self.constant_controlled_solution_file = 'None'
        self.uncontrolled_solution_path = 'None'
        self.optimal_controlled_solution_path = 'None'
        self.run_tag = 'None'
        self.bocop_parameters_json_file = bocop_parameters_json_file
        self.bocop_solution_file = bocop_solution_file

    def load_pkl_data(self,
                      uncontrolled_solution_path=
                      './vaccination_pkl_solutions/' +
                      'no_vaccination_data_solution.pkl',
                      optimal_controlled_solution=
                      './vaccination_pkl_solutions/' +
                      'bocop_solution.pkl',
                      constant_controlled_solution=
                      './vaccination_pkl_solutions/' +
                      'constant_vaccination_' +
                      'data_solution.pkl',
                      vaccination_parameters_json_file=
                      './vaccination_model' + '_parameters/' +
                      'vaccination_parameters' + '.json'):
        with open(vaccination_parameters_json_file) as json_file:
            vaccination_prm_ = json.load(json_file)
        self.parameters = vaccination_prm_
        self.df_uncontrolled_sol = pd.read_pickle(uncontrolled_solution_path)
        self.df_controlled_sol = pd.read_pickle(optimal_controlled_solution)
        self.constant_controlled_solution_file = constant_controlled_solution
        self.uncontrolled_solution_path = uncontrolled_solution_path
        self.optimal_controlled_solution_path = optimal_controlled_solution

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
                # # # print(str_headers_names.string, i)
            if str_footer_names:
                # # # print(str_footer_names.string, i)
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
                #            # # # print('\t', str_state, '\t\tat line:', j)
                state_keys.append(str_state)
            if x_str_control:
                str_control = x_str[3]
                #            # # # print('\t', str_control, '\t\tat line:', j)
                control_keys.append(str_control)
            if x_str_parameter:
                str_parameter = x_str[3]
                #            # # # print('\t', str_parameter, '\t\tat line:', j)
                parameter_keys.append(str_parameter)
            if x_str_boundarycond:
                str_boundarycond = x_str[3]
                # # # # print('\t', str_boundarycond, '\t\tat line:', j)
                boundarycond_keys.append(str_boundarycond)
            if x_str_constraint:
                str_constraint = x_str[3]
                # # # # print('\t', str_constraint, '\t\tat line:', j)
                constraint_keys.append(str_constraint)
            if x_str_constant:
                str_constant = x_str[3]
                #            # # # print('\t', str_constant, '\t\tat line:', j)
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
                          out_put_file_prefix_name=
                          './vaccination_pkl_solutions' + '/bocop'
                          ):
        """

        Parameters
        ----------
        out_put_file_prefix_name : object
        """
        bocop_file_name = self.bocop_solution_file
        state_keys, control_keys, parameter_keys, boundarycond_keys, \
        constraint_keys, constant_keys = \
            self.read_bocop_parameters_keys(parameters_prefix_file_name=
                                            'vaccination_parameters')

        bocop_parameters_file = self.bocop_parameters_json_file
        bocop_prm = self.bocop_run_parameters
        # pd.read_json(bocop_parameters_file, typ='series')
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
        pattern_sci_float = re.compile('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?')
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
                # # # print(pointer, all_lines[pointer], 'n_time_steps:',
                # n_time_steps)
                time_str_block = np.array(all_lines[pointer: pointer +
                                                             n_time_steps + 1])
                for j in np.arange(time_str_block.shape[0]):
                    x_val_sci = pattern_sci_float.search(time_str_block[j])
                    x_val = pattern_value.search(time_str_block[j])
                    if x_val_sci:
                        value = x_val_sci.group()
                        # # # print(j + 1, value)
                    elif x_val:
                        value = x_val.group()
                        # # print(j + 1, value)
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
                        value = x_val_sci.group()
                        # # print(j + 1, value)

                    elif x_val:
                        value = x_val.group()
                        # # print(j + 1, value)
                    data_i.append(np.float(value))
                data_block = np.array(data_i)
                data.append(data_block)
            if x_str_optimized_time_header:
                pointer = i
                final_time = all_lines[pointer + 1]
                x_val_sci = pattern_sci_float.search(final_time)
                x_val = pattern_value.search(final_time)
                if x_val_sci:
                    value = x_val_sci.group()
                    # # print(j + 1, value)
                elif x_val:
                    value = x_val.group()
                    # # print(j + 1, value)
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
                        value = x_val_sci.group()
                        # # print(j + 1, value)
                    elif x_val:
                        value = x_val.group()
                        # # print(j + 1, value)
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
        path_csv = out_put_file_prefix_name + '_' + dt_string + '.csv'
        col_names = ['time', 's', 'e',
                     'i_s', 'i_a', 'r', 'd', 'v',
                     'x_vac', 'y_inc', 'x_zero',
                     'u_V']
        df_bocop_data = pd.DataFrame(data=data, columns=col_names)
        file_name = out_put_file_prefix_name + '_solution.pkl'
        file_name_csv = out_put_file_prefix_name + '_solution.csv'
        df_bocop_data.to_pickle(file_name)
        df_bocop_data.to_csv(file_name_csv)
        df_bocop_data.to_pickle(path_pkl)
        df_bocop_data.to_csv(path_csv)
        return df_bocop_data

    def simulation_run_tagger(self):
        prm = self.parameters
        cov_tag = str(prm['coverage'])
        delta_v_tag = str(prm['delta_v'])
        epsilon_tag = str(prm['epsilon'])
        run_tag = 'run_' + \
                  'coverage_' + cov_tag + \
                  'delta_v_' + delta_v_tag + \
                  'epsilon_' + epsilon_tag
        self.run_tag = run_tag

    def reproductive_number(self):
        beta_s = self.parameters['beta_s']
        beta_a = self.parameters['beta_a']
        delta_e = self.parameters['delta_e']
        alpha_s = self.parameters['alpha_s']
        alpha_a = self.parameters['alpha_a']
        lambda_v = self.parameters['lambda_v']
        epsilon = self.parameters['epsilon']
        if self.parameters['delta_v'] == 0:
            delta_v = 0.0
        else:
            delta_v = 1.0 / self.parameters['delta_v']
        mu = self.parameters['mu']
        p = self.parameters['p']

        f_1 = (p * delta_e * beta_s) / ((alpha_s + mu) * (mu + delta_e))
        f_2 = ((1.0 - p) * delta_e * beta_a) / ((mu + alpha_a) * (mu + delta_e))
        r_00 = f_1 + f_2

        fac = 1.0 - epsilon * lambda_v / (mu + delta_v + lambda_v)
        u_v_0 = self.df_controlled_sol['u_V'][0]
        opt_fac = (mu + delta_v + (1.0 - epsilon) * (lambda_v + u_v_0)) / \
                  (mu + delta_v + (lambda_v + u_v_0))

        r_v_0 = fac * r_00
        opt_r_v_0 = opt_fac * r_00
        return r_00, r_v_0, opt_r_v_0

    def sage_scenarios_plot(self):
        def _hex_to_rgb_(hex_color: str) -> tuple:
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
            # i_s = df_solution['i_s']
            y_inc = df_solution['y_inc']
            d = df_solution['d']
            f = a_S * y_inc + a_D * d  # + 0.5 * c_V * u_V ** 2
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

        #
        prm = self.parameters
        lambda_v_base = prm["lambda_v"]
        bocop_sol_path = self.optimal_controlled_solution_path
        constant_vaccination_path = self.constant_controlled_solution_file
        no_vaccination_sol_path = self.uncontrolled_solution_path
        df_oc = pd.read_pickle(bocop_sol_path)
        df_not_vaccination = pd.read_pickle(no_vaccination_sol_path)
        df_constant_vaccination = pd.read_pickle(constant_vaccination_path)
        border_color_pallet = px.colors.sequential.ice
        fill_color_pallet = bokeh_palettes.all_palettes['Category20'][20]
        # Here we start to draw the dynamics.
        ########################################################################
        # Solution without vaccination
        ########################################################################
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
                            ),
            horizontal_spacing=0.17,
            vertical_spacing=0.3,
        )
        n_cdmx = 100000.0
        fig.add_trace(
            go.Scatter(
                x=df_not_vaccination['time'],
                y=n_cdmx * df_not_vaccination['i_s'],
                line=dict(color=border_color_pallet[1], width=1.0, dash='dot'),
                legendgroup='Prevalence',
                name='Without<br>vaccination',
                showlegend=True
            ),
            row=1, col=1
        )
        fig.add_trace(
            go.Scatter(
                x=df_not_vaccination['time'],
                y=0.25 ** 2 * n_cdmx * df_not_vaccination['i_s'],
                line=dict(color=border_color_pallet[1], width=1.0, dash='dot'),
                showlegend=False,
                legendgroup='Optimized resources',
                name='OP'
            ),
            row=1, col=2
        )
        fig.add_trace(
            go.Scatter(
                x=df_not_vaccination['time'],
                y=n_cdmx * df_not_vaccination['i_a'],
                line=dict(color=border_color_pallet[1], width=1.0, dash='dot'),
                legendgroup='Mitigation',
                showlegend=False
            ),
            row=2, col=1
        )
        fig.add_trace(
            go.Scatter(
                x=df_not_vaccination['time'],
                y=n_cdmx * df_not_vaccination['d'],
                line=dict(color=border_color_pallet[1], width=1.0, dash='dot'),
                legendgroup='Saved lives',
                name='(OP)',
                showlegend=False
            ),
            row=1, col=3
        )
        ########################################################################
        # constant_vaccination
        ########################################################################
        trace_constant_vac_i_s = go.Scatter(
            x=df_constant_vaccination['time'],
            y=n_cdmx * df_constant_vaccination['i_s'],
            fill='tonexty',
            fillcolor=f"rgba{(*_hex_to_rgb_(fill_color_pallet[2]), .7)}",
            line=dict(color=fill_color_pallet[2], width=1.0),
            showlegend=True,
            legendgroup='Mitigation',
            name='(CP)'
        )
        trace_constant_vac_i_a = go.Scatter(
            x=df_constant_vaccination['time'],
            y=n_cdmx * df_constant_vaccination['i_a'],
            fill='tonexty',
            fillcolor=f"rgba{(*_hex_to_rgb_(fill_color_pallet[2]), 0.7)}",
            line=dict(color=fill_color_pallet[2], width=1.0),
            showlegend=False
        )
        trace_constant_vac_hospitalization = go.Scatter(
            x=df_constant_vaccination['time'],
            y=0.25 ** 2 * n_cdmx * df_constant_vaccination['i_s'],
            fill='tonexty',
            fillcolor=f"rgba{(*_hex_to_rgb_(fill_color_pallet[6]), 0.7)}",
            line=dict(color=fill_color_pallet[6], width=0.7),
            legendgroup='Optimized resources',
            name='(CP)',
            showlegend=True
        )
        trace_constant_vac_d = \
            go.Scatter(
                x=df_constant_vaccination['time'],
                y=n_cdmx * df_constant_vaccination['d'],
                fillcolor=f"rgba{(*_hex_to_rgb_(fill_color_pallet[4]), 0.8)}",
                fill='tonexty',
                line=dict(color=fill_color_pallet[4], width=1.0),
                legendgroup='Saved lives',
                name='(CP)'
            )
        coverage_t = constant_vaccination_coverage(df_constant_vaccination)
        trace_constant_vac_coverage = go.Scatter(
            x=df_constant_vaccination['time'],
            y=coverage_t,
            fill='tozeroy',
            fillcolor=f"rgba{(*_hex_to_rgb_(fill_color_pallet[18]), 0.7)}",
            line=dict(color=fill_color_pallet[18], width=0.7),
            legendgroup='Coverage',
            name='(CP)',
            showlegend=False
        )
        cost_t = constant_vaccination_cost(df_constant_vaccination)
        trace_constant_vac_cost = go.Scatter(
            x=df_constant_vaccination['time'],
            y=cost_t,
            fill='tonexty',
            fillcolor=f"rgba{(*_hex_to_rgb_(fill_color_pallet[16]), 0.7)}",
            line=dict(color=fill_color_pallet[16], width=1.0),
            legendgroup='Cost',
            name='(CP)',
            showlegend=True
        )
        #######################################################################
        # Optimal solution trace
        #######################################################################
        trace_optimal_prevalence = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * df_oc['i_s'],
            line=dict(color=border_color_pallet[1],
                      width=1.0,
                      dash='solid'),
            legendgroup='Prevalence',
            name='vaccination',
            showlegend=True
        )
        trace_optimal_i_s = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * df_oc['i_s'],
            fill='tonexty',
            fillcolor=f"rgba{(*_hex_to_rgb_(fill_color_pallet[3]), 0.7)}",
            legendgroup='Mitigation',
            line=dict(color=border_color_pallet[1], width=.7),
            name='(OP)',
            showlegend=True
        )
        #
        trace_optimal_hospitalization = go.Scatter(
            x=df_oc['time'],
            y=0.25 ** 2 * n_cdmx * df_oc['i_s'],
            fill='tonexty',
            fillcolor=f"rgba{(*_hex_to_rgb_(fill_color_pallet[7]), 0.7)}",
            line=dict(color=border_color_pallet[1], width=.7),
            legendgroup='Optimized resources',
            showlegend=True,
            name='(OP)'
        )
        trace_optimal_d = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * df_oc['d'],
            fillcolor=f"rgba{(*_hex_to_rgb_(fill_color_pallet[5]), 0.4)}",
            fill='tonexty',
            legendgroup='Saved lives',
            line=dict(color=border_color_pallet[1], width=1),
            showlegend=True,
            name='(OP)'
        )
        #
        trace_optimal_i_a = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * df_oc['i_a'],
            fillcolor=f"rgba{(*_hex_to_rgb_(fill_color_pallet[3]), 0.7)}",
            fill='tonexty',
            line=dict(color=border_color_pallet[1], width=.7),
            showlegend=False
        )
        trace_optimal_vac_coverage = go.Scatter(
            x=df_oc['time'],
            y=df_oc['x_vac'],
            line=dict(color=fill_color_pallet[19], width=.7),
            fill='tozeroy',
            fillcolor=f"rgba{(*_hex_to_rgb_(fill_color_pallet[19]), 0.5)}",
            legendgroup='Coverage',
            name="(OP)",
            showlegend=False
        )
        #
        lambda_base_v_t = lambda_v_base * np.ones(df_oc["u_V"].shape[0])
        trace_optimal_vaccination_policy = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * (df_oc['u_V'] + lambda_base_v_t),
            fill='tozeroy',
            fillcolor=f"rgba{(*_hex_to_rgb_(fill_color_pallet[19]), 0.5)}",
            line=dict(color=fill_color_pallet[0], width=0.7),
            name="(OP)",
            legendgroup='Policy',
            showlegend=False
        )
        trace_optimal_cost = go.Scatter(
            x=df_oc['time'],
            y=df_oc['x_zero'],
            line=dict(color=fill_color_pallet[16], width=1.0),
            fillcolor=f"rgba{(*_hex_to_rgb_(fill_color_pallet[17]), 0.4)}",
            fill='tozeroy',
            legendgroup='Cost',
            name='(OP)',
            showlegend=True
        )
        #
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
        fig.append_trace(trace_constant_vac_i_s, 1, 1)
        fig.append_trace(trace_optimal_i_s, 1, 1)
        fig.append_trace(trace_optimal_prevalence, 1, 1)
        fig.append_trace(trace_constant_vac_hospitalization, 1, 2)
        fig.append_trace(trace_optimal_hospitalization, 1, 2)
        fig.append_trace(trace_constant_vac_i_a, 2, 1)
        fig.append_trace(trace_optimal_i_a, 2, 1)
        fig.append_trace(trace_constant_vac_coverage, 2, 2)
        fig.append_trace(trace_optimal_vac_coverage, 2, 2)
        fig.append_trace(trace_constant_vac_d, 1, 3)
        fig.append_trace(trace_optimal_d, 1, 3)
        fig.append_trace(trace_vaccination_base, 3, 1)
        fig.append_trace(trace_vaccination_max, 3, 1)
        fig.append_trace(trace_optimal_vaccination_policy, 3, 1)
        fig.append_trace(trace_optimal_cost, 3, 3)
        fig.append_trace(trace_constant_vac_cost, 3, 3)
        #
        # Axis labels
        #
        fig.add_annotation(
            dict(
                text='Cases per 100,000 inhabitants',
                align='left',
                textangle=-90,
                font=dict(family="Arial", size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=-0.065,
                y=1.01)
        )
        fig.update_xaxes(title_text="Time (days)",
                         title_font=dict(size=10, family='Courier'),
                         tickfont=dict(size=10, family='Courier'),
                         row=1, col=3
                         )
        fig.update_xaxes(title_text="Time (days)",
                         title_font=dict(size=10, family='Courier'),
                         tickfont=dict(size=10, family='Courier'),
                         row=2, col=1)
        fig.update_xaxes(title_text="Time (days)",
                         title_font=dict(size=10, family='Courier'),
                         tickfont=dict(size=10, family='Courier'),
                         row=2, col=2)
        fig.update_xaxes(title_text="Time (days)",
                         title_font=dict(size=10, family='Courier'),
                         tickfont=dict(size=10, family='Courier'),
                         row=3, col=1)
        fig.update_xaxes(title_text="Time (days)",
                         title_font=dict(size=10, family='Courier'),
                         tickfont=dict(size=10, family='Courier'),
                         row=3, col=3)
        #
        fig.update_yaxes(tickfont=dict(size=10, family='Courier'),
                         row=1, col=1)
        fig.update_yaxes(tickfont=dict(size=10, family='Courier'),
                         row=1, col=2)
        fig.add_annotation(
            dict(
                text="Cases (per 100,000)",
                align='left',
                textangle=-90,
                font=dict(family="Arial", size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=0.735,
                y=.95
            )
        )
        fig.update_yaxes(tickfont=dict(size=10, family='Courier'),
                         row=1, col=3)
        fig.update_yaxes(tickfont=dict(size=10, family='Courier'),
                         row=2, col=1)
        fig.add_annotation(
            dict(
                text="Population<br>fraction",
                align='left',
                textangle=-90,
                font=dict(family="Arial", size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=0.275,
                y=.50
            )
        )
        fig.update_yaxes(tickfont=dict(size=10, family='Courier'),
                         row=2, col=2)
        fig.add_annotation(
            dict(
                text="Doses per<br>100,000<br>inhabitants",
                align='left',
                textangle=-90,
                font=dict(family="Arial", size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=-0.15,
                y=-0.1)
        )
        fig.update_yaxes(tickfont=dict(size=10, family='Courier'),
                         row=3, col=1)
        fig.add_annotation(
            dict(
                text="Cost per<br>100,000<br>inhabitants",
                align='left',
                textangle=-90,
                font=dict(family="Arial", size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=0.74,
                y=-0.1)
        )
        #
        fig.update_yaxes(tickfont=dict(size=10, family='Courier'),
                         row=3, col=3)
        #
        r_0, r_v_0, r_opt_v_0 = self.reproductive_number()
        c_r_0 = decimal.Decimal(r_0)
        c_r_v_0 = decimal.Decimal(r_v_0)
        c_r_opt_v_0 = decimal.Decimal(r_opt_v_0)
        #
        # Legend annotations
        # ----------------------------------------------------------------------
        # Mitigation
        # ----------------------------------------------------------------------
        fig.add_annotation(
            dict(
                text='Mitigation',
                align='left',
                font=dict(family="Arial",
                          size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=1.21,
                y=0.78)
        )
        # ----------------------------------------------------------------------
        # Saved Beds
        # ----------------------------------------------------------------------
        fig.add_annotation(
            dict(
                text='Beds',
                align='left',
                font=dict(family="Arial",
                          size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=1.17,
                y=0.56)
        )
        # ----------------------------------------------------------------------
        # Coverage
        # ----------------------------------------------------------------------
        # fig.add_annotation(
        #     dict(
        #         text='Coverage',
        #         align='left',
        #         font=dict(family="Arial", size=10),
        #         showarrow=False,
        #         xref='paper',
        #         yref='paper',
        #         x=1.11,
        #         y=.27)
        # )
        # ----------------------------------------------------------------------
        # Saved Lives
        # ----------------------------------------------------------------------
        fig.add_annotation(
            dict(
                text='Saved Lives',
                align='left',
                font=dict(family="Arial",
                          size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=1.23,
                y=0.385
            )
        )
        fig.add_annotation(
            dict(
                text='Cost',
                align='left',
                font=dict(family="Arial", size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=1.165,
                y=0.17)
        )
        data_label = {'eps': prm["epsilon"],
                      'delta_v': round(prm["delta_v"], 1),
                      'time_unit': 'days'}
        str_vaccination_par = \
            r'$\epsilon={:1.2f}, \quad \delta_V ^{{-1}}={:1.1f} \ \mathtt{{{' \
            r':>5}}}$'.format(
                data_label['eps'],
                data_label['delta_v'],
                data_label['time_unit'])
        fig.add_annotation(
            dict(
                text=str_vaccination_par,
                align='left',
                font=dict(family="Arial",
                          size=12),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=1.085,
                y=1.24,
                bordercolor="Black",
                borderwidth=1
            )
        )
        fig.update_layout(
            font=dict(family="Arial", size=10),
            template="simple_white",
            showlegend=True,
            title_text="R0: " + str(round(c_r_0, 6))
                       + '\t\tRv: ' + str(round(c_r_v_0, 6)),
            legend=dict(xanchor='left',
                        yanchor='top',
                        x=1.1,
                        y=1.0,
                        # traceorder="normal",
                        bordercolor="Black",
                        borderwidth=2
                        )
        )
        # ----------------------------------------------------------------------
        # figure alpha number
        # ----------------------------------------------------------------------
        fig.add_annotation(
            dict(
                text='<b>A)',
                align='left',
                font=dict(family="Arial",
                          size=14),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=-.034,
                y=1.06,
            )
        )
        fig.add_annotation(
            dict(
                text='<b>B)',
                align='left',
                font=dict(family="Arial",
                          size=14),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=0.37,
                y=1.06,
            )
        )
        fig.add_annotation(
            dict(
                text='<b>C)',
                align='left',
                font=dict(family="Arial",
                          size=14),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=0.785,
                y=1.06,
            )
        )
        fig.add_annotation(
            dict(
                text='<b>D)',
                align='left',
                font=dict(family="Arial",
                          size=14),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=-0.032,
                y=0.6,
            )
        )
        fig.add_annotation(
            dict(
                text='<b>E)',
                align='left',
                font=dict(family="Arial",
                          size=14),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=0.37,
                y=.6,
            )
        )
        fig.add_annotation(
            dict(
                text='<b>F)',
                align='left',
                font=dict(family="Arial",
                          size=14),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=-0.032,
                y=.14,
            )
        )
        fig.add_annotation(
            dict(
                text='<b>G)',
                align='left',
                font=dict(family="Arial",
                          size=14),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=0.785,
                y=.14,
            )
        )
        for i in fig['layout']['annotations']:
            i['font'] = dict(family="Arial", size=10)
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
        # fig.show()

    # ------------------------------------------------------------------------------
    # Accorded figures
    #

    def constant_vaccination_cost(self, df_solution):
        """

        """
        prm = self.parameters
        a_D = prm['a_D']
        a_S = prm['a_S']
        c_V = prm['c_V']
        u_V = prm['lambda_v']
        t = df_solution['time']
        i_s = df_solution['i_s']
        d = df_solution['d']
        f = a_S * i_s + a_D * d + 0.5 * c_V * u_V ** 2
        cost_t_ = integrate.cumtrapz(f, t, initial=0)
        return cost_t_

    def constant_vaccination_coverage(self, df_solution):
        """

        """
        prm = self.parameters
        lambda_v = prm['lambda_v']
        t = df_solution['time']
        s = df_solution['s']
        e = df_solution['e']
        i_a = df_solution['i_a']
        r = df_solution['r']
        f = lambda_v * (s + e + i_a + r)
        coverage_t_ = integrate.cumtrapz(f, t, initial=0)
        return coverage_t_

    def daly_vaccination_cost(self, df_solution) -> float:
        """
            Units in DALY, see
            https://www.who.int/healthinfo/global_burden_disease/metrics_daly/en/
        """
        prm = self.parameters
        mu = prm['mu']
        delta_e = prm['delta_e']
        alpha_s = prm['alpha_s']
        mu_s = prm['mu_s']
        a_D = prm['a_D']
        a_E = prm['a_E']
        t = df_solution['time']
        # TODO: Chane label for incidence
        y_inc = df_solution['y_inc']
        d = df_solution['d']
        d_0 = prm['D_0']
        y_inc_0 = prm['Y_inc_0']
        n_cdmx = 100000
        life_expectancy = mu ** (-1.0) / 365.0
        disability_weight = 0.28425
        # see Jo 2020 for reference and regarding values for our mean
        # calculation
        # TODO: fix daly function
        case_average_duration = alpha_s ** (-1.0) / 365
        years_of_life_lost = a_D * (d - d_0)
        years_of_life_disability = a_E * (y_inc - y_inc_0)
        daly = n_cdmx * (years_of_life_lost + years_of_life_disability)
        return daly

    def sage_plot_prevalence(self, fig_file_name_prefix='prevalence_fig_'):

        prm = self.parameters
        optimal_controlled_sol_path = self.optimal_controlled_solution_path
        uncontrolled_solution_path = self.uncontrolled_solution_path
        constant_controlled_solution_file = \
            self.constant_controlled_solution_file
        df_oc = pd.read_pickle(optimal_controlled_sol_path)
        df_not_vaccination = pd.read_pickle(uncontrolled_solution_path)
        df_constant_vaccination = \
            pd.read_pickle(constant_controlled_solution_file)
        border_color_pallet = px.colors.sequential.ice
        fill_color_pallet = bokeh_palettes.all_palettes['Category20'][20]
        prevalence_fig = make_subplots(
            rows=1, cols=2,
            specs=[
                [{}, {}]
            ],
            subplot_titles=("Symptomatic",
                            "Death"
                            ),
            horizontal_spacing=0.17
        )
        n_cdmx = 100000
        prevalence_fig.add_trace(
            go.Scatter(
                x=df_not_vaccination['time'],
                y=n_cdmx * df_not_vaccination['i_s'],
                line=dict(color=border_color_pallet[1], width=.7, dash='dot'),
                fill='tozeroy',
                fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[6]), 1.0)}",
                legendgroup='Prevalence',
                name='Without<br>vaccination',
                showlegend=True
            ),
            row=1, col=1
        )
        prevalence_fig.add_trace(
            go.Scatter(
                x=df_not_vaccination['time'],
                y=n_cdmx * df_not_vaccination['d'],
                line=dict(color=border_color_pallet[1], width=.7, dash='dot'),
                legendgroup='Saved lives',
                name='(OP)',
                showlegend=False
            ),
            row=1, col=2
        )
        ########################################################################
        # constant_vaccination
        ########################################################################
        trace_constant_vac_i_s = go.Scatter(
            x=df_constant_vaccination['time'],
            y=n_cdmx * df_constant_vaccination['i_s'],
            fill='tozeroy',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[0]), .7)}",
            line=dict(color=border_color_pallet[0], width=.7),
            showlegend=True,
            legendgroup='Mitigation',
            name='(CP)'
        )
        trace_constant_vac_d = \
            go.Scatter(
                x=df_constant_vaccination['time'],
                y=n_cdmx * df_constant_vaccination['d'],
                fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[4]), 0.8)}",
                fill='tonexty',
                line=dict(color=fill_color_pallet[4], width=.7),
                legendgroup='Saved lives',
                name='(CP)'
            )
        #######################################################################
        # Optimal solution trace
        #######################################################################
        trace_optimal_prevalence = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * df_oc['i_s'],
            line=dict(color=border_color_pallet[1],
                      width=.7,
                      dash='solid'),
            legendgroup='Prevalence',
            name='vaccination',
            showlegend=True
        )
        trace_optimal_i_s = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * df_oc['i_s'],
            fill='tozeroy',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[3]), 0.5)}",
            legendgroup='Mitigation',
            line=dict(color=border_color_pallet[1], width=.7),
            name='(OP)',
            showlegend=True
        )
        #
        trace_optimal_d = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * df_oc['d'],
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[5]), 0.4)}",
            fill='tonexty',
            legendgroup='Saved lives',
            line=dict(color=border_color_pallet[1], width=1),
            showlegend=True,
            name='(OP)'
        )

        prevalence_fig.append_trace(trace_constant_vac_i_s, 1, 1)
        prevalence_fig.append_trace(trace_optimal_i_s, 1, 1)
        prevalence_fig.append_trace(trace_optimal_prevalence, 1, 1)
        prevalence_fig.append_trace(trace_constant_vac_d, 1, 2)
        prevalence_fig.append_trace(trace_optimal_d, 1, 2)
        #
        # Axis labels
        #
        prevalence_fig.add_annotation(
            dict(
                text='Cases per 100,000 inhabitants',
                align='left',
                textangle=-90,
                font=dict(family="Arial", size=11),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=-0.095,
                y=0.8)
        )
        prevalence_fig.update_xaxes(
            title_text="Time (days)",
            title_font=dict(size=10, family='Arial'),
            tickfont=dict(size=10, family='Arial'),
            row=1, col=1
        )
        prevalence_fig.update_xaxes(
            title_text="Time (days)",
            title_font=dict(size=10, family='Arial'),
            tickfont=dict(size=10, family='Arial'),
            row=1, col=2
        )
        #
        prevalence_fig.update_yaxes(
            tickfont=dict(size=10, family='Arial'),
            row=1, col=1
        )
        prevalence_fig.update_yaxes(
            tickfont=dict(size=10, family='Arial'),
            row=1, col=2
        )
        prevalence_fig.update_yaxes(
            tickfont=dict(size=10, family='Arial'),
            row=1, col=1
        )
        prevalence_fig.update_yaxes(
            tickfont=dict(size=10, family='Arial'),
            row=1, col=2
        )
        #
        #
        #
        #
        r_0, r_v_0, r_opt_v_0 = self.reproductive_number()
        c_r_0 = decimal.Decimal(r_0)
        c_r_v_0 = decimal.Decimal(r_v_0)
        c_r_opt_v_0 = decimal.Decimal(r_opt_v_0)
        #
        # Legend annotations
        # ----------------------------------------------------------------------
        # Mitigation
        # ----------------------------------------------------------------------
        prevalence_fig.add_annotation(
            dict(
                text='Vaccination',
                align='left',
                font=dict(family="Arial",
                          size=10),
                showarrow=False,
                xref='paper',
                yref=''
                     'paper',
                x=1.23,
                y=0.82)
        )
        # ----------------------------------------------------------------------
        # Saved Beds
        # ----------------------------------------------------------------------
        prevalence_fig.add_annotation(
            dict(
                text='Saved Lives',
                align='left',
                font=dict(family="Arial",
                          size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=1.23,
                y=0.67
            )
        )
        delta_v_label = round(prm["delta_v"], 2)
        if delta_v_label == 0.0:
            delta_v_label = 'lifelong'
        data_label = {'eps': prm["epsilon"],
                      'delta_v': delta_v_label,
                      'time_unit': 'days'}

        if delta_v_label == 'lifelong':
            str_vaccination_par = \
                r'$\epsilon={:1.2f}, \quad 1 / \delta_V= \mathsf{{{' \
                r':>9}}}$'.format(
                    data_label['eps'],
                    data_label['delta_v'])
        else:
            str_vaccination_par = \
                r'$\epsilon={:1.2f}, \quad 1 / \delta_V={:1.1f} \ \mathsf{{{' \
                r':>5}}}$'.format(
                    data_label['eps'],
                    data_label['delta_v'],
                    data_label['time_unit'])
        prevalence_fig.add_annotation(
            dict(
                text=str_vaccination_par,
                align='left',
                font=dict(family="Arial",
                          size=12),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=1.085,
                y=1.24,
                bordercolor="Black",
                borderwidth=1
            )
        )
        str_policy_legend = '<b>Vaccination Policy:</b><br>' + \
                            '    Optimal (OP)<br>' + \
                            '    Constant (CP)'
        prevalence_fig.add_annotation(
            dict(
                text=str_policy_legend,
                align='left',
                font=dict(family="Arial",
                          size=12),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=1.295,
                y=0.05
            )
        )
        prevalence_fig.update_layout(
            font=dict(family="Arial", size=10),
            template="simple_white",
            showlegend=True,
            title_text="R0: " + str(round(c_r_0, 6))
                       + '\t\tRv: ' + str(round(c_r_v_0, 6)),
            legend=dict(xanchor='left',
                        yanchor='top',
                        x=1.1,
                        y=1.0,
                        # traceorder="normal",
                        bordercolor="Black",
                        borderwidth=2
                        )
        )
        # ----------------------------------------------------------------------
        # prevalence_figure alpha number
        # ----------------------------------------------------------------------
        prevalence_fig.add_annotation(
            dict(
                text='<b>A)',
                align='left',
                font=dict(family="Arial",
                          size=14),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=-.034,
                y=1.06,
            )
        )
        prevalence_fig.add_annotation(
            dict(
                text='<b>B)',
                align='left',
                font=dict(family="Arial",
                          size=14),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=0.57,
                y=1.06,
            )
        )
        for i in prevalence_fig['layout']['annotations']:
            i['font'] = dict(family="Arial", size=10)
        if not os.path.exists("images"):
            os.mkdir("images")
        # prevalence_fig.write_image("images/prevalence_fig1.pdf")
        golden_width = 718  # width in px
        golden_ratio = 1.618
        # pio.kaleido.scope.default_format = "pdf"
        # pio.kaleido.scope.default_width = golden_width
        # pio.kaleido.scope.default_height = golden_width / golden_ratio
        # pio.kaleido.scope.default_scale = .50
        # prevalence_fig.to_image(format="pdf", engine="kaleido")
        # prevalence_fig.write_image("images/prevalence_fig.pdf")
        # TODO:  Edit legend respect to groups
        time_now = datetime.now()
        dt_string = time_now.strftime("%b-%d-%Y_%H_")
        run_tag = self.run_tag
        path_fig = fig_file_name_prefix + \
                   '_' + \
                   dt_string + \
                   run_tag + '.pdf'
        # # print(path_fig)
        pio.write_image(prevalence_fig, "./images/" + path_fig)
        pio.write_image(prevalence_fig, "./images/prevalence_fig.pdf")
        pio.write_image(prevalence_fig, "./images/prevalence_fig.png")
        # prevalence_fig.show()

    # --------------------------------------------------------------------------
    # Figure 02
    # --------------------------------------------------------------------------

    def sage_plot_optimal_signal(self,
                                 fig_file_name_prefix="optimal_signal_fig"):
        n_cdmx = 100000
        prm = self.parameters
        lambda_v_base = prm["lambda_v"]
        uncontrolled_solution_path = self.uncontrolled_solution_path
        optimal_controlled_sol_path = self.optimal_controlled_solution_path
        constant_controlled_solution_file = \
            self.constant_controlled_solution_file
        df_not_vaccination = pd.read_pickle(uncontrolled_solution_path)
        df_oc = pd.read_pickle(optimal_controlled_sol_path)
        df_constant_vaccination = pd.read_pickle(
            constant_controlled_solution_file)
        fill_color_pallet = bokeh_palettes.all_palettes['Category20c'][20]
        optimal_signal_fig = make_subplots(
            rows=4, cols=2,
            specs=[
                [{"rowspan": 2}, {"rowspan": 2}],
                [None, None],
                [{"rowspan": 2, "colspan": 2}, None],
                [None, None]
            ],
            subplot_titles=(
                "Coverage",
                "Cost",
                "Vaccination Policy"
            ),
            horizontal_spacing=0.17,
            vertical_spacing=0.3,
        )
        ########################################################################
        # constant_vaccination
        ########################################################################
        coverage_t = self.constant_vaccination_coverage(df_constant_vaccination)
        trace_constant_vac_coverage = go.Scatter(
            x=df_constant_vaccination['time'],
            y=coverage_t,
            fill='tozeroy',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[13]), 1.0)}",
            line=dict(color=fill_color_pallet[12], width=0.7),
            legendgroup='Coverage',
            name='(CP) Coverage',
            showlegend=False
        )
        trace_optimal_vac_coverage = go.Scatter(
            x=df_oc['time'],
            y=df_oc['x_vac'],
            line=dict(color=fill_color_pallet[16], width=.7),
            fill='tozeroy',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[19]), 0.7)}",
            legendgroup='Coverage',
            name="(OP) Coverage",
            showlegend=False
        )
        #
        no_vaccination_cost_t = self.daly_vaccination_cost(df_not_vaccination)
        constant_cost_t = self.daly_vaccination_cost(df_constant_vaccination)
        self.daly_constant_cost = constant_cost_t.values[-1]
        optimal_cost_t = self.daly_vaccination_cost(df_oc)
        self.daly_optimal_cost = optimal_cost_t.values[-1]
        trace_no_vac_cost = go.Scatter(
            x=df_not_vaccination['time'],
            y=no_vaccination_cost_t,
            fill='tonexty',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[9]), 0.7)}",
            line=dict(color=fill_color_pallet[8], width=1.0),
            legendgroup='Cost',
            name='No vaccination',
            showlegend=True
        )
        trace_constant_vac_cost = go.Scatter(
            x=df_constant_vaccination['time'],
            y=constant_cost_t,
            fill='tonexty',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[13]), 0.7)}",
            line=dict(color=fill_color_pallet[12], width=1.0),
            legendgroup='Cost',
            name='(CP)',
            showlegend=True
        )
        trace_optimal_cost = go.Scatter(
            x=df_oc['time'],
            y=optimal_cost_t,
            line=dict(color=fill_color_pallet[16], width=1.0),
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[19]), 1.0)}",
            fill='tozeroy',
            legendgroup='Cost',
            name='(OP)',
            showlegend=True
        )
        #######################################################################
        # Optimal solution trace
        #######################################################################
        opt_x_vac_prevalence = df_oc['s'] + df_oc['e'] \
                               + df_oc['i_a'] + df_oc['r']
        cst_x_vac_prevalence = df_constant_vaccination['s'] + \
                                df_constant_vaccination['e'] + \
                               df_constant_vaccination['i_a'] + \
                               df_constant_vaccination['r']

        lambda_base_v_t = lambda_v_base * np.ones(df_oc["u_V"].shape[0])
        trace_optimal_vaccination_policy = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * (df_oc['u_V'] + lambda_base_v_t),
            fill='tozeroy',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[19]), 0.5)}",
            line=dict(color=fill_color_pallet[0], width=0.5),
            name="Optimal<br>vaccination<br>policy",
            legendgroup='Policy',
            showlegend=True
        )
        trace_vaccination_base = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * lambda_base_v_t,
            fill='tozeroy',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[19]), 0.5)}",
            line=dict(color=fill_color_pallet[0], width=0.7),
            name="Optimal<br>vaccination<br>policy",
            legendgroup='Policy',
            showlegend=True
        )
        #
        trace_vaccination_base = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * lambda_v_base * np.ones(df_oc["u_V"].shape[0]),
            line=dict(color=fill_color_pallet[18], width=0.7, dash='dot'),
            name='Constant<br>policy',
            showlegend=True
        )
        trace_optimal_vaccination_action = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * (df_oc['u_V'] + lambda_base_v_t) * opt_x_vac_prevalence,
            fill='tozeroy',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[19]), 0.5)}",
            line=dict(color=fill_color_pallet[16], width=0.7),
            name="Optimal<br>vaccination<br>policy",
            legendgroup='Policy',
            showlegend=True
        )
        trace_vaccination_base_action = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * lambda_v_base * cst_x_vac_prevalence,
            line=dict(color=fill_color_pallet[0], width=0.7, dash='dot'),
            fill='tozeroy',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[0]), 0.2)}",
            name='Constant<br>policy',
            showlegend=True
        )
        #
        unitary_vaccine_price = 6.1942
        base_vaccine_price = integrate.cumtrapz(
            x=df_oc['time'],
            y=(n_cdmx * lambda_v_base * unitary_vaccine_price) * \
              cst_x_vac_prevalence
        )
        total_base_vaccination_price = \
            round(decimal.Decimal(base_vaccine_price[-1]), 0)
        opt_vaccine_price = integrate.cumtrapz(
            x=df_oc['time'],
            y=n_cdmx * unitary_vaccine_price * \
              (df_oc['u_V'] + lambda_base_v_t) * opt_x_vac_prevalence
        )
        total_opt_vaccination_price = \
            round(decimal.Decimal(opt_vaccine_price[-1]), 0)
        str_prices = '<b>Cost:</b><br>' + \
                    '(CP): ' + str(total_base_vaccination_price) + \
                    ' dls<br>' + \
                    '(OP): ' + str(total_opt_vaccination_price) + \
                    ' dls<br>' +\
                    'price/dose: 6.20 dls'
        #
        optimal_signal_fig.add_annotation(
            dict(
                text=str_prices,
                align='left',
                font=dict(family="Arial",
                          size=12,
                          color=fill_color_pallet[16]
                          ),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=1.3,
                y=0.01)
        )
#
#
        optimal_signal_fig.append_trace(trace_constant_vac_coverage, 1, 1)
        optimal_signal_fig.append_trace(trace_optimal_vac_coverage, 1, 1)
        optimal_signal_fig.append_trace(trace_no_vac_cost, 1, 2)
        optimal_signal_fig.append_trace(trace_optimal_cost, 1, 2)
        optimal_signal_fig.append_trace(trace_constant_vac_cost, 1, 2)
        optimal_signal_fig.append_trace(trace_vaccination_base_action, 3, 1)
        optimal_signal_fig.append_trace(trace_optimal_vaccination_action, 3, 1)
        # optimal_signal_fig.append_trace(trace_vaccination_base, 3, 1)
        # optimal_signal_fig.
        # append_trace(trace_optimal_vaccination_policy, 3, 1)
        #
        # Axis labels
        #
        optimal_signal_fig.add_annotation(
            dict(
                text='Population proportion',
                align='left',
                textangle=-90,
                font=dict(family="Arial", size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=-0.09,
                y=1.03)
        )
        optimal_signal_fig.update_xaxes(
            title_text="Time (days)",
            title_font=dict(size=10, family='Arial'),
            tickfont=dict(size=10, family='Arial'),
            row=1, col=1
        )
        optimal_signal_fig.update_xaxes(
            title_text="Time (days)",
            title_font=dict(size=10, family='Arial'),
            tickfont=dict(size=10, family='Arial'),
            row=1, col=2
        )
        optimal_signal_fig.update_xaxes(
            title_text="Time (days)",
            title_font=dict(size=10, family='Arial'),
            tickfont=dict(size=10, family='Arial'),
            row=3, col=1
        )
        #
        optimal_signal_fig.update_yaxes(
            tickfont=dict(size=10, family='Arial'),
            row=1, col=1
        )
        optimal_signal_fig.update_yaxes(
            tickfont=dict(size=10, family='Arial'),
            row=1, col=2
        )
        #
        optimal_signal_fig.update_yaxes(
            tickfont=dict(size=10, family='Arial'),
            row=1, col=3)
        optimal_signal_fig.update_yaxes(
            tickfont=dict(size=10, family='Arial'),
            row=2, col=1
        )

        optimal_signal_fig.update_yaxes(
            tickfont=dict(size=10, family='Arial'),
            row=2, col=2
        )
        optimal_signal_fig.add_annotation(
            dict(
                text="Doses per<br>100,000<br>inhabitants",
                align='left',
                textangle=-90,
                font=dict(family="Arial", size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=-0.15,
                y=0.07)
        )
        optimal_signal_fig.update_yaxes(
            tickfont=dict(size=10, family='Arial'),
            row=3, col=1
        )
        optimal_signal_fig.add_annotation(
            dict(
                # text="Cost per 100,000 inhabitants",
                text="DALYs per<br>100,000<br>inhabitants",
                align='left',
                textangle=-90,
                font=dict(family="Arial", size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=0.47,
                # y=1.05,
                y=.95
            )
        )
        #
        optimal_signal_fig.update_yaxes(tickfont=dict(size=10, family='Arial'),
                                        row=3, col=3)
        #
        r_0, r_v_0, r_opt_v_0 = self.reproductive_number()
        c_r_0 = decimal.Decimal(r_0)
        c_r_v_0 = decimal.Decimal(r_v_0)
        c_r_opt_v_0 = decimal.Decimal(r_opt_v_0)
        delta_v_label = round(prm["delta_v"], 2)
        if delta_v_label == 0.0:
            delta_v_label = 'lifelong'
        data_label = {'eps': prm["epsilon"],
                      'delta_v': delta_v_label,
                      'time_unit': 'days'}

        if delta_v_label == 'lifelong':
            str_vaccination_par = \
                r'$\epsilon={:1.2f}, \quad 1 / \delta_V= \mathsf{{{' \
                r':>9}}}$'.format(
                    data_label['eps'],
                    data_label['delta_v'])
        else:
            str_vaccination_par = \
                r'$\epsilon={:1.2f}, \quad 1 / \delta_V={:1.1f} \ \mathsf{{{' \
                r':>5}}}$'.format(
                    data_label['eps'],
                    data_label['delta_v'],
                    data_label['time_unit'])
        str_policy_legend = '<b>Vaccination Policy:</b><br>' + \
                            '    Optimal (OP)<br>' + \
                            '    Constant (CP)'
        optimal_signal_fig.add_annotation(
            dict(
                text=str_policy_legend,
                align='left',
                font=dict(family="Arial",
                          size=12),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=1.325,
                y=0.25
            )
        )
        optimal_signal_fig.add_annotation(
            dict(
                text=str_vaccination_par,
                align='left',
                font=dict(family="Arial",
                          size=12),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=1.085,
                y=1.24,
                bordercolor="Black",
                borderwidth=1
            )
        )
        optimal_signal_fig.update_layout(
            font=dict(family="Arial", size=10),
            template="simple_white",
            showlegend=True,
            title_text="R0: " + str(round(c_r_0, 6))
                       + '\t\tRv: ' + str(round(c_r_v_0, 6)),
            legend=dict(xanchor='left',
                        yanchor='top',
                        x=1.1,
                        y=1.0,
                        # traceorder="normal",
                        bordercolor="Black",
                        borderwidth=2
                        )
        )
        # ----------------------------------------------------------------------
        # optimal_signal_figure alpha number
        # ----------------------------------------------------------------------
        optimal_signal_fig.add_annotation(
            dict(
                text='<b>A)',
                align='left',
                font=dict(family="Arial",
                          size=14),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=-.034,
                y=1.06,
            )
        )
        optimal_signal_fig.add_annotation(
            dict(
                text='<b>B)',
                align='left',
                font=dict(family="Arial",
                          size=14),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=0.57,
                y=1.06,
            )
        )
        #
        optimal_signal_fig.add_annotation(
            dict(
                text='<b>C)',
                align='left',
                font=dict(family="Arial",
                          size=14),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=-0.032,
                y=.4,
            )
        )
        #
        for i in optimal_signal_fig['layout']['annotations']:
            i['font'] = dict(family="Arial", size=10)
        if not os.path.exists("images"):
            os.mkdir("images")
        golden_width = 718  # width in px
        golden_ratio = 1.618
        time_now = datetime.now()
        dt_string = time_now.strftime("%b-%d-%Y_%H_")
        run_tag = self.run_tag
        path_fig = fig_file_name_prefix + \
                   '_' + \
                   dt_string + \
                   run_tag + '.png'
        path_fig_pdf = fig_file_name_prefix + \
                   '_' + \
                   dt_string + \
                   run_tag + '.pdf'
        pio.write_image(optimal_signal_fig, "images/" + path_fig)
        pio.write_image(optimal_signal_fig, "images/" + path_fig_pdf)
        pio.write_image(optimal_signal_fig, "images/optimal_signal_fig.pdf")
        pio.write_image(optimal_signal_fig, "images/optimal_signal_fig.png")
