import numpy as np
from scipy.integrate import solve_ivp
import json
import re
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime

bocop_main_scene_data = 'bocop_data/'\
                  + 'deltaR_05_Coverage_50_TH_10M_14M/delta_V_05_years/' \
                  + 'epsilon-7.000000e-01.sol'


class NumericsCovid19:
    def __init__(self,
                 uncontrolled_solution_path='not_vaccination_solution.pkl',
                 optimal_controlled_solution='bocop_solution.pkl',
                 bocop_parameters_json_file='bocop_run_parameters' +
                                            '.json',
                 vaccination_parameters_json_file='vaccination_parameters' +
                                                  '.json',
                 bocop_solution_file=bocop_main_scene_data,
                 n_steps=2000):
        self.constant_control = False
        with open(bocop_parameters_json_file) as json_file:
            bocop_prm_ = json.load(json_file)
        with open(vaccination_parameters_json_file) as json_file:
            vaccination_prm_ = json.load(json_file)
        self.parameters = vaccination_prm_
        self.n_steps = n_steps
        self.time_state_solution_names = ['time', 's', 'e',
                                          'i_s', 'i_a',
                                          'r', 'd', 'v',
                                          'x_vac']

        # TODO: 'cost'
        self.state_dim = self.time_state_solution_names.__len__()
        self.t_eval = np.linspace(0, self.parameters['T'],
                                  num=self.n_steps,
                                  endpoint=True)
        self.x_sol_vaccine = np.zeros([self.n_steps, self.state_dim])
        self.x_sol_not_vaccine = np.zeros([self.n_steps, self.state_dim])
        self.df_sol_cte_vaccine = pd.DataFrame(self.x_sol_vaccine)
        self.df_sol_not_vaccine = pd.DataFrame(self.x_sol_not_vaccine)
        self.u_v = np.zeros([self.state_dim, 1])
        self.vaccination_parameters_json_file = vaccination_parameters_json_file
        self.bocop_solution_file = bocop_solution_file
        self.state_keys = []
        self.control_keys = []
        self.parameter_keys = []
        self.boundarycond_keys = []
        self.constraint_keys = []
        self.control_keys = []

    def read_bocop_parameters_keys(self, parameters_prefix_file_name):
        # TODO: Recover parameters key from bocop_run_parameters.json
        star_pointer_block = 0
        end_pointer_block = 0
        state_keys = []
        control_keys = []
        parameter_keys = []
        boundarycond_keys = []
        constraint_keys = []
        constant_keys = []
        path = self.vaccination_parameters_json_file
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
        #
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
        self.control_keys = control_keys
        return (state_keys, control_keys, parameter_keys, boundarycond_keys,
                constraint_keys, constant_keys)

    def read_bocop_parameters_values(self,
                                     parameters_prefix_file_name=
                                     'vaccination_parameters'):
        [state_keys, control_keys, parameter_keys,
         boundarycond_keys, constraint_keys, constant_keys] = \
            self.read_bocop_parameters_keys(parameters_prefix_file_name=
                                            'vaccination_parameters')
        path = self.bocop_solution_file
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

    def data_solution_save(self, file_name_prefix='vaccination_data_solution'):
        """
        Save scipy.integrate.solve_ivp output as pandas DataFrame
        Parameters
        ----------
        Returns
        -------
        df.plk: a pkl pandas DataFrame file.
        """
        time_now = datetime.now()
        dt_string = '_' + time_now.strftime("%b-%d-%Y_%H_%M")
        if not self.constant_control:
            file_name_prefix = 'no_' + file_name_prefix
            df = pd.DataFrame(self.x_sol_not_vaccine)
        else:
            df = pd.DataFrame(self.x_sol_vaccine)
            file_name_prefix = 'constant_' + file_name_prefix
        df.columns = self.time_state_solution_names
        file_name = file_name_prefix + dt_string + ".pkl"
        df.to_pickle(file_name)
        df.to_pickle(file_name_prefix + '.pkl')
        return

    def vaccination_constant_control(self, t):
        lambda_v = self.parameters['lambda_v']
        u_v_ = 0 * 1.0 * lambda_v
        self.u_v = u_v_
        return u_v_

    def get_lsoda_solution(self, constant_control=False):
        self.constant_control = constant_control

        def rhs_vaccination(t, y_, beta_s, beta_a, epsilon,
                            delta_e, delta_v, delta_r,
                            p,
                            alpha_a, alpha_s,
                            mu, theta,
                            lambda_v, constant_control_):
            """"
                Version designed in May 22-2020,
                according to the normalized version and
                lockdown process of the above Kermack Model.
                @param beta_a:
                    rate of asymptomatic contact
                @param beta_s:
                    rate of asymptomatic contact
                @param epsilon:
                    vaccination failure rate
                @param delta_e:
                    the inverse denote  latency time
                @param delta_v:
                    the inverse denote immunity time
                @param delta_r:
                    the inverse denote immunity time
                @param p:
                @param alpha_s:
                @param alpha_a:
                @type mu:
                @param theta:
                @param lambda_v:
                @param constant_control_:
                @param y_:
                @type t: object

            Parameters
            ----------
            constant_control
            """
            s, e, i_s, i_a, r, d, v, x_vaccine_counter = y_
            #
            if not constant_control:
                lambda_v = 0
            u_v = self.vaccination_constant_control(t)
            n_bar = s + e + i_s + i_a + r + v
            force_infection = (beta_s * i_s + beta_a * i_a) / n_bar
            rhs_s = mu * n_bar - force_infection * s - (
                    mu + (lambda_v + u_v)) * s \
                    + (delta_v ** (-1.0)) * v + delta_r * r
            rhs_e = force_infection * ((1.0 - epsilon) * v + s) - (
                    mu + delta_e) * e
            rhs_i_s = p * delta_e * e - (mu + alpha_s) * i_s
            rhs_i_a = (1 - p) * delta_e * e - (mu + alpha_a) * i_a
            rhs_r = theta * alpha_s * i_s + alpha_a * i_a - (mu + delta_r) * r
            rhs_d = (1 - theta) * alpha_s * i_s
            rhs_v = (lambda_v + u_v) * s - \
                    (1.0 - epsilon) * force_infection * v - \
                    (mu + (delta_v ** (-1.0))) * v
            rhs_x_vaccine_counter = (lambda_v + u_v) * (s + e + i_a + r)
            rhs = np.array([rhs_s, rhs_e,
                            rhs_i_s, rhs_i_a,
                            rhs_r, rhs_d,
                            constant_control_ * rhs_v,
                            rhs_x_vaccine_counter])
            return rhs

        prm = self.parameters
        T = prm["T"]
        n_whole = prm["n_pop"]
        s_zero = prm["S_0"]
        e_zero = prm["E_0"]
        i_s_zero = prm["I_S_0"]
        i_a_zero = prm["I_A_0"]
        r_zero = prm["R_0"]
        d_zero = prm["D_0"]
        v_zero = prm["V_0"]
        x_vaccination_zero = prm["X_0"]
        args = (prm["beta_s"], prm["beta_a"], prm["epsilon"],
                prm["delta_e"], prm["delta_v"], prm["delta_r"],
                prm["p"],
                prm["alpha_a"], prm["alpha_s"],
                prm["mu"], prm["theta"],
                prm["lambda_v"], self.constant_control)
        z_0 = np.array([s_zero,
                        e_zero,
                        i_s_zero,
                        i_a_zero,
                        r_zero,
                        d_zero,
                        v_zero,
                        x_vaccination_zero
                        ])
        sol_not_v = solve_ivp(rhs_vaccination,
                              [0.0, T],
                              z_0,
                              dense_output=False,
                              method='LSODA',
                              t_eval=self.t_eval,
                              events=None,
                              vectorized=False,
                              args=args
                              )
        y = np.c_[sol_not_v.t, sol_not_v.y.T]
        #
        if not self.constant_control:
            self.x_sol_not_vaccine = y
        else:
            self.x_sol_vaccine = y
        self.data_solution_save()
        return y
