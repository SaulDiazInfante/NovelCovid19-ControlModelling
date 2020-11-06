import numpy as np
from scipy.integrate import solve_ivp
import json
import re
# import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime


class NumericsCovid19:
    def __init__(self,
                 bocop_solution_file,
                 uncontrolled_solution_path='./vaccination_pkl_solutions' +
                                            '/not_vaccination_solution.pkl',
                 optimal_controlled_solution='/vaccination_pkl_solutions' +
                                             '/bocop_solution.pkl',
                 bocop_parameters_json_file='./bocop_run_parameters/' +
                                            'bocop_run_parameters' +
                                            '.json',
                 vaccination_parameters_json_file='./vaccination_model_' +
                                                  'parameters'
                                                  + '/vaccination_parameters' +
                                                  '.json'):
        self.run_tag = 'None'
        self.constant_control = False
        self.bocop_solution_file = bocop_solution_file
        self.read_bocop_parameters_values()
        with open(bocop_parameters_json_file) as json_file:
            bocop_prm_ = json.load(json_file)
        with open(vaccination_parameters_json_file) as json_file:
            vaccination_prm_ = json.load(json_file)
        self.bocop_run_parameters = bocop_prm_
        self.parameters = vaccination_prm_
        self.time_state_solution_names = ['time', 's', 'e',
                                          'i_s', 'i_a',
                                          'r', 'd', 'v',
                                          'x_vac', 'y_inc', 'x_super_zero']

        # TODO: 'cost'
        n_steps = self.bocop_run_parameters['discretization.steps']
        self.state_dim = self.time_state_solution_names.__len__()
        self.t_eval = \
            np.linspace(0, self.parameters['T'],
                        num=n_steps,
                        endpoint=True)
        self.x_sol_vaccine = np.zeros([n_steps, self.state_dim])
        self.x_sol_not_vaccine = np.zeros([n_steps, self.state_dim])
        self.df_sol_cte_vaccine = pd.DataFrame(self.x_sol_vaccine)
        self.df_sol_not_vaccine = pd.DataFrame(self.x_sol_not_vaccine)
        self.u_v = np.zeros([self.state_dim, 1])
        self.vaccination_parameters_json_file = vaccination_parameters_json_file
        self.uncontrolled_solution_path = uncontrolled_solution_path
        self.optimal_controlled_solution = optimal_controlled_solution
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
        # time_now = datetime.now()
        # dt_string = time_now.strftime("%b-%d-%Y_%H_%M")
        # path_json = parameters_prefix_file_name + dt_string + '.json'
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
                # print(str_headers_names.string, i)
            if str_footer_names:
                # print(str_footer_names.string, i)
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
                #            # print('\t', str_state, '\t\tat line:', j)
                state_keys.append(str_state)
            if x_str_control:
                str_control = x_str[3]
                #            # print('\t', str_control, '\t\tat line:', j)
                control_keys.append(str_control)
            if x_str_parameter:
                str_parameter = x_str[3]
                #            # print('\t', str_parameter, '\t\tat line:', j)
                parameter_keys.append(str_parameter)
            if x_str_boundarycond:
                str_boundarycond = x_str[3]
                #            # print('\t', str_boundarycond, '\t\tat line:', j)
                boundarycond_keys.append(str_boundarycond)
            if x_str_constraint:
                str_constraint = x_str[3]
                #            # print('\t', str_constraint, '\t\tat line:', j)
                constraint_keys.append(str_constraint)
            if x_str_constant:
                str_constant = x_str[3]
                #            # print('\t', str_constant, '\t\tat line:', j)
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
                                     './vaccination_model_parameters/' +
                                     'vaccination_parameters'):
        [state_keys, control_keys, parameter_keys,
         boundarycond_keys, constraint_keys, constant_keys] = \
            self.read_bocop_parameters_keys(parameters_prefix_file_name=
                                            './vaccination_model_parameters/' +
                                            'vaccination_parameters'
                                            )
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
                # print(i, str_header_bounds.string)
                pointer = i + 11
                for j in np.arange(pointer, pointer + boundary_cond_dim):
                    x_val_sci = pattern_sci_float.search(all_lines[j])
                    x_val = pattern_value.search(all_lines[j])
                    if x_val_sci:
                        # print(j + 1, x_val_sci.string)
                        value = x_val_sci.group()
                    elif x_val:
                        # print(j + 1, x_val.string)
                        value = x_val.group()
                    boundarycond_values.append(np.float(value))
            #
            if str_header_constants:
                # print(str_header_constants.string, i)
                pointer = i + 8
                for j in np.arange(pointer, pointer + constant_dim):
                    x_val_sci = pattern_sci_float.search(all_lines[j])
                    x_val = pattern_value.search(all_lines[j])
                    if x_val_sci:
                        # print(j + 1, x_val_sci.string)
                        value = x_val_sci.group()
                    elif x_val:
                        # print(j + 1, x_val.string)
                        value = x_val.group()
                    constants_values.append(np.float(value))
            #
            if str_header_optimal_time:
                # print(str_header_optimal_time.string, i)
                pointer = i + 1
                x_val_sci = pattern_sci_float.search(all_lines[pointer])
                x_val = pattern_value.search(all_lines[pointer])
                if x_val_sci:
                    # print(pointer, x_val_sci.string)
                    value = x_val_sci.group()
                    # print(value)
                elif x_val:
                    # print(pointer, x_val.string)
                    value = x_val.group()
                constants_values.append(np.float(value))
                constant_keys.append('T')
            i = i + 1
        data_constants = dict(zip(constant_keys, constants_values))
        data_boundarycond = dict(zip(boundarycond_keys, boundarycond_values))
        data = data_constants
        data.update(data_boundarycond)
        data_json = json.dumps(data)
        # print(json.dumps(data, indent=4))
        with open(path_json, 'w') as json_file:
            json.dump(data, json_file, indent=4)
        with open(parameters_prefix_file_name + '.json', 'w') as json_file:
            json.dump(data, json_file, indent=4)
        df_record_parameters = pd.DataFrame(data, index=[1])
        return df_record_parameters

    def get_run_bocop_parameters(self,
                                 prefix_file_name='./bocop_run_parameters/' +
                                                  'bocop_run_parameters'
                                 ):
        bocop_file_name = self.bocop_solution_file

        def check_token(x_str_, x_val_, i_):
            if x_str_ and x_val_:
                string_prop_name_ = x_str.group()
                str_prop_val_ = np.int(x_val[0])
                names.append(string_prop_name_)
                values.append(str_prop_val_)
                # print(i_, string_prop_name_, str_prop_val_)

        def check_token_str(x_str_, x_val_, i_):
            if x_str_ and x_val_:
                string_prop_name_ = x_str.group()
                str_prop_val_ = x_val_[0]
                names.append(string_prop_name_)
                values.append(str_prop_val_)
                # print(i_, string_prop_name_, str_prop_val_)

        path = bocop_file_name
        names = []
        values = []
        with open(path) as f:
            all_lines = f.readlines()
        i = 1
        # ----------------------------------------------------------------------
        # Regular expressions for pattern search
        # ----------------------------------------------------------------------
        # time regular expressions
        pattern_str = re.compile(r'time\.([a-z]*)')
        # pattern_str_none = re.compile(r'none')
        pattern_str_free = re.compile(r'final')
        pattern_int_value = re.compile(r'\d+')
        # Dimension re
        pattern_str_state_dimension = re.compile(r'state.dimension')
        pattern_str_control = re.compile(r'control.dimension')
        pattern_str_algebraic = re.compile(r'algebraic.dimension')
        pattern_str_parameter = re.compile(r'parameter.dimension')
        pattern_str_constant = re.compile(r'constant.dimension')
        pattern_str_boundarycond = re.compile(r'boundarycond.dimension')
        pattern_str_constraint = re.compile(r'constraint.dimension')
        pattern_int_value = re.compile(r'\d+')
        # Discretization re:
        pattern_str_steps = re.compile(r'discretization.steps')
        pattern_str_method = re.compile(r'discretization.method')
        pattern_value_str_method = \
            re.compile(r'discretization.method string\s+(\S+)')
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
        # solution file name
        pattern_str_solution_file = re.compile(r'solution.file')
        pattern_value_str_solution_file = \
            re.compile(r'solution.file string\s+(\S+)')
        # Optimized parameters: final_time
        pattern_str_optimized_time = re.compile(r'Parameters')
        # Search according to above patterns
        pattern_sci_float = re.compile(
            '[-+]?[\d]+\.?[\d]*[Ee](?:[-+]?[\d]+)?')
        #
        pattern_value = re.compile(r'\d+\.\d+|\d+')
        i = 0
        j = 0
        final_time = 1.0
        value = 0.0
        for line in all_lines:
            x_str = pattern_str.search(line)
            x_val = pattern_int_value.findall(line)
            x_str_final = pattern_str_free.findall(line)
            # x_str_none = pattern_str_free.findall(line)
            if x_str and x_val:
                string_prop_name = 'time.' + x_str.group(1)
                str_prop_val = np.int(x_val[0])
                names.append(string_prop_name)
                values.append(str_prop_val)
                # print(i, string_prop_name, str_prop_val)
            elif x_str and x_str_final:
                string_prop_name = 'time.' + 'free'
                str_prop_val = x_str_final[0]
                names.append(string_prop_name)
                values.append(str_prop_val)
                # print(i, string_prop_name, str_prop_val)
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
            x_str = pattern_str_solution_file.search(line)
            x_val = pattern_value_str_solution_file.findall(line)
            check_token_str(x_str, x_val, i)
            x_str_optimized_time_header = pattern_str_optimized_time.search(
                line)
            if x_str_optimized_time_header:
                pointer = i
                final_time = all_lines[pointer + 1]
                x_val_sci = pattern_sci_float.search(final_time)
                x_val = pattern_value.search(final_time)
                if x_val_sci:
                    # print(j + 1, x_val_sci.string)
                    value = x_val_sci.group()
                elif x_val:
                    # print(j + 1, x_val.string)
                    value = x_val.group()
                final_time = np.float(value)
                names.append('final_time')
                values.append(final_time)
            i = i + 1
        time_now = datetime.now()
        dt_string = time_now.strftime("%b-%d-%Y_%H_%M")
        path_json = prefix_file_name + '_' + dt_string + '.json'
        #
        bocop_data = dict(zip(names, values))
        bocop_data_df = pd.DataFrame(bocop_data, index=[1])
        bocop_data_df.to_json(path_json)
        # print(json.dumps(bocop_data, indent=4))
        with open(path_json, 'w') as json_file:
            json.dump(bocop_data, json_file, indent=4)
        with open(prefix_file_name + '.json', 'w') as json_file:
            json.dump(bocop_data, json_file, indent=4)
        with open(prefix_file_name + '.json') as json_file:
            bocop_run_prm_ = json.load(json_file)
            self.bocop_run_parameters = bocop_run_prm_
        return bocop_data_df

    def simulation_run_tagger(self):
        prm = self.parameters
        cov_tag = str(prm['coverage'])
        time_horizon_tag = str(prm['T'])
        delta_v_tag = str(prm['delta_v'])
        delta_r_tag = str(prm['delta_r'])
        epsilon_tag = str(prm['epsilon'])
        run_tag = 'run' + \
                  '_coverage_' + cov_tag + \
                  '_tHor_' + time_horizon_tag + \
                  '_delta_r_' + delta_r_tag + \
                  '_delta_v_' + delta_v_tag + \
                  '_epsilon_' + epsilon_tag
        self.run_tag = run_tag
        return run_tag

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
        pkl_data_folder = './vaccination_pkl_solutions/pkl/'
        csv_data_folder = './vaccination_pkl_solutions/csv/'
        run_tag = self.simulation_run_tagger()
        file_name = pkl_data_folder + file_name_prefix + \
                    run_tag + dt_string + ".pkl"
        file_name_csv = csv_data_folder + file_name_prefix + \
                        run_tag + dt_string + ".csv"
        df.to_pickle(file_name)
        df.to_csv(file_name_csv)
        df.to_pickle(pkl_data_folder + file_name_prefix + '.pkl')
        df.to_csv(csv_data_folder + file_name_prefix + '.csv')
        return

    def vaccination_constant_control(self, t):
        lambda_v = self.parameters['lambda_v']
        u_v_ = 1.0 * lambda_v
        self.u_v = u_v_
        return u_v_

    def get_lsoda_solution(self, constant_control=False):
        self.constant_control = constant_control
        prm = self.parameters
        a_d = prm['a_D']
        a_e = prm['a_E']

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
            s, e, i_s, i_a, r, d, v, x_vaccine_counter, y_inc, x_zero = y_
            #
            if not constant_control:
                lambda_v = 0
            n_bar = s + e + i_s + i_a + r + v
            force_infection = (beta_s * i_s + beta_a * i_a) / n_bar
            if delta_v == 0:
                delta_v_indicator = 0.0
            else:
                delta_v_indicator = 1.0 / delta_v
            rhs_s = mu * n_bar - force_infection * s - (
                    mu + lambda_v) * s \
                    + delta_v_indicator * v + delta_r * r
            rhs_e = force_infection * ((1.0 - epsilon) * v + s) - (
                    mu + delta_e) * e
            rhs_i_s = p * delta_e * e - (mu + alpha_s) * i_s
            rhs_i_a = (1 - p) * delta_e * e - (mu + alpha_a) * i_a
            rhs_r = (1 - theta) * alpha_s * i_s + alpha_a * i_a \
                    - (mu + delta_r) * r
            rhs_d = theta * alpha_s * i_s
            rhs_v = lambda_v * s - \
                    (1.0 - epsilon) * force_infection * v - \
                    (mu + delta_v_indicator) * v
            rhs_x_vaccine_counter = lambda_v * (s + e + i_a + r)
            rhs_y_inc = p * delta_e * e
            rhs_x_zero = 0.0

            rhs = np.array([rhs_s, rhs_e,
                            rhs_i_s, rhs_i_a,
                            rhs_r, rhs_d,
                            constant_control_ * rhs_v,
                            rhs_x_vaccine_counter, rhs_y_inc, rhs_x_zero])
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
        y_inc_0 = prm["Y_inc_0"]
        x_zero_0 = prm["xZero_0"]
#
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
                        x_vaccination_zero,
                        y_inc_0,
                        x_zero_0
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
