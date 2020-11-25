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
from plotly.subplots import make_subplots
from scipy import integrate
from sage_scenarios_simulation import *


class CovidNumericalModelIncidence(CovidNumericalModel):

    def __init__(self,
                 bocop_solution_file,
                 bocop_parameters_json_file):
        super().__init__(bocop_solution_file,
                          bocop_parameters_json_file)

    def incidence_quadrature(self, df_solution):
        """

        Parameters
        ----------
        df_solution: Pandas data frame COVID19 model solution.

        Returns
        -------
        Quadrature of \int_{t_{i-1}}^{t_i} \delta_E E
        """
        prm = self.parameters
        bocop_run_prm = self.bocop_run_parameters
        p = prm["p"]
        delta_e = prm["delta_e"]
        t = df_solution['time']
        e = df_solution['e']
        f_incidence = p * delta_e * e
        initial_incidence = 0
        quad = integrate.cumtrapz(f_incidence, initial=initial_incidence)
        return quad

    def sage_plot_incidence(self, fig_file_name_prefix="incidence_fig_"):
        """

        Returns
        -------
        Deaths and reported infected incidence plots
        """
        prm = self.parameters
        self.load_pkl_data()
        optimal_controlled_sol_path = self.optimal_controlled_solution_path
        uncontrolled_solution_path = self.uncontrolled_solution_path
        constant_controlled_solution_file = \
            self.constant_controlled_solution_file
        df_oc = pd.read_pickle(optimal_controlled_sol_path)
        df_not_vaccination = pd.read_pickle(uncontrolled_solution_path)
        df_constant_vaccination = \
            pd.read_pickle(constant_controlled_solution_file)
        # optimal_i_s_incidence = self.incidence_computing(df_oc)
        # not_vaccination_incidence =
        # self.incidence_computing(df_not_vaccination)
        # constant_i_s_incidence = \
        #    self.incidence_computing(df_constant_vaccination)
        optimal_i_s_incidence = df_oc["y_inc"]
        not_vaccination_incidence = df_not_vaccination["y_inc"]
        constant_i_s_incidence = df_constant_vaccination["y_inc"]
        border_color_pallet = px.colors.sequential.ice
        fill_color_pallet = bokeh_palettes.all_palettes['Category20'][20]
        incidence_fig = make_subplots(
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
        incidence_fig.add_trace(
            go.Scatter(
                x=df_not_vaccination['time'],
                y=n_cdmx * df_not_vaccination["y_inc"],
                line=dict(color=border_color_pallet[1], width=.7, dash='dot'),
                legendgroup='Incidence_line',
                name='Without<br>vaccination',
                showlegend=True
            ),
            row=1, col=1
        )
        trace_optimal_incidence_line = go.Scatter(
            x=df_oc['time'],
            y=n_cdmx * df_oc["y_inc"],
            line=dict(color=border_color_pallet[1],
                      width=.7,
                      dash='solid'),
            legendgroup='Incidence_line',
            name='vaccination',
            showlegend=True
        )
        incidence_fig.add_trace(
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
        # ----------------------------------------------------------------------
        # constant_vaccination
        # ----------------------------------------------------------------------
        trace_constant_vac_i_s_incidence = go.Scatter(
            x=df_constant_vaccination["time"],
            y=n_cdmx * df_constant_vaccination["y_inc"],
            fill='tonexty',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[2]), .7)}",
            line=dict(color=fill_color_pallet[2], width=.7),
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
        trace_optimal_incidence_line_fill = go.Scatter(
            x=df_oc["time"],
            y=n_cdmx * df_oc["y_inc"],
            fill='tonexty',
            fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[3]), .7)}",
            line=dict(color=border_color_pallet[1],
                      width=.7,
                      dash='solid'),
            legendgroup='Mitigation',
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
        # TODO: Fix Legends
        incidence_fig.append_trace(trace_constant_vac_i_s_incidence, 1, 1)
        incidence_fig.append_trace(trace_optimal_incidence_line_fill, 1, 1)
        incidence_fig.append_trace(trace_optimal_incidence_line, 1, 1)
        incidence_fig.append_trace(trace_constant_vac_d, 1, 2)
        incidence_fig.append_trace(trace_optimal_d, 1, 2)
        #
        # Axis labels
        #
        incidence_fig.add_annotation(
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
        incidence_fig.update_xaxes(
            title_text="Time (days)",
            title_font=dict(size=10, family='Arial'),
            tickfont=dict(size=10, family='Arial'),
            row=1, col=1
        )
        incidence_fig.update_xaxes(
            title_text="Time (days)",
            title_font=dict(size=10, family='Arial'),
            tickfont=dict(size=10, family='Arial'),
            row=1, col=2
        )
        #
        incidence_fig.update_yaxes(
            tickfont=dict(size=10, family='Arial'),
            row=1, col=1
        )
        incidence_fig.update_yaxes(
            tickfont=dict(size=10, family='Arial'),
            row=1, col=2
        )
        incidence_fig.update_yaxes(
            tickfont=dict(size=10, family='Arial'),
            row=1, col=1
        )
        incidence_fig.update_yaxes(
            tickfont=dict(size=10, family='Arial'),
            row=1, col=2
        )
        #
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
        incidence_fig.add_annotation(
            dict(
                text='Mitigation',
                align='left',
                font=dict(family="Arial",
                          size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=1.21,
                y=0.8)
        )
        # ----------------------------------------------------------------------
        # Saved Beds
        # ----------------------------------------------------------------------
        incidence_fig.add_annotation(
            dict(
                text='Saved Lives',
                align='left',
                font=dict(family="Arial",
                          size=10),
                showarrow=False,
                xref='paper',
                yref='paper',
                x=1.23,
                y=0.6
            )
        )
#
        delta_r_label = round(prm["delta_r"], 2)
        delta_v_label = round(prm["delta_v"], 2)

        if delta_v_label == 0.0:
            delta_v_label = 'lifelong'
        if delta_r_label == 0.0:
            delta_r_label = 0.0 # 'lifelong'
        else:
            delta_r_label = delta_r_label ** (-1)
#
        data_label = {'eps': prm["epsilon"],
                      'delta_r': delta_r_label,
                      'delta_v': delta_v_label,
                      'time_unit': 'days'}
#
        if delta_v_label == 'lifelong':
            str_vaccination_par = \
                r'$\epsilon={:1.2f},' \
                r' \quad 1 / \delta_R={:1.1f} \mathsf{{{:>5}}}' \
                r' \quad 1 / \delta_V={:1.1f} \mathsf{{{:>9}}}$'.format(
                    data_label['eps'],
                    data_label['delta_r'],
                    data_label['time_unit'],
                    data_label['delta_v'])
        else:
            str_vaccination_par = \
                r'$\epsilon={:1.2f},' \
                r' \quad 1 / \delta_R={:1.1f} \ \mathsf{{{:>5}}}' \
                r' \quad 1 / \delta_V={:1.1f} \ \mathsf{{{:>5}}}$'.format(
                    data_label['eps'],
                    data_label['delta_r'],
                    data_label['time_unit'],
                    data_label['delta_v'],
                    data_label['time_unit'],
                )
        incidence_fig.add_annotation(
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
        incidence_fig.add_annotation(
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
        incidence_fig.update_layout(
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
        # incidence_figure alpha number
        # ----------------------------------------------------------------------
        incidence_fig.add_annotation(
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
        incidence_fig.add_annotation(
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

        for i in incidence_fig['layout']['annotations']:
            i['font'] = dict(family="Arial", size=10)
        if not os.path.exists("images"):
            os.mkdir("images")
        # incidence_fig.write_image("images/incidence_fig1.pdf")
        golden_width = 748  # 718  # width in px
        golden_ratio = 1.618
        golden_height = golden_width / golden_ratio
        time_now = datetime.now()
        dt_string = time_now.strftime("%b-%d-%Y_%H_")
        run_tag = self.run_tag
        path_fig_pdf = fig_file_name_prefix + \
                   '_' + \
                   dt_string + \
                   run_tag + '.pdf'
        path_fig_png = fig_file_name_prefix + \
                   '_' + \
                   dt_string + \
                   run_tag + '.png'
        # print(path_fig_png)
        pio.write_image(incidence_fig,
                        "./images/pdf/incidence/" + path_fig_pdf,
                        width=golden_width,
                        height=golden_height,
                        scale=10)
        pio.write_image(incidence_fig,
                        file="./images/png/incidence/" + path_fig_png,
                        width=golden_width,
                        height=golden_height,
                        scale=10)
        pio.write_image(incidence_fig,
                        file="./images/pdf/incidence/incidence_fig.pdf",
                        width=golden_width,
                        height=golden_height,
                        scale=10)
        pio.write_image(incidence_fig,
                        file="./images/png/incidence/incidence_fig.png",
                        width=golden_width,
                        height=golden_height,
                        scale=10)
        pio.write_image(incidence_fig,
                        file="./images/pdf/incidence/incidence_fig.pdf",
                        width=golden_width,
                        height=golden_height,
                        scale=10)
        # TODO: Fix  parameters_legend
        # print(path_fig)
        # TODO:  Edit legend respect to groups
        # incidence_fig.show()
