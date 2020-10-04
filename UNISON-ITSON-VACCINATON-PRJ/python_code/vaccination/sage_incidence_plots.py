import decimal
import json
import os
import re
from datetime import datetime
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

    def __init__(self):
        super().__init__()
        self.cumulatively_reported_incidence = np.zeros(300)
        self.reported_incidence = np.zeros(300)
        self.sampler_time = np.zeros(300)

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

    def incidence_computing(self, df_solution):
        """

        Returns
        -------
            An numpy array with the disease incidence as time
            function.
            $$
                \Pi_I(t_i):= \int_{t \in [t_{i-1}, t_{i}]} p \delta_E E
            $$
        """
        prm = self.parameters
        t = df_solution['time']
        i_s = df_solution['i_s']
        operation_time_step_multiple = np.argmin(np.abs(1.0 - t))
        sampler_time = t[0:-1:operation_time_step_multiple]
        self.sampler_time = sampler_time
        incidence = np.zeros(sampler_time.shape[0])
        incidence_increments = np.zeros(sampler_time.shape[0])
        cumulative_incidence_t = self.incidence_quadrature(df_solution)
        cumulative_incidence_per_day = \
            cumulative_incidence_t[0:-1:operation_time_step_multiple]
        incidence_increments[1:] = np.ediff1d(cumulative_incidence_per_day)
        self.cumulatively_reported_incidence = cumulative_incidence_per_day
        self.reported_incidence = incidence_increments
        return incidence

    def sage_plot_incidence(self):
        """

        Returns
        -------
        Deaths and reported infected incidence plots
        """
        prm = self.parameters
        optimal_controlled_sol_path = self.optimal_controlled_solution_path
        uncontrolled_solution_path = self.uncontrolled_solution_path
        constant_controlled_solution_file = \
            self.constant_controlled_solution_file
        df_oc = pd.read_pickle(optimal_controlled_sol_path)
        df_not_vaccination = pd.read_pickle(uncontrolled_solution_path)
        df_constant_vaccination = pd.read_pickle(
            constant_controlled_solution_file)
        optimal_controlled_sol_path = self.optimal_controlled_solution_path
        uncontrolled_solution_path = self.uncontrolled_solution_path
        constant_controlled_solution_file = \
            self.constant_controlled_solution_file
        df_oc = pd.read_pickle(optimal_controlled_sol_path)
        df_not_vaccination = pd.read_pickle(uncontrolled_solution_path)
        df_constant_vaccination = pd.read_pickle(
            constant_controlled_solution_file)
        optimal_i_s_incidence = self.incidence_computing(df_oc)
        not_vaccination_incidence = self.incidence_computing(df_not_vaccination)
        constant_i_s_incidence = \
            self.incidence_computing(df_constant_vaccination)
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
        n_cdmx = prm["n_pop"] / 100000
        incidence_fig.add_trace(
            go.Scatter(
                x=self.sampler_time,
                y=n_cdmx * not_vaccination_incidence,
                line=dict(color=border_color_pallet[1], width=.7, dash='dot'),
                legendgroup='Incidence_line',
                name='Without<br>vaccination',
                showlegend=True
            ),
            row=1, col=1
        )
        trace_optimal_incidence_line = go.Scatter(
            x=self.sampler_time,
            y=n_cdmx * optimal_i_s_incidence,
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
            x=self.sampler_time,
            y=n_cdmx * constant_i_s_incidence,
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
            x=self.sampler_time,
            y=n_cdmx * optimal_i_s_incidence,
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
                y=0.78)
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
                y=0.57
            )
        )
#
        data_label = {'eps': prm["epsilon"],
                      'delta_v': round(prm["delta_v"], 1),
                      'time_unit': 'days'}
        str_vaccination_par = \
            r'$\epsilon={:1.2f}, \quad \delta_V ^{{-1}}={:1.1f} \ \mathtt{{{' \
            r':>5}}}$'.format(
                data_label['eps'],
                data_label['delta_v'],
                data_label['time_unit'])
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
                yref=''
                     'paper',
                x=1.295,
                y=0.05
            )
        )
        incidence_fig.update_layout(
            font=dict(family="Arial", size=10),
            template="simple_white",
            showlegend=True,
            title_text="R0: " + str(round(c_r_0, 6))
                       + '\t\tRv: ' + str(round(c_r_v_0, 6))
                       + '\t\tOC-Rv: '
                       + str(round(c_r_opt_v_0, 6)),
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
        golden_width = 718  # width in px
        golden_ratio = 1.618
        pio.kaleido.scope.default_format = "eps"
        pio.kaleido.scope.default_width = golden_width
        pio.kaleido.scope.default_height = golden_width / golden_ratio
        pio.kaleido.scope.default_scale = .50
        incidence_fig.to_image(format="pdf", engine="kaleido")
        incidence_fig.write_image("images/incidence_fig.pdf")
        # TODO:  Edit legend respect to groups
        # incidence_fig.show()
