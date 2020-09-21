import numpy as np
import re
import json
import pandas as pd
import os
import plotly.express as px
import plotly.graph_objects as go
import decimal
import bokeh.palettes as bokeh_palettes
import plotly.io as pio
from datetime import datetime
from plotly.subplots import make_subplots
from scipy import integrate


def reproductive_number(**kwargs):
    beta_s = kwargs['beta_s']
    beta_a = kwargs['beta_a']
    delta_e = kwargs['delta_e']
    alpha_s = kwargs['alpha_s']
    alpha_a = kwargs['alpha_a']
    lambda_v = kwargs['lambda_v']
    epsilon = kwargs['epsilon']
    delta_v = kwargs['delta_v']
    mu = kwargs['mu']
    p = kwargs['p']

    f_1 = (p * delta_e * beta_s) / ((alpha_s + mu) * (mu + delta_e))
    f_2 = ((1.0 - p) * delta_e * beta_a) / ((mu + alpha_a) * (mu + delta_e))
    r_00 = f_1 + f_2

    fac = (mu + delta_v + (1.0 - epsilon) * lambda_v) / \
          (mu + delta_v + lambda_v)
    u_v_0 = df_oc['u_V'][0]
    opt_fac = (mu + delta_v + (1.0 - epsilon) * (lambda_v + u_v_0)) / \
              (mu + delta_v + (lambda_v + u_v_0))

    r_v_0 = fac * r_00
    opt_r_v_0 = opt_fac * r_00
    return r_00, r_v_0, opt_r_v_0


def hex_to_rgb(hex_color: str) -> tuple:
    hex_color = hex_color.lstrip("#")
    if len(hex_color) == 3:
        hex_color = hex_color * 2
    return int(hex_color[0:2], 16), \
           int(hex_color[2:4], 16), int(hex_color[4:6], 16)


def constant_vaccination_cost(df_solution):
    """

    """
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


def constant_vaccination_coverage(df_solution):
    """

    """
    lambda_v = prm['lambda_v']
    t = df_solution['time']
    s = df_solution['s']
    e = df_solution['e']
    i_a = df_solution['i_a']
    r = df_solution['r']
    f = lambda_v * (s + e + i_a + r)
    coverage_t_ = integrate.cumtrapz(f, t, initial=0)
    return coverage_t_


def load_parameters(file_name='/vaccination_model_parameters/' +
                    'vaccination_parameters.json '):
    """ Load the parameters given a JSON file.
    Parameters
    ----------
    file_name: filename with parameters values in JSON format.

    Returns
    -------
    prm: dictionary
    """
    with open(file_name) as json_file:
        prm_ = json.load(json_file)
    return prm_


prm = load_parameters('../vaccination_model_parameters/' +
                      'vaccination_parameters.json')
lambda_v_base = prm["lambda_v"]
bocop_sol_path = '../vaccination_pkl_solutions/bocop_solution.pkl'
constant_vaccination_path = '../vaccination_pkl_solutions/' + \
                            'constant_vaccination_data_solution.pkl'
no_vaccination_sol_path = '../vaccination_pkl_solutions/' + \
                          'no_vaccination_data_solution.pkl'
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
n_cdmx = prm['n_pop']
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
            df_oc["u_V"].max() \
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
kwargs = {"beta_s": prm["beta_s"],
          "beta_a": prm["beta_a"],
          "delta_e": prm["delta_e"],
          "p": prm["p"],
          "alpha_a": prm["alpha_a"],
          "alpha_s": prm["alpha_s"],
          "mu": prm["mu"],
          "theta": prm["theta"],
          "lambda_v": prm["lambda_v"],
          "delta_v": prm["delta_v"],
          "epsilon": prm["epsilon"]
          }

r_0, r_v_0, r_opt_v_0 = reproductive_number(**kwargs)
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
