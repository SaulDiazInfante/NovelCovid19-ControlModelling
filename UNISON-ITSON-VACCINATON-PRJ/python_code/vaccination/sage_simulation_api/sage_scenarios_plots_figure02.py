import numpy as np
import json
import pandas as pd
import os
import plotly.express as px
import plotly.graph_objects as go
import decimal
import bokeh.palettes as bokeh_palettes
import plotly.io as pio
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
    r_v_0_ = fac * r_00
    opt_r_v_0 = opt_fac * r_00
    return r_00, r_v_0_, opt_r_v_0


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
########################################################################
# Solution without vaccination
########################################################################
fig02 = make_subplots(
    rows=3, cols=2,
    specs=[
        [{"rowspan": 2}, {"rowspan": 2}],
        [None, None],
        [{"colspan": 2}, None]
    ],
    subplot_titles=(
        "Coverage",
        "Cost",
        "Vaccination Policy"
    ),
    horizontal_spacing=0.17,
    vertical_spacing=0.3,
)
n_cdmx = 100000.0 / prm['n_pop']
########################################################################
# constant_vaccination
########################################################################
coverage_t = constant_vaccination_coverage(df_constant_vaccination)
trace_constant_vac_coverage = go.Scatter(
    x=df_constant_vaccination['time'],
    y=coverage_t,
    fill='tozeroy',
    fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[18]), 0.7)}",
    line=dict(color=fill_color_pallet[18], width=0.7),
    legendgroup='Coverage',
    name='(CP) Coverage',
    showlegend=True
)
cost_t = constant_vaccination_cost(df_constant_vaccination)
trace_constant_vac_cost = go.Scatter(
    x=df_constant_vaccination['time'],
    y=n_cdmx * cost_t,
    fill='tonexty',
    fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[16]), 0.7)}",
    line=dict(color=fill_color_pallet[16], width=1.0),
    legendgroup='Cost',
    name='(CP) Cost',
    showlegend=True
)
#######################################################################
# Optimal solution trace
#######################################################################
#
trace_optimal_vac_coverage = go.Scatter(
    x=df_oc['time'],
    y=df_oc['x_vac'],
    line=dict(color=fill_color_pallet[19], width=.7),
    fill='tozeroy',
    fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[19]), 0.5)}",
    legendgroup='Coverage',
    name="(OP) Coverage",
    showlegend=True
)
#
lambda_base_v_t = lambda_v_base * np.ones(df_oc["u_V"].shape[0])
trace_optimal_vaccination_policy = go.Scatter(
    x=df_oc['time'],
    y=n_cdmx * (df_oc['u_V'] + lambda_base_v_t),
    fill='tozeroy',
    fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[19]), 0.5)}",
    line=dict(color=fill_color_pallet[0], width=0.7),
    name="Optimal<br>vaccination<br>policy",
    legendgroup='Policy',
    showlegend=True
)
trace_optimal_cost = go.Scatter(
    x=df_oc['time'],
    y=n_cdmx * df_oc['x_zero'],
    line=dict(color=fill_color_pallet[16], width=1.0),
    fillcolor=f"rgba{(*hex_to_rgb(fill_color_pallet[17]), 0.4)}",
    fill='tozeroy',
    legendgroup='Cost',
    name='(OP) Cost',
    showlegend=True
)
#
trace_vaccination_base = go.Scatter(
    x=df_oc['time'],
    y=n_cdmx * lambda_v_base * np.ones(df_oc["u_V"].shape[0]),
    line=dict(color=fill_color_pallet[0], width=0.7, dash='dot'),
    name='Constant<br>policy',
    showlegend=True
)
fig02.append_trace(trace_constant_vac_coverage, 1, 1)
fig02.append_trace(trace_optimal_vac_coverage, 1, 1)
fig02.append_trace(trace_optimal_cost, 1, 2)
fig02.append_trace(trace_constant_vac_cost, 1, 2)
fig02.append_trace(trace_vaccination_base, 3, 1)
fig02.append_trace(trace_optimal_vaccination_policy, 3, 1)
#
# Axis labels
#
fig02.add_annotation(
    dict(
        text='Population proportion',
        align='left',
        textangle=-90,
        font=dict(family="Arial", size=10),
        showarrow=False,
        xref='paper',
        yref='paper',
        x=-0.075,
        y=.9)
)
fig02.update_xaxes(
    title_text="Time (days)",
    title_font=dict(size=10, family='Arial'),
    tickfont=dict(size=10, family='Arial'),
    row=1, col=1
)
fig02.update_xaxes(
    title_text="Time (days)",
    title_font=dict(size=10, family='Arial'),
    tickfont=dict(size=10, family='Arial'),
    row=1, col=2
)
fig02.update_xaxes(
    title_text="Time (days)",
    title_font=dict(size=10, family='Arial'),
    tickfont=dict(size=10, family='Arial'),
    row=3, col=1
)
#
fig02.update_yaxes(
    tickfont=dict(size=10, family='Arial'),
    row=1, col=1
)
fig02.update_yaxes(
    tickfont=dict(size=10, family='Arial'),
    row=1, col=2
)
#
fig02.update_yaxes(
    tickfont=dict(size=10, family='Arial'),
    row=1, col=3)
fig02.update_yaxes(
    tickfont=dict(size=10, family='Arial'),
    row=2, col=1
)

fig02.update_yaxes(
    tickfont=dict(size=10, family='Arial'),
    row=2, col=2
)
fig02.add_annotation(
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
fig02.update_yaxes(
    tickfont=dict(size=10, family='Arial'),
    row=3, col=1
)
fig02.add_annotation(
    dict(
        text="Cost per 100,000 inhabitants",
        align='left',
        textangle=-90,
        font=dict(family="Arial", size=10),
        showarrow=False,
        xref='paper',
        yref='paper',
        x=0.50,
        y=0.95)
)
#
fig02.update_yaxes(tickfont=dict(size=10, family='Arial'),
                   row=3, col=3)
#
r_0, r_v_0, r_opt_v_0 = reproductive_number(**prm)
c_r_0 = decimal.Decimal(r_0)
c_r_v_0 = decimal.Decimal(r_v_0)
c_r_opt_v_0 = decimal.Decimal(r_opt_v_0)

data_label = {'eps': prm["epsilon"],
              'delta_v': round(prm["delta_v"], 1),
              'time_unit': 'days'}
str_vaccination_par = \
    r'$\epsilon={:1.2f}, \quad \delta_V ^{{-1}}={:1.1f} \ \mathtt{{{' \
    r':>5}}}$'.format(
        data_label['eps'],
        data_label['delta_v'],
        data_label['time_unit'])
str_policy_legend = '<b>Vaccination Policy:</b><br>' + \
                    '    Optimal (OP)<br>' + \
                    '    Constant (CP)'
fig02.add_annotation(
    dict(
        text=str_policy_legend,
        align='left',
        font=dict(family="Arial",
                  size=12),
        showarrow=False,
        xref='paper',
        yref=''
             'paper',
        x=1.325,
        y=0.05
    )
)

fig02.add_annotation(
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
fig02.update_layout(
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
# fig02ure alpha number
# ----------------------------------------------------------------------
fig02.add_annotation(
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
fig02.add_annotation(
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
fig02.add_annotation(
    dict(
        text='<b>C)',
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
#
for i in fig02['layout']['annotations']:
    i['font'] = dict(family="Arial", size=10)
if not os.path.exists("images"):
    os.mkdir("images")
# fig02.write_image("images/fig021.pdf")
golden_width = 718  # width in px
golden_ratio = 1.618
pio.kaleido.scope.default_format = "eps"
pio.kaleido.scope.default_width = golden_width
pio.kaleido.scope.default_height = golden_width / golden_ratio
pio.kaleido.scope.default_scale = 0.5
fig02.to_image(format="pdf", engine="kaleido")
fig02.write_image("images/fig02.pdf")
