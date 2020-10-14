# Simple SEIR model
# by: Travis N. Vaught
#
# This is a bokeh version of the model shown here:
#   https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/
#   and here:
#   https://towardsdatascience.com/social-distancing-to-slow-the-coronavirus-768292f04296
#
# Dependencies:
#   numpy, scipy, bokeh
#
# Run by typing ```bokeh serve sir.py```
#   then point your browser to: http://localhost:5006/sir
#
#

# Major Package Imports
import numpy as np
from scipy.integrate import odeint

from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import Slider, Div
from bokeh.io import curdoc
from bokeh.layouts import column, row


# Total population, N.
N = 1000
# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 1, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
E0 = I0 - R0

# Incubation period, 1/alpha, contact rate, beta, and mean recovery rate, 1/gamma.
alpha_inv = Slider(title="Incubation (latency) Period, 1/alpha (in days)",
                   value=3, start=1, end=10, step=1)
beta = Slider(title="Contact rate, beta",
              value=0.35, start=.01, end=0.80, step=0.01)
gamma_inv = Slider(title="mean recovery rate, 1/gamma, (in days)",
                   value=10, start=4, end=25, step=1)

# A grid of time points (in days)
t = np.linspace(0, 200, 200)

# The SIR model differential equations.
def deriv(y, t, N, alpha, beta, gamma):
    S, E, I, R = y
    dSdt = -beta * S * I / N
    dEdt = beta * S * I / N - alpha * E
    dIdt = alpha * E - gamma * I
    dRdt = gamma * I
    return dSdt, dEdt, dIdt, dRdt

# Initial plot
# Slider Values
a = 1/alpha_inv.value
b = beta.value
g = 1/gamma_inv.value

# Initial conditions vector
y0 = S0, E0, I0, R0

# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, a, b, g))
S, E, I, R = ret.T
Rsub0 = b/g
title_text = "SEIR Model R₀ = {0:.2f}".format(Rsub0)
source = ColumnDataSource(dict(t=t, S=S/N, E=E/N, I=I/N, R=R/N))


def update_data(attrname, old, new):
    # Get the current slider values
    a = 1/alpha_inv.value
    b = beta.value
    g = 1/gamma_inv.value
    # Initial conditions vector
    y0 = S0, E0, I0, R0
    # Integrate the SIR equations over the time grid, t.
    ret = odeint(deriv, y0, t, args=(N, a, b, g))
    S, E, I, R = ret.T
    Rsub0 = b/g
    title_text = "SEIR Model R₀ = {0:.2f}".format(Rsub0)
    source.data = dict(t=t, S=S/N, E=E/N, I=I/N, R=R/N)
    fig.title.text = title_text


for w in [alpha_inv, beta, gamma_inv]:
    w.on_change('value', update_data)


# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = figure(width=600, height=400, title=title_text)
fig.line('t', 'S', line_color='blue', alpha=0.5,
         legend_label='Susceptible', source=source)
fig.line('t', 'E', line_color='orange',
         legend_label='Exposed', source=source)
fig.line('t', 'I', line_color='red', alpha=0.5,
         legend_label='Infected', source=source)
fig.line('t', 'R', line_color='darkgreen', alpha=0.5,
         legend_label='Recovered/Removed', source=source)
fig.legend.location = 'center_right'
fig.xaxis.axis_label = 'time (days)'
fig.yaxis.axis_label = '% of population'

# Set up layouts and add to document
input_label ="<h3>SEIR Model Input Values </h3> "
inputs = column(Div(text="<h3>SEIR Model Input Values</h3>"),
                alpha_inv, beta, gamma_inv)

curdoc().add_root(row(inputs, fig, width=800))
curdoc().title = "Simple SIR Model"