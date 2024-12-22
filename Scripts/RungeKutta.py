import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import AutoMinorLocator
from collections import defaultdict
import os

name = 'RK_test'
aspect_ratio = 'square'
filenames = [r"D:\Univerzita\Numerické metody\Radial_problem\Task 1\RK_test.txt"]
x_column = 1
y_columns = [2,3,4]
show_legend = True
graph_title = r'Runge-Kutta test for linear harmonic oscillator'
legend_title = 0
legend_frame = False
white = False
positive_filter = False

# Size settings
SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 24
LEGEND_SIZE = 12
SIZE = 6

# Legends settings
X_axis = r'Number of points'
Y_axis = r'Absolute difference from exact solution'
legends = [[r'RK1',r'RK2',r'RK4']]
labels = [[0,0,0]]
labels_position = [[(0,0)]]
colors = [['blue','red','green']]
markers = [[]]
line_styles = [['--','-.','-']]
line_widths = [['','','']]
ylog = True
xlog = True

# Step settings
limit_x = False
limit_y = False
x_range = [1.24, 1.7]
y_range = [10 ** (-4), 10 ** (4)]
minstep_x = False
minstep_y = False
step_x = 0.1
step_y = 10

# LaTeX Settings
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Latin Modern Roman']

plt.rcParams['text.usetex'] = True
plt.rcParams["font.family"] = ["Latin Modern Roman"]
plt.rcParams['axes.titlepad'] = 10
plt.rcParams['axes.labelpad'] = 10
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.edgecolor'] = "#000000"
plt.rcParams["figure.autolayout"] = True
plt.rcParams["legend.handlelength"] = 3
plt.rcParams["legend.framealpha"] = 1
plt.rcParams["legend.borderpad"] = 0.8

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=LEGEND_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Fig. ratio settings
if aspect_ratio == 'golden_ratio':
    golden_ratio = (1 + 5 ** 0.5) / 2
    fig, ax = plt.subplots(figsize=(SIZE * golden_ratio, SIZE))
else:
    fig, ax = plt.subplots(figsize=(SIZE, SIZE))


for file_idx, filename in enumerate(filenames):
    data = np.loadtxt(filename, comments='#', delimiter=None)
    y_columns_indices = [col - 1 for col in y_columns]

    x = data[:, x_column-1]

    for col_idx, y_col in enumerate(y_columns):
        y = data[:, y_col-1]
        label = legends[file_idx][col_idx]
        if positive_filter:
            valid_indices = y > 0
            y = y[valid_indices]
            x = x[valid_indices]
        if label:
            if line_widths[file_idx][col_idx] == '':
                line, = ax.plot(x, y, label=label, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx])
            else:
                line, = ax.plot(x, y, label=label, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx], linewidth=line_widths[file_idx][col_idx])
        else:
            if line_widths[file_idx][col_idx] == '':
                line, = ax.plot(x, y, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx])
            else:
                line, = ax.plot(x, y, color=colors[file_idx][col_idx], linestyle=line_styles[file_idx][col_idx], linewidth=line_widths[file_idx][col_idx])
        if labels[file_idx][col_idx] != 0:
            abs_x_pos = labels_position[file_idx][col_idx][0]
            abs_y_pos = labels_position[file_idx][col_idx][1]
            ax.text(abs_x_pos, abs_y_pos, labels[file_idx][col_idx], color=colors[file_idx][col_idx], verticalalignment='bottom')



# Graph settings
ax.set_xlabel(X_axis)
ax.set_ylabel(Y_axis)
if limit_x:
    ax.set_xlim(x_range)
if limit_y:
    ax.set_ylim(y_range)
ax.set_title('')
if show_legend:
    if legend_title != 0:
        if white:
            ax.legend(title=legend_title, facecolor='white', edgecolor='None', loc='lower right')
        else:
            ax.legend(frameon=legend_frame, title=legend_title)
    else:
        if white:
            ax.legend(facecolor='white', edgecolor='None')
        else:
            ax.legend(frameon=legend_frame)
ax.tick_params(axis='x', direction='in', which='both', top=True, bottom=True, labelbottom=True, length=6)
ax.tick_params(axis='y', direction='in', which='both', left=True, right=True, labelleft=True, length=6)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='x', direction='in', which='minor', top=True, bottom=True, labelbottom=True, length=2)
ax.tick_params(axis='y', direction='in', which='minor', left=True, right=True, labelleft=True, length=2)
ax.grid(False)
if ylog:
    ax.set_yscale('log')
if xlog:
    ax.set_xscale('log')
if limit_x:
    ax.set_xlim(x_range)
    if minstep_x:
        ax.xaxis.set_major_locator(plt.MultipleLocator(step_x))
if limit_y:
    ax.set_ylim(y_range)
    if minstep_y:
        ax.yaxis.set_major_locator(plt.MultipleLocator(step_y))
if graph_title != 0 and graph_title != '' and graph_title != ' ':
    ax.set_title(graph_title, size=MEDIUM_SIZE)

fig.savefig(f'{name}.pdf', format='pdf', bbox_inches='tight')
plt.show()
