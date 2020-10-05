import openpyxl as xl
from openpyxl.utils import get_column_letter
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter, AutoMinorLocator
from matplotlib import rc, rcParams

rcParams['ytick.direction'] = 'out'
plt.rcParams['text.usetex'] = False

rc('font', **{'family': 'serif', 'serif': ['Arial']})
for axis in ['xtick', 'ytick']:
    rc(axis, labelsize=13)

wb = xl.load_workbook(r'C:\Users\edanb\Desktop\Bi2WO6_data.xlsx')
sheet = wb['Fig5']

x = [
    r'Bi$_2$WO$_6$',
    r'Sb 1',
    r'Sb 2',
    r'Sb 3',
    r'Sb 4',
    r'Sb 5',
    r'Sb 6',
    r'Sb$_2$WO$_6$',
]

plt.close('all')

y_axes = {'labelpad': 16, 'fontsize': 18}

# define axes limits
y_lim = (-3, 1)

fig, axes = plt.subplots(1,
                         1,
                         sharex=True,
                         figsize=(7.5, 3),
                         gridspec_kw={'hspace': 0.6})

axes = [axes]

# set spines
axes[0].spines['top'].set_visible(False)
for i in range(1):
    axes[i].spines['right'].set_visible(False)

# set x-axis ticks
axes[0].tick_params(axis='x', bottom=False, pad=14)
plt.setp(plt.gca().get_xticklabels(), rotation=0, va='baseline')

# set y-axis limits
ylim = [(1.8, 2)]

# set y-axis ticks
for i in range(1):
    axes[i].set_ylim(*ylim[i])
    axes[i].set_yticks(np.linspace(*(ylim[i] + (5, ))))
    axes[i].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axes[i].yaxis.set_minor_locator(AutoMinorLocator())

# set y-axis label
big_ax = fig.add_subplot(111, frameon=False)
big_ax.set_ylabel("Partial Charge (e)", **y_axes, c='k')
big_ax.tick_params(labelcolor='none', left=False, bottom=False)

handles = []

# Bi
y = [cell.value for cell in sheet['B':'B'] if cell.value is not None]
axes[0].plot(x[:-1], y[:-1], 'k', ls='--', lw=0.75, zorder=0)
handles.append(axes[0].scatter(x[:-1],
                               y[:-1],
                               marker='D',
                               s=40,
                               c='violet',
                               zorder=1,
                               label='Bi'))

# Sb
y = [cell.value for cell in sheet['C':'C'] if cell.value is not None]
axes[0].plot(x[1:], y[1:], 'k', ls='--', lw=0.75, zorder=0)
handles.append(axes[0].scatter(x[1:],
                               y[1:],
                               marker='s',
                               s=40,
                               c='y',
                               zorder=1,
                               label='Sb'))

labels = []
for h in handles:
    labels.append(h.get_label())
big_ax.legend(handles,
              labels,
              bbox_to_anchor=(0.25, 1),
              ncol=len(handles),
              frameon=False,
              columnspacing=0.25,
              fontsize=13,
              borderpad=0,
              borderaxespad=0,
              handletextpad=0.1)

plt.tight_layout()
plt.subplots_adjust(
    left=0.2,
    right=0.9,
)

plt.show()

plt.savefig(r'C:\Users\edanb\Desktop\Fig5.tiff', dpi=600, format='tiff')
