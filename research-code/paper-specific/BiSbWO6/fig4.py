import openpyxl as xl
from openpyxl.utils import get_column_letter
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter, AutoMinorLocator
from matplotlib import rc, rcParams

rcParams['ytick.direction'] = 'in'
rc('font', **{'family': 'serif', 'serif': ['Arial']})
for axis in ['xtick', 'ytick']:
    rc(axis, labelsize=13)

wb = xl.load_workbook(r'C:\Users\edanb\Desktop\Bi2WO6_data.xlsx')
sheet = wb['Fig4']

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

y_axes = {'labelpad': 10, 'fontsize': 18}

# define axes limits
y_lim = (-15.65, -15.35)

fig = plt.figure(figsize=(7.5, 3))
ax = plt.subplot(111)

# primary horizontal
ax.tick_params(axis='x', pad=14, bottom=False)
ax.set_xlim(-0.5, len(x) - 0.5)
plt.setp(plt.gca().get_xticklabels(), rotation=0, va='baseline')

# primary vertical
ax.set_ylabel('eV', **y_axes)
ax.set_ylim(*y_lim)
ax.set_yticks(np.linspace(*(y_lim + (4, ))))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y', pad=6)

# secondary vertical
ax_sec = ax.twinx()
ax_sec.set_ylim(*y_lim)
ax_sec.set_yticks(np.linspace(*(y_lim + (4, ))))
ax_sec.yaxis.set_minor_locator(AutoMinorLocator())
ax_sec.set_yticklabels(['' for i in range(len(ax.get_yticklabels()))])

# delH
y = [cell.value for cell in sheet['B':'B'] if cell.value is not None]
ax.plot(x, y, 'k', ls='--', lw=0.75, zorder=1)
ax.scatter(x, y, marker='o', s=30, c='steelblue', zorder=2, label='ΔH')
ax.hlines((y[4], y[6]), (-0.5, -0.5), (len(y) - 0.5, len(y) - 0.5),
          lw=0.5,
          ls='-',
          color='gray',
          zorder=0)

# delH + TS
y = [cell.value for cell in sheet['C':'C'] if cell.value is not None]
ax.scatter(x, y, marker=11, s=20, c='olivedrab', zorder=1, label='ΔH + TS')

ax.legend(bbox_to_anchor=((0.23, 1.)),
          fontsize=13,
          frameon=False,
          handletextpad=0.1,
          borderpad=0)
labels = ax.get_legend().get_texts()
labels[0].set_position((0, -2))
labels[1].set_position((0, -5))

plt.tight_layout()
plt.subplots_adjust(
    left=0.2,
    right=0.9,
)

plt.show()

plt.savefig(r'C:\Users\edanb\Desktop\Fig4.tiff', dpi=600, format='tiff')
