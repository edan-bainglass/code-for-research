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

wb = xl.load_workbook(r'C:\Users\edanb\Desktop\defense.xlsx')
sheet = wb['bandAlign']

x = [r'Cu$_2$O', r'Fe$_2$O$_3$', r'WO$_3$', r'TiO$_2$', r'ZnO', r'Ga$_2$O$_3$']

plt.close('all')

y_axes = {'labelpad': 12, 'fontsize': 18}

# define axes limits
y_lim = (-4, 2)

fig = plt.figure(figsize=(8.5, 4))
ax = plt.subplot(111)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)

# primary horizontal
ax.tick_params(axis='x', pad=14, bottom=False)
ax.set_xlim(-0.5, len(x) - 0.5)
plt.setp(plt.gca().get_xticklabels(), rotation=0, va='baseline')

# primary vertical
# ax.set_ylabel(r'E$_{HSE}$ (eV)', **y_axes)
ax.set_ylim(*y_lim)
ax.set_yticks(np.linspace(*(y_lim + (7, ))))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y', pad=6)

# define marker
verts = [
    (-1.5, 0.1),
    (1.5, 0.1),
    (1.5, -0.1),
    (-1.5, -0.1),
]

# VBM
y = [cell.value for cell in sheet['B':'B'] if cell.value is not None]
ax.plot(x, y, 'k', marker=verts, markersize=50, lw=0)

# CBM
y = [cell.value for cell in sheet['C':'C'] if cell.value is not None]
ax.plot(x, y, 'k', marker=verts, markersize=50, lw=0)

# redox potentials
ax.hlines(0, -0.5, len(x) - 0.5, 'purple', linestyles='--', linewidth=1)
ax.hlines(-1.23, -0.5, len(x) - 0.5, 'g', linestyles='--', linewidth=1)

plt.tight_layout()

plt.show()

plt.savefig(r'C:\Users\edanb\Desktop\fig.tiff', dpi=600, format='tiff')
