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
sheet = wb['Fig8']

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

y_axes = {'labelpad': 12, 'fontsize': 18}

# define axes limits
y_lim = (-3, 1)

fig = plt.figure(figsize=(7.5, 4))
ax = plt.subplot(111)

# primary horizontal
ax.tick_params(axis='x', pad=14, bottom=False)
ax.set_xlim(-0.5, len(x) - 0.5)
plt.setp(plt.gca().get_xticklabels(), rotation=0, va='baseline')

# primary vertical
ax.set_ylabel('eV', **y_axes)
ax.set_ylim(*y_lim)
ax.set_yticks(np.linspace(*(y_lim + (5, ))))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y', pad=6)

# secondary vertical
ax_sec = ax.twinx()
ax_sec.set_ylim(*y_lim)
ax_sec.set_yticks(np.linspace(*(y_lim + (5, ))))
ax_sec.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax_sec.yaxis.set_minor_locator(AutoMinorLocator())
ax_sec.set_yticklabels(['' for i in range(len(ax.get_yticklabels()))])

# VBM
y = [cell.value for cell in sheet['D':'D'] if cell.value is not None]
ax.bar(x, y_lim[0], .7, color='steelblue', label='VBM')
ax.bar(x, y, .7, color='w', label='VBM_bg')

# CBM
y = [cell.value for cell in sheet['E':'E'] if cell.value is not None]
ax.bar(x, y_lim[1], .7, color='tab:olive', label='CBM')
ax.bar(x, y, .7, color='w', label='CBM_bg')

# Bi2WO6 level
l = 0.3
ax.hlines(-0.3 + l, -0.5, len(x) - 0.5, 'gray', ls='-', lw=0.5)
ax.hlines((0 + l, -1.23 + l), (-0.5, -0.5), (len(x) - 0.5, len(x) - 0.5),
          'k',
          ls='--',
          lw=0.75)

plt.tight_layout()

plt.show()

plt.savefig(r'C:\Users\edanb\Desktop\Fig8.tiff', dpi=600, format='tiff')
