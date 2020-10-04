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
sheet = wb['Fig6']

x = [
    r'Sb 1',
    r'Sb 2',
    r'Sb 3',
    r'Sb 4',
    r'Sb 5',
    r'Sb 6',
]

plt.close('all')

y_axes = {'labelpad': 10, 'fontsize': 18}

# define axes limits
y_lim = (2, 2.6)
sec_y_lim = (0, 0.2)

fig = plt.figure(figsize=(5, 3.6))
ax = plt.subplot(111)

# primary horizontal
ax.tick_params(axis='x', pad=14, bottom=False)
ax.set_xlim(-0.5, len(x) - 0.5)
plt.setp(plt.gca().get_xticklabels(), rotation=0, va='baseline')

# primary vertical
ax.set_ylabel('Bond Length (Å)', **y_axes)
ax.set_ylim(*y_lim)
ax.set_yticks(np.linspace(*(y_lim + (7, ))))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(axis='y', pad=6)

# secondary vertical
ax_sec = ax.twinx()
ax_sec.set_ylabel('D', **y_axes)
ax_sec.set_ylim(*sec_y_lim)
ax_sec.set_yticks(np.linspace(*(sec_y_lim + (5, ))))
ax_sec.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax_sec.yaxis.set_minor_locator(AutoMinorLocator())
ax_sec.tick_params(axis='y', pad=6)

handles = []

h = 0.1
w = 0.8
verts = [(-w, h), (w, h), (w, -h), (-w, -h)]

# B1
y = [cell.value for cell in sheet['B':'B'] if cell.value is not None]
handles.append(
    ax.scatter(x,
               y,
               marker=verts,
               s=100,
               c='olivedrab',
               zorder=0,
               label=r'B$\mathregular{_1}$'))

# B2
y = [cell.value for cell in sheet['C':'C'] if cell.value is not None]
handles.append(
    ax.scatter(x,
               y,
               marker=verts,
               s=100,
               c='royalblue',
               zorder=0,
               label=r'B$\mathregular{_2}$'))

# B3
y = [cell.value for cell in sheet['D':'D'] if cell.value is not None]
handles.append(
    ax.scatter(x,
               y,
               marker=verts,
               s=100,
               c='goldenrod',
               zorder=0,
               label=r'B$\mathregular{_3}$'))

# D
y = [cell.value for cell in sheet['E':'E'] if cell.value is not None]
ax_sec.plot(x, y, 'k', ls='--', lw=0.75, zorder=0)
handles.append(
    ax_sec.scatter(x, y, marker='X', s=50, c='maroon', zorder=1, label='D'))

labels = []
for h in handles:
    labels.append(h.get_label())
ax.legend(handles,
          labels,
          ncol=4,
          loc=9,
          frameon=False,
          columnspacing=0.1,
          fontsize=13,
          borderpad=1,
          borderaxespad=0,
          handletextpad=0.1)
labels = ax.get_legend().get_texts()
labels[0].set_position((0, -3))
labels[1].set_position((0, -3))
labels[2].set_position((0, -3))
labels[3].set_position((0, -2))

plt.tight_layout()

plt.show()

plt.savefig(r'C:\Users\edanb\Desktop\Fig6.tiff', dpi=600, format='tiff')
