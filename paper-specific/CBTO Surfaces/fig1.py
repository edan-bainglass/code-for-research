import openpyxl as xl
from openpyxl.utils import get_column_letter
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter, AutoMinorLocator
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Arial']})

wb = xl.load_workbook(r'C:\Users\edanb\Desktop\data.xlsx')
sheet = wb['110']

plt.close('all')

x_axes = {'labelpad': 12, 'fontsize': 18}
y_axes = {'labelpad': 12, 'fontsize': 18}

# define axes limits
x_lim = (3, 21)
y_lim = (0.285, 0.305)
sec_y_lim = (0, 0.25)

fig = plt.figure(figsize=(7.4, 5.8))
ax = plt.subplot(111)

# primary horizontal
ax.set_xlabel('N', **x_axes)
ax.tick_params(axis='x', pad=12)
ax.set_xlim(*x_lim)
ax.set_xticks(np.linspace(*(x_lim + (10, ))))
ax.tick_params(axis='both', labelsize=12)

# primary vertical
ax.set_ylabel('σ (J/m²)', **y_axes)
ax.set_ylim(*y_lim)
ax.set_yticks(np.linspace(*(y_lim + (5, ))))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
ax.yaxis.set_minor_locator(AutoMinorLocator())

# secondary vertical
sec = ax.twinx()
sec.set_ylabel(r'Average Ion Displacement (Å)', **y_axes)
sec.set_ylim(*sec_y_lim)
sec.set_yticks(np.linspace(*(sec_y_lim + (6, ))))
sec.yaxis.set_minor_locator(AutoMinorLocator())
sec.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
sec.tick_params(axis='y', labelsize=12)

x = [cell.value for cell in sheet['A':'A'] if cell.value is not None]

y = [cell.value for cell in sheet['B':'B'] if cell.value is not None]
ax.plot(x[2:],
        y,
        ls='--',
        lw=1,
        marker='x',
        mfc='k',
        mec='k',
        mew=1,
        label='Linear')

y = [cell.value for cell in sheet['C':'C'] if cell.value is not None]
ax.plot(x[2:],
        y,
        ls='--',
        lw=1,
        marker='x',
        mfc='k',
        mec='k',
        mew=1,
        label='Bulk')

y = [cell.value for cell in sheet['D':'D'] if cell.value is not None]
ax.plot(x[2:],
        y,
        ls='--',
        lw=1,
        marker='x',
        mfc='k',
        mec='k',
        mew=1,
        label='Ceder')

y = [cell.value for cell in sheet['E':'E'] if cell.value is not None]
sec.plot(x,
         y,
         'k',
         ls='--',
         lw=1,
         marker='.',
         markersize=10,
         mfc='purple',
         mec='purple',
         mew=1,
         label='Displacement')

ax.legend(loc=2, fontsize=14)

sec.legend(loc=1, fontsize=12)

plt.tight_layout()

plt.show()

plt.savefig(r'C:\Users\edanb\Desktop\dispFig.jpeg', dpi=600, format='jpeg')
