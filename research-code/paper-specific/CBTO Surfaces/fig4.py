import openpyxl as xl
from openpyxl.utils import get_column_letter
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter, AutoMinorLocator
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['Arial']})

wb = xl.load_workbook(r'C:\Users\edanb\Desktop\data.xlsx')
sheet = wb['Energies']

x = [
    r'$(0\bar12)$',
    r'$(101)$',
    r'$(\bar111)$',
    r'$(1\bar21)$',
    r'$(2\bar10)$',
    # r'$(\bar230)$',
    r'$(012)$',
    # r'$(0\bar32)$',
    # r'$(210)$',
    r'$(11\bar1)$',
]

plt.close('all')

x_axes = {'labelpad': 12, 'fontsize': 18}
y_axes = {'labelpad': 12, 'fontsize': 18}

# define axes limits
y_lim = (0, 3)
sec_y_lim = (0, 0.06)

fig = plt.figure(figsize=(7.4, 5.8))
ax = plt.subplot(111)

# primary horizontal
ax.set_xlabel('Surface Planes', **x_axes)
ax.tick_params(axis='x', pad=12)
plt.setp(plt.gca().get_xticklabels(), rotation=0, va='center')

# primary vertical
# ax.set_ylabel('Surface Energy (eV/Å²)', **y_axes)
ax.set_ylabel('σ (J/m²)', **y_axes)
ax.set_ylim(*y_lim)
ax.set_yticks(np.linspace(*(y_lim + (7, ))))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(bottom=False, axis='both', labelsize=12)

# secondary vertical
sec = ax.twinx()
sec.set_ylabel(r'S (Å$^-$³)', **y_axes)
sec.set_ylim(*sec_y_lim)
sec.set_yticks(np.linspace(*(sec_y_lim + (8, ))))
sec.yaxis.set_minor_locator(AutoMinorLocator())
sec.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
sec.tick_params(axis='y', labelsize=12)

# unrelaxed
y = [cell.value for cell in sheet['E':'E'] if cell.value is not None]
ax.bar(x, y, .7, color='y', label='unrelaxed')

# relaxed
y = [cell.value for cell in sheet['F':'F'] if cell.value is not None]
ax.bar(x, y, .7, color='k', label='relaxed')

ax.legend(loc=2, fontsize=14)

# S/A
y = [cell.value for cell in sheet['C':'C'] if cell.value is not None]
sec.scatter(x, y, 100, 'darkred', '_', label='S')

# calculate and plot trendline to S/A
x = np.arange(7)
p = np.poly1d(np.polyfit(x, y, 1))
sec.plot(x, p(x), 'k', ls=':', lw=0.75)

sec.legend(loc=1, fontsize=12)

# Ceder
# sheet = wb['Ceder']

# # unrelaxed
# y = [cell.value for cell in sheet['A':'A'] if cell.value is not None]
# ax.scatter(x, y, color='b', zorder=5)

# # relaxed
# y = [cell.value for cell in sheet['B':'B'] if cell.value is not None]
# ax.scatter(x, y, color='r', zorder=5)

plt.tight_layout()

plt.show()

plt.savefig(r'C:\Users\edanb\Desktop\dispFig.jpeg', dpi=600, format='jpeg')
