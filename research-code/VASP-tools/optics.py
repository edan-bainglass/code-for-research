from numpy import array, linspace
from pymatgen.io.vasp.outputs import Vasprun
import matplotlib.pyplot as plt

# plot properties
plt.figure(figsize=(6, 5))
# gs = plt.GridSpec(1,2,wspace=0.3)
ticks = {
    'minor.width': 1.25,
    'major.width': 1.25,
    'major.pad': 6,
    'labelsize': 14
}
axes = {'lw': 1.25, 'labelweight': 'bold', 'labelpad': 10, 'labelsize': 14}
font = {'family': 'times new roman'}
opt = {'font': font, 'axes': axes, 'xtick': ticks, 'ytick': ticks}
for i in opt:
    plt.rc(i, **opt[i])
plt.subplots_adjust(bottom=.2, top=.9, left=.1, right=.9)

# get optical data
v = Vasprun(r'C:\Users\edanb\Desktop\vasprun.xml')
eps = v.dielectric  # eps 0,1,2 = energies, real, imaginary

# calculate absorption coefficient
k = [(0.5 * ((array(i)**2 + array(j)**2)**0.5 - array(i)))**0.5
     for i, j in zip(eps[1], eps[2])]
a = [2.0 / 3e10 * array(i) / 6.582e-16 * j for i, j in zip(eps[0], k)]  # 1/cm

r = 1 / 2  # for Tauc plot

xAxis = 'energy'
# xAxis = 'wavelength'

legend = ['xx', 'yy', 'zz', 'xy', 'yz', 'xz']
style1 = 'k.'
style2 = 'k.'
indices = [0]
alpha = [array(i)[indices] for i in a]  # select directions for plotting

# y-axis scales
scale1 = -5  # for absorption coefficient
scale2 = 10  # for Tauc plot

# plot absorption coefficient
plt.subplot(111)
if xAxis == 'energy':
    x1 = eps[0]
    y1 = [array(i) * 10**scale1 for i in alpha]
    plt.xlim(0, 2.5)
    plt.xlabel('eV')
    plt.ylim(0, 1)
    x = linspace(0, 2.5, 1000)
else:
    # remove 0.0 eV entry to avoid conversion issues
    alpha_wave = alpha.copy()
    alpha_wave.pop(0)
    x1 = [4.136e-15 * 3e17 / i for i in eps[0] if i != 0.0]
    y1 = [array(i) / 10**scale1 for i in alpha_wave]
    #    y1 = [array(i) for i in alpha_wave]
    plt.xlim(350, 950)
    plt.xlabel('Wavelength (nm)')
    #    plt.yscale('log')
    #    plt.ylim(1e0,1e5)
    plt.ylim(0, 1)
    x = linspace(350, 950, 1000)

i = 28
plt.plot(x1, y1, style1)
plt.plot(x, y1[i] + (y1[i + 1] - y1[i]) / (x1[i + 1] - x1[i]) * (x - x1[i]),
         'r')
plt.minorticks_on()
plt.ylabel('Absorption Coefficient (10$^{%d}$ cm$^{-1}$)' % scale1)
# plt.legend([legend[i] for i in indices])

plt.tight_layout()

# x2 = eps[0]
# y2 = [(array(i)*array(j))**(1/r)/10**scale2 for i,j in zip(alpha,eps[0])]

# plt.subplot(gs[1])
# i = 28
# x = linspace(1, 3.5, 1000)
# plt.plot(x2,y2,style2)
# plt.plot(x,y2[i]+(y2[i+1]-y2[i])/(x2[i+1]-x2[i])*(x-x2[i]),'r')
# plt.minorticks_on()
# plt.xlim(1,3.5)
# plt.xlabel('eV')
# plt.ylim(0,8)
# plt.ylabel(r'($\alpha$h$\nu$)$^{%d}$ (10$^{%d}$ eV$^{%d}$/cm$^{%d}$)'
#                                           % (1/r,scale2,1/r,1/r))
# plt.legend([legend[i] for i in indices])

plt.show()
