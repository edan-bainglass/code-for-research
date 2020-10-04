import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from scipy.interpolate import interp1d
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin

# TODO - comment out the line below when running. This code is only used when a band structure plot is running, so fig and bands are stored in memory
fig = bands = None

c = 7.6305025868
edge = 'CBM'
gap = True
band = 51
diff = -1
distances = bands.distance
energies = np.copy(bands.bands[Spin.up][band if edge == 'VBM' else band + 1])

if gap:
    energies -= bands.get_vbm()['energy']
    fig.axes[0].collections.clear()
else:
    energies -= bands.efermi
    fig.axes[0].collections.clear()

# fig.axes[0].lines = fig.axes[0].lines[:-1]
plt.draw()

if edge == 'VBM':
    op = 'max'
elif edge == 'CBM':
    op = 'min'

i = np.where(energies == eval(op + '(energies)'))[0][0]  # central point
# i = 0  # central point
l = 5  # number of points away from central point
n = 4  # number of interpolated points

if diff == 1:
    x, y = distances[i:i + l + 1], energies[i:i + l + 1]
elif diff == -1:
    x, y = distances[i - l:i + 1], energies[i - l:i + 1]

p = interp1d(x, y, kind='quadratic')
xi = np.linspace(x[0], x[-1], l * (n + 1) + 1)
yi = p(xi)

plt.scatter(x, y, s=40, c='y', label='points')
plt.scatter(xi, yi, c='k', s=5, label='interpol')

# %% Parabolic fit

# s = 1
# step = xi[-1] - xi[-2]
# xp = np.array(xi.tolist() + [distances[i] + j *
#                              step for j in range(1, l*(n+1) + 1)])
# yp = np.array(yi.tolist() + sorted([yi[j]
#                                     for j in range(l*(n+1))], reverse=True))
# p = np.poly1d(np.polyfit(xp[s:-s], yp[s:-s], 2))
# plt.plot(xi, b(xi))

# %% Calculate effM

i = np.where(yi == eval(op +
                        '(yi)'))[0][0]  # central point of finite difference
# i = 0  # central point of finite difference
h = 2  # index step size for finite difference (h)

if diff == 1:  # forward
    curv = (yi[i + 2 * h] - 2 * yi[i + h] + yi[i]) / (xi[i + h] - xi[i])**2
    print("forward  = [ {0:0.3f} eV ; {1:0.3f} 2pi/A ]".format(
        abs(yi[i + h] - yi[i]), xi[i + h] - xi[i]))
elif diff == -1:  # backward
    curv = (yi[i] - 2 * yi[i - h] + yi[i - 2 * h]) / (xi[i - h] - xi[i])**2
    print("backward = [ {0:0.3f} eV ; {1:0.3f} 2pi/A ]".format(
        abs(yi[i] - yi[i - h]), xi[i] - xi[i - h]))

# central difference

# 3-point stencil
# curv = (yi[i+h] - 2*yi[i] + yi[i-h]) / (xi[i] - xi[i-h])**2

# 5-point stencil
# curv = (-1 * yi[i+2*h] + 16*yi[i+h] + -30*yi[i] + 16 *
#         yi[i-h] - yi[i-2*h]) / 12 / (xi[i] - xi[i-h])**2

# print("forward  = [ {0:0.3f} eV ; {1:0.3f} 2pi/A ]".format(
#     abs(yi[i + h] - yi[i]), xi[i + h] - xi[i]))
# print("backward = [ {0:0.3f} eV ; {1:0.3f} 2pi/A ]".format(
#     abs(yi[i] - yi[i - h]), xi[i] - xi[i - h]))

print("\neff_m = {0:0.3f}".format(1 / curv * c))

plt.show()
