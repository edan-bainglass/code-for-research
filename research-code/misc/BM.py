import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import numpy as np

vols = []
energies = []
with open('EV', 'r') as f:
    for line in f.readlines():
        data = line.split()
        vols.append(float(data[0]))
        energies.append(float(data[1]))


def Birch_Murnaghan(parameters, V):

    E0, B0, BP, V0 = parameters

    E = E0 + 9 * V0 * B0 / 16 * (((V0 / V)**(2 / 3) - 1)**3 * BP +
                                 ((V0 / V)**(2 / 3) - 1)**2 *
                                 (6 - 4 * (V0 / V)**(2 / 3)))

    return E


def objective(pars, y, x):
    # we will minimize this function
    return y - Birch_Murnaghan(pars, x)


# initial guess of parameters
x0 = [min(energies), 1, 2, vols[energies.index(min(energies))]]

plsq = leastsq(objective, x0, args=(energies, vols))

print('E0 = {0:0.3f} eV'.format(plsq[0][0]))
print('V0 = {0:0.3f} Å³'.format(plsq[0][3]))
print('B0 = {0:0.3f} GPa'.format(plsq[0][1] * 160.2176621))
print('BP = {0:0.3f}'.format(plsq[0][2]))

plt.close('all')
plt.plot(vols, energies, 'ro')

# plot the fitted curve on top
x = np.linspace(min(vols), max(vols), 50)
y = Birch_Murnaghan(plsq[0], x)
plt.plot(x, y, 'k-')
plt.xlabel('Volume (Å³)')
plt.ylabel('Energy (eV)')
plt.tight_layout()
