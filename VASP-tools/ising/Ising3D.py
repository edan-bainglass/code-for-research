import numpy as np
from numpy.random import random, randint

from math import exp

import matplotlib.pyplot as plt

# Globals
S, E = None, None


# randomly assign spin values to lattice sites
def init(L, spin):
    if spin == 1:
        return 2 * randint(2, size=(L, L, L)) - 1
    else:
        return np.full([L, L, L], 1)


# sweep lattice and perform requested operation
def sweep(func, lat):
    # define subroutines
    def nearSum():
        return lat[(i - 1) % L, j, k] + lat[(i + 1) % L, j, k] + lat[
            i, (j - 1) % L, k] + lat[i, (j + 1) % L,
                                     k] + lat[i, j,
                                              (k - 1) % L] + lat[i, j,
                                                                 (k + 1) % L]

    def flip_or_not():
        def calcE():
            return -1 * lat[i, j, k] * nearSum()

        e1 = calcE()  # calculate local energy
        lat[i, j, k] *= -1  # flip spin at current site
        e2 = calcE()  # recalculate local energy
        eDif = e2 - e1  # calculate energy difference
        if eDif < 0:  # determine if flip is energetically favorable
            pass
        elif random() <= exp(-eDif / (1. * t)):
            pass
        else:
            lat[i, j, k] *= -1

    def calc():
        global E, S
        E += -lat[i, j, k] * nearSum()
        S += lat[i, j, k]

    # perform sweep
    for i in range(L):
        for j in range(L):
            for k in range(L):
                # execute requested operation
                eval(func)


# display lattice (FIX THIS!!!)
def showLat(lat, t1, t2):
    if t1 <= t <= t2:
        fig_lat = plt.figure()
        for i in range(3):
            ax = fig_lat.add_subplot(1, 3, i + 1)
            ax.imshow(lat[:, :, i])
            ax.axis('off')
            ax.title(f'{i} @ T = {round(t,2)}', fontsize=8)
            ax.tight_layout()


# execute Ising algorithm
def Ising(t, init_spin):
    global S, E
    # initialize lattice
    lat = init(L, init_spin)
    # thermalize lattice
    for n in range(intswp):
        sweep('flip_or_not()', lat)
    # calculate magnetization and energy for ensemble
    avgMag = avgE = 0.0
    for n in range(sweeps):
        sweep('flip_or_not()', lat)
        S = E = 0.0
        sweep('calc()', lat)
        avgMag = avgMag + S
        avgE = avgE + E / 2
        del S, E
    # normalize magnetization and energy averages
    Mag.append(abs(avgMag) / sweeps / L**3)
    Enrg.append(avgE / sweeps)
    del avgMag, avgE

    return lat


# Main Program #

# set parameters
L = 4  # lattice size
intswp = 1  # number of sweep to equilibriate lattice
sweeps = 100  # number of sweeps for calculations
TR1 = np.linspace(0.1, 3.5, 10)
TR2 = np.linspace(3.6, 5.5, 20)
TR3 = np.linspace(5.6, 8.0, 10)
Mag, Enrg = [], []
T = np.concatenate((TR1, TR2, TR3))

for t in T:
    lat = Ising(t,
                0)  # Initial lattice: 0 (all spin up); 1 (random spin mixture)

# plot results
fig = plt.figure(figsize=(7, 9))
ax1 = fig.add_subplot(211)
plt.title("3D Ising Model", fontsize=20, pad=20)
plt.scatter(T, Mag, color='IndianRed')
plt.ylabel("Magnetization ", fontsize=20, labelpad=10)
plt.axis('tight')
ax2 = fig.add_subplot(212)
plt.scatter(T, Enrg)
plt.xlabel("Temperature (T)", fontsize=20, labelpad=10)
plt.ylabel("Energy ", fontsize=20, labelpad=10)
plt.axis('tight')
