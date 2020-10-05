import numpy as np
from numpy.random import random, randint

from math import exp

import matplotlib.pyplot as plt

# Globals
S, E = None, None


# randomly assign spin values to lattice sites
def init(L, spin):
    if spin == 1:
        return 2 * randint(2, size=(L, L)) - 1
    else:
        return np.full([L, L], 1)


# sweep lattice and perform requested operation
def sweep(func, lat):
    # define subroutines
    def nearSum():
        return lat[(i - 1) % L, j] + lat[
            (i + 1) % L, j] + lat[i, (j - 1) % L] + lat[i, (j + 1) % L]

    def flip_or_not():
        def calcE():
            return -1 * lat[i, j] * nearSum()

        e1 = calcE()  # calculate local energy
        lat[i, j] *= -1  # flip spin at current site
        e2 = calcE()  # recalculate local energy
        eDif = e2 - e1  # calculate energy difference
        if eDif >= 0:  # determine if flip is energetically favorable
            if random() > exp(-eDif / (1 * t)): lat[i, j] *= -1

    def calc():
        global E, S
        E = E - lat[i, j] * nearSum()
        S = S + lat[i, j]

    # perform sweep
    for i in range(L):
        for j in range(L):
            # execute requested operation
            eval(func)


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
    Mag.append(abs(avgMag) / sweeps / L**2)
    Enrg.append(avgE / sweeps)
    del avgMag, avgE


# Main Program #

# set parameters
L = 5  # lattice size
intswp = 1  # number of sweep to equilibriate lattice
sweeps = 100  # number of sweeps for calculations
T1 = [0.1, 1.5, 10]
T2 = [1.6, 3.0, 30]
T3 = [3.1, 5.0, 10]
T, Mag, Enrg = [], [], []

for r in [T1, T2, T3]:
    temp = np.linspace(r[0], r[1], r[2])  # define temperature range
    T.extend(temp)
    # sweep through the lattice for each temperature step
    for t in temp:
        Ising(t,
              0)  # Initial lattice: 0 (all spin up); 1 (random spin mixture)

# plot results
fig = plt.figure(figsize=(7, 9))
ax1 = fig.add_subplot(211)
plt.title("2D Ising Model", fontsize=20, pad=20)
plt.scatter(T, Mag, color='IndianRed')
plt.ylabel("Magnetization ", fontsize=20, labelpad=10)
plt.axis('tight')
ax2 = fig.add_subplot(212)
plt.scatter(T, Enrg)
plt.xlabel("Temperature (T)", fontsize=20, labelpad=10)
plt.ylabel("Energy ", fontsize=20, labelpad=10)
plt.axis('tight')
