# thingies
import os
import time
import numpy as np
from multiprocessing import Pool
from sympy.solvers import solveset
from sympy import Symbol, S

starttime = time.time()
BGNDfilename = "\\0.25V_energy_nosmooth_HOPG.txt"
outputfilename = "output.txt"
energy = []
inverseenergy = []
ecounts = []
samplebiasBGND = 0.25
samplebiasDATA = 3.0
factor = samplebiasDATA - samplebiasBGND
intensityratio = 0.5

# Find inifinity channel
x = Symbol("x")
Energyfxn = 34.8699 - 0.131095 * x \
    + 0.0002061 * x**2 - 1.63837 * 10**-7 * x**3 \
    + 6.50396 * 10**-11 * x**4 - 1.0294 * 10**-14 * x**5  # HOPG
sols = solveset(Energyfxn, x, domain=S.Reals)
inf = int(sols.args[0])

# Grab BGND data
with open(os.getcwd() + BGNDfilename, "r") as file:
    line = file.readline()
    cnt = 0
    while cnt <= inf:
        line = file.readline()
        col1, col2 = line.split()
        energy.append(float(col1))
        inverseenergy.append(float(col1))
        ecounts.append(float(col2))
        cnt += 1

energy = np.array(energy) + factor
inverseenergy = (np.array(inverseenergy) + factor)**-0.5
ecounts = np.array(ecounts) * intensityratio


def calculateY(i, func=Energyfxn):
    return solveset(func - i, x, S.Reals).args[0]


if __name__ == '__main__':

    pool = Pool()
    channelnumbers = pool.map(calculateY, inverseenergy)
    pool.close()
    pool.join()

    print()
    print(channelnumbers)

    endtime = time.time()
    livetime = endtime - starttime
    print("\nIt took:", livetime, "seconds")
