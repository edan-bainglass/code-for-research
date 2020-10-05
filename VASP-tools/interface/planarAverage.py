import numpy as np

from ase.calculators.vasp import VaspChargeDensity
import matplotlib.pyplot as plt
import sys

from widgets import Slider, PremiumCheckButtons, Button

# sys.path.append(r"E:\Research\Code\Python\VaspTools")

# Globals
average, latticelength, ngridpts, idir, direction, ax, x, macro = None


def macroscopic_average(potential, periodicity, resolution):

    macro_average = np.zeros(shape=(len(potential)))
    period_points = int((periodicity / resolution))
    # Period points must be even
    if period_points % 2 != 0:
        period_points = period_points + 1

    length = len(potential)
    for i in range(length):
        start = i - int(period_points / 2)  # - L/2
        end = i + int(period_points / 2)  # + L/2
        if start < 0:
            start = start + length
            macro_average[i] = macro_average[i] \
                + sum(potential[0:end]) \
                + sum(potential[start:length])
        elif end >= length:
            end = end - length
            macro_average[i] = macro_average[i] \
                + sum(potential[start:length]) \
                + sum(potential[0:end])
        else:
            macro_average[i] = macro_average[i] + sum(potential[start:end])

        macro_average[i] = macro_average[i] / period_points

    print(f"Macroscopic Average = {np.average(macro_average):.3f}")

    return macro_average


def load(loc, axis="Z"):

    global average, latticelength, ngridpts, idir, direction

    LOCPOTfile = 'E:/Research/VASP Data/' + loc + '/LOCPOT'
    direction = axis

    vasp_charge = VaspChargeDensity(filename=LOCPOTfile)
    atoms = vasp_charge.atoms[-1]
    cell = atoms.cell
    potl = vasp_charge.chg[-1] * atoms.get_volume()
    latticelength = np.dot(cell, cell.T).diagonal()**0.5
    ngridpts = np.array(potl.shape)

    if direction == "X":
        idir = 0
    elif direction == "Y":
        idir = 1
    else:
        idir = 2

    average = np.zeros(ngridpts[idir], np.float)

    for ipt in range(ngridpts[idir]):
        if direction == "X":
            average[ipt] = potl[ipt, :, :].sum()
        elif direction == "Y":
            average[ipt] = potl[:, ipt, :].sum()
        else:
            average[ipt] = potl[:, :, ipt].sum()

    average /= ngridpts[(idir + 1) % 3] * ngridpts[(idir + 2) % 3]

    tot_avg_V = potl.sum() / ngridpts[0] / ngridpts[1] / ngridpts[2]
    print(f"\nTotal Average Electrostatic Potential = {tot_avg_V:.3f}\n")


def update(val):
    ax.lines.clear()
    ax.plot(x, average, 'k-')
    macro = macroscopic_average(average, val,
                                latticelength[idir] / ngridpts[idir])
    ax.plot(x, macro, 'r--')
    plt.draw()


def plot(periods=[], save=False):

    p = periods if periods else [latticelength[idir] / 2]

    global macro

    macro = macroscopic_average(average, p[0],
                                latticelength[idir] / ngridpts[idir])

    if interface:
        macro = macroscopic_average(macro, p[1],
                                    latticelength[idir] / ngridpts[idir])

    if printAverages:
        V_a = macro[int(a / latticelength[idir] * ngridpts[idir])]
        V_b = macro[int(b / latticelength[idir] * ngridpts[idir])]

        print(V_a)
        print(V_b)
        print(V_b - V_a)

    global x
    x = [
        i * latticelength[idir] / float(ngridpts[idir] - 1)
        for i in range(ngridpts[idir])
    ]

    fig = plt.figure(figsize=(12, 6))

    global ax

    ax = fig.add_subplot(111)
    ax.set_xlim(min(x), max(x))
    ax.set_xlabel(r'%s-axis ($\mathrm{\AA}$)' % direction)
    ax.set_ylabel('Electrostatic Potential (eV)')
    ax.plot(x, average, 'k-', x, macro, 'r--')
    plt.subplots_adjust(bottom=0.28)

    # if len(periods) > 1:
    #     macro = macroscopic_average(average, p[1],
    #                                 latticelength[idir] / ngridpts[idir])
    #     plt.plot(x, macro, 'b--')

    if save:
        plt.savefig(loc + '/macAvg', dpi=1000)


# MAIN

title = {'fontsize': 20, 'pad': 20}
axes = {'labelpad': 16, 'labelsize': 24}
ticks = {'labelsize': 20}
arrows = {'fontsize': 12, 'labelpad': 20}
font = {'family': 'times new roman'}
opt = {'font': font, 'xtick': ticks, 'ytick': ticks, 'axes': axes}
for i in opt:
    plt.rc(i, **opt[i])

locs = []

interface = True
printAverages = True
a = 13
b = 35

periods = [
    2.5,  # BVO
    4.4,  # VO
    # 8.4,
    # 16.7,
]

plt.close('all')

for loc in locs:
    load(loc, axis="Z")
    plot(periods=periods, save=True)

axSlider = plt.axes([0.15, 0.05, 0.75, 0.02])  # [left,bottom,width,height]
sPeriod = Slider(axSlider,
                 'Period',
                 0.001,
                 latticelength[idir],
                 valinit=0.001,
                 valstep=0.001,
                 valfmt="%1.3f")
sPeriod.on_changed(update)
sPeriod.label.set_size(14)

plt.show()
