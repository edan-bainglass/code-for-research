import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.patches import Rectangle
from matplotlib.ticker import AutoMinorLocator

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType
from pymatgen.electronic_structure.dos import Dos

import numpy as np
from scipy.interpolate import interp1d, splev, splrep
from collections import OrderedDict

locs = [
    r'new CBTO/Slabs/2/021/O7/dos_bands/',
    r'new CBTO/Slabs/1/110/O4/dos_bands/',
    r'new CBTO/Slabs/1/1-11/O8/dos_bands/',
    r'new CBTO/Slabs/2/112/O3/dos_bands/',
    r'new CBTO/Slabs/2/201/O8/dos_bands/',
]

dosruns = []
for loc in locs:
    dosruns.append(Vasprun("E:/Research/VASP Data/" + loc + "dos.xml"))

figs = []


def plotDOS(fig, ax, dosrun, orient='y', **kargs):
    def add_densities(density1, density2):
        return {
            spin: np.array(density1[spin]) + np.array(density2[spin])
            for spin in density1.keys()
        }

    def _get_orb_type(orb):
        try:
            return orb.orbital_type
        except AttributeError:
            return orb

    def add_site_spd_dos(spd, site):
        for orb, pdos in dosrun.complete_dos.pdos[site].items():
            orbital_type = _get_orb_type(orb)
            if orbital_type in spd:
                spd[orbital_type] = add_densities(spd[orbital_type], pdos)
            else:
                spd[orbital_type] = pdos

    def formatDos(x):

        l = np.array([0 if i == 0 else 1 for i in interp1d(e, x)(y)])
        start = []
        stop = []
        prevL = None
        for i, j in enumerate(l):
            if j == 1 and prevL == 0:
                start.append(i)
            if len(start) > 0 and (j == 0 or i == len(l) - 1) \
                    and prevL == 1:
                stop.append(i)
            prevL = j

        if smooth:
            func = splev(y, splrep(e, x)) * l
        else:
            func = x * l

        return func, [[i, j + 1] for i, j in zip(start, stop)]

    def plotIntervals(ax, x, y, o, intervals, t, c, a, lw, ls='-'):

        if orient != 'y':
            x, y = y, x

        if ax.get_label() == 'zoom':
            lw *= .75

        if (shade or t == 'Total' and ax.get_label() != 'zoom'):
            ax.fill_between(x,
                            y,
                            edgecolor=c / 255,
                            facecolor=c / 255,
                            alpha=0.5,
                            lw=lw,
                            ls=ls)
        else:
            if smooth:
                for i, j in intervals:
                    ax.plot(x[i:j][eval(o)],
                            y[i:j][eval(o)],
                            color=c / 255,
                            lw=lw,
                            ls=ls)
            else:
                ax.plot(x, y, color=c / 255, lw=lw, label=t, ls=ls)

        handles.append(
            Rectangle((0.0, 0.0),
                      .35,
                      .15,
                      ls=ls,
                      lw=lw,
                      facecolor=c / 255,
                      edgecolor=c / 255))

        labels.append(t)

    ###########################################################################

    _axes = [ax]
    if 'zoom' in kargs and kargs['zoom']:
        axz = ax.inset_axes(kargs['zoom'][2])
        axz.set_xlim(kargs['zoom'][0])
        axz.set_ylim(kargs['zoom'][1])
        axz.set_xticklabels('')
        axz.set_yticklabels('')
        axz.tick_params(bottom=False, left=False)
        axz.set_label('zoom')
        _axes.append(axz)
        ax.indicate_inset_zoom(axz)

    if 'spin' in kargs:
        spin = kargs['spin']
    else:
        if dosrun.is_spin:
            spin = [Spin.up, Spin.down]
        else:
            spin = [Spin.up]
    shade = None
    smooth = None
    if 'shade' in kargs:
        shade = kargs['shade']
    if 'smooth' in kargs:
        smooth = kargs['smooth']
    if 'ldos' in kargs:
        bulk = kargs['ldos'][2]
        ldos = True
        lower = kargs['ldos'][0]
        upper = kargs['ldos'][1]
    else:
        ldos = False
    legendOn = kargs['legend'] if 'legend' in kargs else False
    ymax = kargs['ymax'] if 'ymax' in kargs else False
    ldos = kargs['ldos'] if 'ldos' in kargs else False
    shift = kargs['shift'] if 'shift' in kargs else 0.
    plotTotal = kargs['plotTotal'] if 'plotTotal' in kargs else False
    plotOrbs = kargs['plotOrbs'] if 'plotOrbs' in kargs else False

    elements = kargs['elements'] if 'elements' in kargs else list(
        dict.fromkeys(dosrun.atomic_symbols))
    lw = kargs['lw'] if 'lw' in kargs else 1

    handles = []
    labels = []

    e = dosrun.tdos.energies - dosrun.tdos.get_cbm_vbm()[1] - shift
    y = e
    if smooth:
        y = np.arange(e[0], e[-1], 0.0005)

    for sp in spin:

        i = int(sp)
        if orient == 'y':
            o = 'x[i:j]>=0' if i == 1 else 'x[i:j]<=0'
        else:
            o = 'y[i:j]>=0' if i == 1 else 'y[i:j]<=0'

        # total DOS
        if plotTotal and not ldos:
            c = np.array([0, 0, 0])
            a = 0.1
            x, inter = formatDos(dosrun.tdos.densities[sp] * i)
            plotIntervals(ax, x, y, o, inter, 'Total', c, a, 0.5, ls=':')

        # partial DOS
        if plotOrbs:
            for name in elements:

                # spd projected DOS
                comp = dosrun.complete_dos
                if ldos:
                    if bulk:
                        sites = [
                            s for s in comp.structure.sites
                            if (s.c > lower and s.c < upper)
                            and s.species_string == name
                        ]
                    else:
                        sites = [
                            s for s in comp.structure.sites
                            if (s.c < lower or s.c > upper)
                            and s.species_string == name
                        ]
                    spd = dict()
                    for site in sites:
                        add_site_spd_dos(spd, site)
                    spd_dos = {
                        orb: Dos(comp.efermi, comp.energies, densities)
                        for orb, densities in spd.items()
                    }
                else:
                    spd_dos = comp.get_element_spd_dos(name)

                # plot selected orbitals
                for orb, props in orbitals[name].items():
                    c = np.array(props[1])
                    a = 0.25
                    x, inter = formatDos(
                        spd_dos[eval('OrbitalType.' + orb)].densities[sp] * i)

                    for ax in _axes:
                        plotIntervals(ax, x, y, o, inter, f"{name} {props[0]}",
                                      c, a, lw)

    if ymax:
        ax.set_ylim(top=ymax)

    if legendOn:
        legend = OrderedDict(zip(labels, handles))
        ax.legend(legend.values(), legend.keys(), loc=1, fontsize=16)


# %% setup

plt.close('all')

x_axes = {'labelpad': 14, 'fontsize': 14}
y_axes = {'labelpad': 0, 'fontsize': 16}

title_format = {'pad': 16, 'fontsize': 40}

rc('font', **{'family': 'serif', 'serif': ['Arial']})

fig, axes = plt.subplots(
    4,
    2,
    sharex='col',  # sharey='row',
    gridspec_kw={
        'wspace': 0.15,
        'hspace': 0.
    },
    figsize=(9, 8.5))

figs.append(fig)

ax_big = fig.add_subplot(111, frameon=False)
ax_big.set_xlabel(r"E - E$\mathrm{_f}$ (eV)",
                  labelpad=34,
                  fontsize=40,
                  color='k')
ax_big.set_ylabel("Density of States", labelpad=0, fontsize=40)
ax_big.tick_params(labelcolor='none',
                   top=False,
                   bottom=False,
                   left=False,
                   right=False,
                   pad=0)

xlim = (-0.5, 2)
ylim = (0, 100)

for ax, title in zip(axes[0], ['Bulk', 'Surface']):
    ax.set_title(title, color='w', **title_format)

for ax in axes[-1]:
    ax.set_xlim(xlim)
    ax.set_xticks(np.linspace(xlim[0], xlim[1], 11), minor=True)

for ax, z in zip(axes.flatten(), reversed(list(range(8)))):
    ax.set_zorder(z)
    ax.set_ylim(ylim)
    ax.tick_params(left=False,
                   labelsize=12,
                   labelcolor='k',
                   pad=8,
                   direction='out',
                   length=6)
    ax.tick_params(axis='x', which='minor', length=4)
    ax.set_yticklabels([])
    ax.axvline(0, 0, ylim[1], c='k', ls='--', lw=1.)

orbitals = {
    'Cu': {
        'd': ['3d', [180, 85, 15]]
    },
    'Bi': {
        'p': ['6p', [0, 100, 200]]
    },
    'W': {
        'd': ['5d', [100, 100, 100]]
    },
    'O': {
        'p': ['2p', [154, 0, 0]]
    }
}

elements = ['Cu', 'Bi', 'W', 'O']
bulk = [True, False]
legend = [False, True]
ymax = [ylim[1] * 3 / 5, ylim[1]]

for j in range(2):
    for i in range(4):

        dosopt = {
            # 'shade': True,
            'smooth': True,
            # 'zoom': [(-0.5, 0), (0, 20), [0.8, .5, .2, .5]],
            'legend': legend[j],
            'lw': 1.5,
            'ymax': ymax[j],
            'spin': [Spin.up],  # comment out for both spins
            'plotOrbs': True,
            'elements': [elements[i]],
            'ldos': [.196, .427, bulk[j]],
            # 'ldos': [.203, .439, bulk[j]],
            # 'ldos': [.202, .441, bulk[j]],
            # 'ldos': [.193, .419, bulk[j]],
            # 'ldos': [.182, .396, bulk[j]],
        }

        plotDOS(fig, axes[i, j], dosruns[0], orient='x', **dosopt)

        # break  # for testing

plt.subplots_adjust(
    top=0.9,
    bottom=0.2,
    left=0.12,
    right=0.96,
)

plt.show()
plt.savefig(r'C:\Users\edanb\Desktop\dispFig.jpeg', dpi=600, format='jpeg')
