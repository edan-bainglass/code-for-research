import numpy as np

from pymatgen.electronic_structure.core import Spin, OrbitalType

from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga

from pymatgen.io.vasp.outputs import BSVasprun, Vasprun

from pymatgen.electronic_structure.dos import Dos

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import MultipleLocator, FuncFormatter, AutoMinorLocator
from matplotlib.collections import LineCollection

from scipy.interpolate import splev, splrep

from matplotlib import cm

# define global namespace variables
origLocs = bandruns = dosruns = bands = fig = ax1 = _axes = None
prim, path, newPath = None

# FUNCTIONS


class My_Axes_Y(mpl.axes.Axes):
    name = 'My_Axes_Y'

    def drag_pan(self, button, key, x, y):
        mpl.axes.Axes.drag_pan(self, button, 'y', x, y)


class My_Axes_X(mpl.axes.Axes):
    name = 'My_Axes_X'

    def drag_pan(self, button, key, x, y):
        mpl.axes.Axes.drag_pan(self, button, 'x', x, y)


def add_densities(density1, density2):
    return {
        spin: np.array(density1[spin]) + np.array(density2[spin])
        for spin in density1.keys()
    }


def colorPick(event):

    global clickColor

    if event.inaxes:
        im = event.inaxes.get_images()[0]
        colors = im.cmap(im.norm(im.get_array().data))
        clickColor = colors[int(event.ydata)][0][0:3]
        clickColor = [int(i * 255) for i in clickColor.tolist()]
    plt.draw()


def show_cmaps(names=None):

    f = plt.figure(figsize=(3, 5))

    a = np.outer(np.arange(0, 1, 0.001), np.ones(10))

    maps = [m for m in cm.datad if not m.endswith("_r")]
    maps.sort()

    if names:
        maps = names if isinstance(names, list) else [names]
    l = len(maps)

    for i, m in enumerate(maps):
        plt.subplot(1, l, i + 1, projection='My_Axes_Y')
        plt.axis("off")
        plt.imshow(a, aspect='auto', cmap=cm.get_cmap(m), origin="lower")
        # plt.title(m,rotation=90,fontsize=10)

    f.canvas.mpl_connect('button_press_event', colorPick)


def plot(locs,
         parse=False,
         figsize=None,
         dpi=300,
         fmt='jpeg',
         lw=1.0,
         save=False):

    for p in (My_Axes_X, My_Axes_Y):
        mpl.projections.register_projection(p)

    def setDos(axis, minor, spin=[Spin.up], lim=None):
        """ Set up plotting area for DOS. If bands are plotted,
            set y-axis to vertical"""
        def spines(zero, hide):
            for z in zero:
                axD.spines[z].set_position('zero')
            for h in hide:
                axD.spines[h].set_visible(False)

        if len(spin) == 1:
            if spin[0].value == 1:
                # axD.set_title(u'DOS $\u2191$', **title)
                if axis == 'y':
                    axD.set_xlim(0, lim[1])
                    spines(zero=['left'], hide=['right', 'top', 'bottom'])
                    axD.yaxis.tick_left()
                else:
                    axD.set_ylim(0, lim[1])
                    spines(zero=['bottom'], hide=['left', 'right', 'top'])
                    axD.xaxis.tick_bottom()
                    # axD.tick_params(bottom=False)
            else:
                axD.set_title(u'DOS $\u2193$', **title)
                if axis == 'y':
                    axD.set_xlim(lim[0], 0)
                    spines(zero=['right', 'top'], hide=['left', 'bottom'])
                    axD.yaxis.tick_right()
                else:
                    axD.set_ylim(lim[0], 0)
                    spines(zero=['right', 'top'], hide=['left', 'bottom'])
                    axD.xaxis.tick_top()
                    axD.xaxis.set_label_coords(0.5, -0.1)
        else:
            axD.set_title('DOS', **title)
            if axis == 'y':
                axD.set_xlim(lim)
                spines(zero=['left', 'top'], hide=['right', 'bottom'])
                axD.set_xlabel(u'$\u2193$     $\u2191$', **arrows)
            else:
                axD.set_ylim(lim)
                spines(zero=['right', 'bottom'], hide=['left', 'top'])
                axD.set_ylabel(u'$\u2190$     $\u2192$', **arrows)
            axD.xaxis.set_label_coords(0.5, -0.04)

        # define tick labels for density (in arbitrary units)
        if 'DOSticks' in dosopt:

            def func(x, pos):
                return "" if np.isclose(x, 0) else "%d" % x \
                    if np.equal(np.mod(x, 1), 0) else x
        else:

            def func(x, pos):
                return ""

        if axis == 'y':  # if bands are plotted...
            axD.axhline(0, lw=1.5, color='k', ls='-')
            axD.xaxis.set_visible(False)
            axD.yaxis.set_minor_locator(MultipleLocator(minor))
            # axD.yaxis.set_visible(False)
        else:
            # axD.yaxis.tick_right()
            axD.set_xlabel(r"E - E$\mathrm{_f}$ (eV)", **axes)
            axD.yaxis.set_major_formatter(FuncFormatter(func))
            axD.xaxis.set_minor_locator(MultipleLocator(minor))
            plt.tight_layout()

        # axD.minorticks_on()
        axD.tick_params(axis='y', which='major', length=0)
        axD.tick_params(axis='y', which='minor', length=0)
        # axD.tick_params(axis='x', which='major', length=6)
        # axD.tick_params(axis='x', which='minor', length=3)

        global _axes
        _axes = [axD]
        if 'zoom' in dosopt and dosopt['zoom']:
            if len(dosopt['zoom']) > 2:
                axz = axD.inset_axes(dosopt['zoom'][2])
            elif axis == 'y':
                axz = axD.inset_axes([0.5, 0.7, 0.5, 0.3])
            else:
                axz = axD.inset_axes([0.8, 0.55, 0.2, 0.45])
            axz.set_xlim(dosopt['zoom'][0])
            axz.set_ylim(dosopt['zoom'][1])
            axz.set_xticklabels('')
            axz.set_yticklabels('')
            axz.set_label('zoom')
            _axes.append(axz)
            axD.indicate_inset_zoom(axz)

    ###########################################################################

    # set general linewidth for plot borders
    for i in [
            'axes.linewidth', 'xtick.major.width', 'xtick.minor.width',
            'ytick.major.width', 'ytick.minor.width'
    ]:
        mpl.rcParams[i] = lw

    # define plot axes
    global origLocs, bandruns, dosruns, bands, clickWidth, fig

    if parse:
        bandruns = []
        dosruns = []
        origLocs = locs.copy()

    if colorChoice[0]:
        clickWidth = colorChoice[1]
        show_cmaps('gist_ncar')

    for i, loc in enumerate(locs):

        floc = "E:/Research/VASP Data/" + loc

        proj = bandopt['projections'] if 'projections' in bandopt else False

        if 'bands' in plots or 'BZ' in plots:

            if parse:
                bandruns.append(
                    BSVasprun(floc + "/bands.xml", parse_projected_eigen=proj))
            else:
                i = origLocs.index(loc)

            bands = bandruns[i].get_band_structure(floc + "/KPOINTS",
                                                   line_mode=True,
                                                   efermi=bandruns[i].efermi)

        if 'bands' in plots:
            figsize = figsize if figsize else (5, 7)
            fig = plt.figure(figsize=figsize)
            gs = fig.add_gridspec(1, 5, wspace=0.) \
                if 'dos' in plots else fig.add_gridspec(1, 2)

            axB = fig.add_subplot(gs[0:2], projection='My_Axes_Y')
            # axB.set_title('Band Structure', **title)
            # axB.set_xlabel(r'Symmetry Points', c='k', **axes)
            axB.set_ylabel(r"E - E$\mathrm{_f}$ (eV)", **axes)
            # axB.yaxis.tick_right()
            # axB.yaxis.set_label_position("right")
            axB.tick_params(axis='x', which='both', length=0, pad=10)
            axB.tick_params(axis='y',
                            which='both',
                            direction='out',
                            pad=6,
                            labelsize=ticks['fontsize'],
                            length=6)
            axB.set_yticks(np.linspace(-1.5, 2.5, 5), minor=True)
            axB.yaxis.set_minor_locator(AutoMinorLocator())
            axB.tick_params(axis='y', which='minor', length=3)

            plotBands(axB, bandruns[i], **bandopt)

        if 'dos' in plots:

            if parse:
                dosruns.append(Vasprun(floc + r"\..\dos.xml"))
            else:
                i = origLocs.index(loc)

            if 'spin' in dosopt:
                spin = dosopt['spin']
            else:
                if dosruns[i].is_spin:
                    spin = [Spin.up, Spin.down]
                else:
                    spin = [Spin.up]

            xlim = dosopt['xlim'] if 'xlim' in dosopt else (-4, 4)
            ylim = dosopt['ylim'] if 'ylim' in dosopt else (-15, 15)
            shift = dosopt['shift'] if 'shift' in dosopt else 0.
            if 'ylim' in dosopt:
                ylim = dosopt['ylim']
            else:
                e = dosruns[-1].tdos.energies - \
                    dosruns[-1].tdos.get_cbm_vbm()[1] - shift
                ind = [
                    e.tolist().index(i)
                    for i in e[(e > xlim[0]) & (e < xlim[1])]
                ]
                lim = max([
                    max(a[ind[0]:ind[-1]])
                    for a in dosruns[-1].tdos.densities.values()
                ])
                ylim = (-lim, lim)

            vertical = dosopt['vertical'] if 'vertical' in dosopt else None

            if 'bands' in plots:
                axD = fig.add_subplot(gs[2:],
                                      sharey=axB,
                                      projection='My_Axes_Y')
                axD.tick_params(labelleft=False)
                setDos('y', 0.5, spin=spin, lim=ylim)
                plotDOS(fig, dosruns[i], **dosopt)
            else:
                if figsize:
                    figsize = figsize
                elif vertical:
                    figsize = (5, 10)
                else:
                    figsize = (14, 5)
                fig = plt.figure(figsize=figsize)

                if vertical:
                    # axD = fig.add_subplot(111,projection='My_Axes_Y')
                    axD = fig.add_subplot(111)
                    axD.set_ylim(xlim)
                    setDos('y', 0.25, spin=spin, lim=ylim)
                    plotDOS(fig, dosruns[i], orient='y', **dosopt)
                else:
                    # axD = fig.add_subplot(111, projection='My_Axes_X')
                    axD = fig.add_subplot(111)
                    axD.set_xlim(xlim)
                    setDos('x', 0.25, spin=spin, lim=ylim)
                    plotDOS(fig, dosruns[i], orient='x', **dosopt)

        plt.tight_layout()

        if 'BZ' in plots:
            plotBZ()

        if save:
            folder = r'C:/Users/edanb/Desktop'
            # folder = floc
            plt.savefig(folder + '/ ' + str(i) + '_%s' % '_'.join(plots) +
                        '.' + str(fmt),
                        format='tiff',
                        dpi=dpi)

        plt.show()


def plotBands(ax, bandrun, **kargs):
    """ Plot band structure """
    def setTicks(ax, dk):
        """ Collect symmetry k-point labels and positions """

        if k.label:
            if k.label != prevK[0]:
                if k.label == branch['name'].split('-')[0]:
                    ax.axvline(prevK[1], color='k', ls='-', lw=1.5)
                    if k.label != k0:
                        dk += segWidth
                    ax.axvline(d + dk, color='k', ls='-', lw=1.5)
                if k.label.startswith("\\") or k.label.find("_") != -1:
                    kLabels['label'].append('$' + k.label + '$')
                else:
                    kLabels['label'].append(k.label)
                kLabels['distance'].append(d + dk)
            else:
                ax.axvline(d + dk, color='k', ls='-', lw=0.4)

        return dk

    def plotBranch(ax, d, e, tol, i, c, ls):
        """ Plot each band over each k-point branch """

        global ax1
        ax1 = ax

        def contrib(species, data, el, orb):
            """ Calculate orbital contributions """

            # SPEED THIS UP BY UTILIZING NUMPY'S VECTORIZED CODE

            p = [i for i, s in enumerate(species) if s == el]

            contrib = []
            c = np.array(orbitals[el][orb][1]) / 255
            o = eval('OrbitalType.' + orb + '.value')
            for b in range(*bandRange):
                for k in range(len(d) - 1):
                    contrib.append(
                        np.concatenate([
                            c,
                            [
                                np.average([
                                    data[b, k + i, o, p[0]:p[-1] + 1].sum()
                                    for i in range(2)
                                ])
                            ]
                        ]))
                    if smooth:
                        for i in range(int((tol - len(d)) / len(d))):
                            contrib.append(contrib[-1])
                if smooth:
                    for i in range(int((tol - len(d)) / len(d))):
                        contrib.append(contrib[-1])
            return contrib

        segs = []  # store segmented bands if plotting projections

        if smooth:

            if not tol:
                tol = len(d) - 1

            for b in range(*bandRange):

                tck = splrep(d, e[:, b])
                step = (d[-1] - d[0]) / (tol - 1)
                xs = np.array([x * step + d[0] for x in range(tol)])
                ys = np.array(
                    [splev(x * step + d[0], tck, der=0) for x in range(tol)])
                for y in ys:
                    if np.isnan(y):
                        print(" -> Smoothing failure :( ")
                        break

                if projections:
                    pts = np.array([xs, ys]).T.reshape(-1, 1, 2)
                    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
                    segs.append(seg)
                else:
                    # if 44 <= b <= 45:
                    #     ax.plot(xs, ys, 'dodgerblue', lw=2, alpha=1.0,
                    #             picker=2, label=b + 1)
                    # if 46 <= b <= 47:
                    #     ax.plot(xs, ys, 'darkorange', lw=2, alpha=1.0,
                    #             picker=2, label=b + 1)
                    # if 116 <= b <= 117:
                    #     ax.plot(xs, ys, 'forestgreen', lw=2, alpha=1.0,
                    #             picker=2, label=b + 1)
                    # else:
                    ax.plot(xs,
                            ys,
                            c,
                            lw=lw,
                            alpha=alpha,
                            picker=2,
                            label=b + 1)

        else:
            if projections:
                for b in range(*bandRange):
                    pts = np.array([d, e[:, b]]).T.reshape(-1, 1, 2)
                    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)
                    segs.append(seg)
            else:
                ax.plot(d,
                        e[:, bandRange[0]:bandRange[1]],
                        c,
                        lw=lw,
                        alpha=alpha,
                        picker=2)

        if projections:
            species = [str(i) for i in bands.structure.species]
            data = bands.projections[s][:, i[0]:i[1] + 1]
            for el in elements:
                for orb in orbitals[el].keys():
                    colors = contrib(species, data, el, orb)
                    ax.add_collection(
                        LineCollection(np.array(segs).reshape(-1, 2, 2),
                                       colors=colors,
                                       lw=lw,
                                       alpha=alpha,
                                       linestyles=ls))

        if showVBM:
            ax.plot([d[0], d[-1]], [0, 0], color='k', ls='-', lw=1.5)
        if showCBM:
            ax.plot([d[0], d[-1]], [Eg, Eg], color='k', ls='-', lw=1.5)

    ###########################################################################

    bandRange = kargs['bandRange'] \
        if 'bandRange' in kargs else [0, bands.nb_bands]
    xlim = kargs['xlim'] if 'xlim' in kargs else (0, len(bands.branches))
    ylim = kargs['ylim'] if 'ylim' in kargs else [-4, 4]
    spin = kargs['spin'] if 'spin' in kargs else [Spin.up, Spin.down]
    zero_to_fermi = kargs['zero_to_fermi'] \
        if 'zero_to_fermi' in kargs else None
    showVBM = kargs['showVBM'] if 'showVBM' in kargs else None
    showCBM = kargs['showCBM'] if 'showCBM' in kargs else None
    bandMarks = kargs['bandMarks'] if 'bandMarks' in kargs else None
    projections = kargs['projections'] if 'projections' in kargs else None
    lw = kargs['lw'] if 'lw' in kargs else 1
    alpha = kargs['alpha'] if 'alpha' in kargs else 1
    if 'smooth' in kargs:
        smooth = kargs['smooth'][0]
        tol = kargs['smooth'][1]
    else:
        smooth = None
        tol = None

    ###########################################################################

    def onpick(event):
        band = event.artist.get_label()
        if event.mouseevent.button == 1:
            bandColor(ax, band, clickColor, clickWidth, alpha=clickAlpha)
        elif event.mouseevent.button == 3:
            bandColor(ax, band, 'k', lw, alpha)

        print("\n -> ", clickColor, ' @ band ', band, sep='')

        plt.draw()

    fig = ax.get_figure()
    fig.canvas.mpl_connect('pick_event', onpick)

    ###########################################################################

    def bandColor(ax, band, color='k', width=lw, alpha=alpha):
        """ Change color of a band by band number"""

        if len(list(color)) > 1:
            color = [i / 255 for i in color]

        for l in ax.lines:
            if l.get_label() == str(band):
                l.set_color(color)
                l.set_linewidth(width)
                l.set_alpha(alpha)

    ###########################################################################

    Eg = bands.get_band_gap()['energy']
    try:
        VBM = bands.get_vbm()['energy']
        CBM = bands.get_cbm()['energy']
        print('\nBand Gap = %0.3f eV\nVBM = %0.3f eV\nCBM = %0.3f eV' %
              (Eg, VBM, CBM))
    except:
        print('\nNo gap - metallic?')

    ###########################################################################

    if no_shift:
        zero_energy = 0.
    elif bands.is_metal():
        zero_energy = bands.efermi
    elif zero_to_fermi:
        zero_energy = bands.efermi
    else:
        zero_energy = bands.get_vbm()['energy']

    elements = bands.structure.symbol_set

    ###########################################################################

    k0 = bands.kpoints[0].label
    prevK = ('', 0)
    dk = 0
    kLabels = {'label': [], 'distance': []}
    totD = list()
    for branch in bands.branches[xlim[0]:(xlim[1] + 1)]:

        i = (branch['start_index'], branch['end_index'])

        # calculate x-axis data
        distances = []
        for k, d in zip(bands.kpoints[i[0]:i[1] + 1],
                        bands.distance[i[0]:i[1] + 1]):
            dk = setTicks(ax, dk)
            distances.append(d + dk)

        # calculate y-axis data
        e = bandrun.eigenvalues
        energies = {'1': None, '-1': None}
        for s in spin:
            energies[str(s)] = e[s][i[0]:i[1] + 1, :, 0]
            plotBranch(ax, np.array(distances), energies[str(s)] - zero_energy,
                       tol, i, 'k' if s.value == 1 else 'r',
                       '-' if s.value == 1 else '--')

        prevK = (k.label, d + dk)
        totD.append(distances)

    if bandMarks:
        shift = bands.branches[xlim[0]]['start_index']
        totD = np.array(totD).flatten()
        cMarks = 'k' if projections else bandMarkColor

        # automatic
        for i in [bands.get_vbm(), bands.get_cbm()]:
            ax.scatter(totD[i['kpoint_index'][0] - shift],
                       i['energy'] - zero_energy,
                       color=cMarks,
                       s=bandMarkSize)

        # manual
        # ax.scatter(distances[-1], bands.bands[Spin.up][156][-1] - zero_energy,
        #            color=cMarks, s=bandMarkSize)
        # ax.scatter(distances[-1], bands.bands[Spin.up][155][-1] - zero_energy,
        #            color=cMarks, s=bandMarkSize)

    ax.set_xticks(kLabels['distance'])
    ax.set_xticklabels(kLabels['label'], **ticks)
    ax.set_xlim(kLabels['distance'][0], kLabels['distance'][-1])

    if bandRange[0] != 0:
        yMin = 0
        yMax = 0
        for s in spin:
            eMin = (energies[str(s)][:, bandRange[0]] - zero_energy).min()
            eMax = (energies[str(s)][:, bandRange[1]] - zero_energy).max()
            if yMin > eMin:
                yMin = eMin
            if yMax < eMax:
                yMax = eMax
        # ax.set_ylim(np.floor(yMin), np.ceil(yMax))
        ax.set_ylim(ylim)
    else:
        ax.set_ylim(ylim)


def plotDOS(fig, dosrun, orient='y', **kargs):

    from scipy.interpolate import interp1d

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

        # global tempX, tempY
        # tempX = x
        # tempY = y

        if ax.get_label() == 'zoom':
            lw *= .75

        face = c

        if t == 'Total':
            face = c + 245

        if (shade or t == 'Tota') and ax.get_label() != 'zoom':
            var = 1
            zo = 1
            if 'Ti' in t:
                var = 1
            if t == 'Ti 2s':
                var = 4
                zo = 2
            ax.fill_between(x,
                            y * var,
                            edgecolor=c / 255,
                            facecolor=face / 255,
                            alpha=1.,
                            lw=lw,
                            ls=ls,
                            zorder=zo)
            handles.append(ax.collections[-1])
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

            if t == 'Total':
                handles.append(
                    Rectangle((0.0, 0.0),
                              .35,
                              .15,
                              ls=ls,
                              lw=0.5,
                              alpha=1.,
                              facecolor=face / 255,
                              edgecolor=c / 255))
            else:
                handles.append(
                    Rectangle((0.0, 0.0),
                              .35,
                              .15,
                              ls=ls,
                              lw=lw,
                              facecolor=c / 255,
                              edgecolor=c / 255))
            # ax.fill_between(x[i:j][eval(o)], y[i:j][eval(o)],
            #                 edgecolor=c * 0, facecolor=c / 255,
            #                 alpha=a, lw=0)

        ax.tick_params(
            axis='both',
            which='both',
            #    pad=ticks['pad'],
            labelsize=ticks['fontsize'])

        labels.append(t)

    ###########################################################################

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
    ldos = kargs['ldos'] if 'ldos' in kargs else False
    shift = kargs['shift'] if 'shift' in kargs else 0.
    plotTotal = kargs['plotTotal'] if 'plotTotal' in kargs else False
    plotOrbs = kargs['plotOrbs'] if 'plotOrbs' in kargs else False

    elements = list(dict.fromkeys(dosrun.atomic_symbols))
    labels = list()
    handles = list()
    lw = kargs['lw'] if 'lw' in kargs else 1

    e = dosrun.tdos.energies - dosrun.tdos.get_cbm_vbm()[1] - shift
    y = e
    if no_shift:
        e = dosrun.tdos.energies
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
            for ax in _axes:
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

    if legendOn:
        from collections import OrderedDict
        legend = OrderedDict(zip(labels, handles))

        _axes[0].legend(
            legend.values(),
            legend.keys(),
            loc=legendOn[2],
            # bbox_to_anchor=(0.5, 0.32, 0.5, 0.0),
            fontsize=legendOn[1])


def plotBZ():

    from pymatgen.electronic_structure.plotter import \
        plot_brillouin_zone_from_kpath as pltBZ

    st = bands.structure
    st_sym = sga(st)

    global prim, path, newPath

    prim = sga(st).get_primitive_standard_structure(
        international_monoclinic=False)

    t = '{0} - {1} ({2})'.format(st_sym.get_lattice_type().capitalize(),
                                 st_sym.get_space_group_symbol(),
                                 st_sym.get_space_group_number())
    path = HighSymmKpath(prim)
    pltBZ(path, title=t)


# %% MAIN #######################################

orbitals = {
    'Al': {
        's': ['3s', [129, 178, 214]],
        'p': ['3p', [64, 128, 128]],
    },
    'Ti': {
        's': ['2s', [180, 200, 80]],
        'd': ['3d', [200, 200, 200]]
    },
    'Cu': {
        # 's': ['4s', [190, 190, 0]],
        'd': ['3d', [0, 0, 150]]
    },
    'Zn': {
        's': ['2s', [255, 215, 0]],
        'd': ['3d', [150, 150, 150]]
    },
    'Bi': {
        's': ['6s', [240, 0, 140]],
        # 'p': ['6p', [0, 0, 200]]
    },
    'In': {
        's': ['5s', [5, 67, 235]],
        'd': ['4d', [0, 162, 232]]
    },
    'W': {
        # 's': ['5d', [255, 0, 0]],
        # 'p': ['5d', [0, 255, 0]],
        'd': ['5d', [120, 120, 120]]
        # 'd': ['5d', [10, 125, 0]]
    },
    'V': {
        'd': ['3d', [215, 165, 0]]
    },
    'Sb': {
        's': ['5s', [140, 180, 80]],
        'p': ['5p', [250, 150, 20]]
    },
    'Ta': {
        # 'd': ['5d', [183, 154, 86]]
        'd': ['5d', [0, 0, 150]]
    },
    'Nb': {
        'd': ['4d', [20, 160, 160]]
    },
    'O': {
        # 's': ['2s', [100, 0, 0]],
        'p': ['2p', [150, 0, 0]]
    },
    'S': {
        # 's': ['3s', [200, 0, 0]],
        # 'p': ['3p', [5, 67, 235]]
        'p': ['3p', [240, 185, 70]]
    },
    'Te': {
        's': ['5s', [0, 200, 0]],
        'p': ['5p', [0, 0, 200]]
    },
    'N': {
        'p': ['3p', 'white']
    },
    'Hg': {
        's': ['6s', [240, 0, 140]],
        'd': ['5d', [75, 0, 130]]
    },
    'Y': {
        # 's': ['4_5s', [240, 0, 140]],
        # 'p': ['4p', [75, 0, 130]],
        'd': ['4d', [120, 120, 120]]
    },
}

title = {'fontsize': 20, 'pad': 20}
axes = {'fontsize': 24, 'labelpad': 8}
ticks = {'fontsize': 18}  # , 'pad': 10
arrows = {'fontsize': 12, 'labelpad': 20}
font = {'family': 'Arial'}
opt = {'font': font}
for i in opt:
    plt.rc(i, **opt[i])

dosopt = {
    'xlim': (-10, 6),
    'ylim': (0, 75),
    'DOSticks': True,
    'shade': True,
    'spin': [Spin.up],  # comment out for both spins
    'zoom': [(-2.5, -1.5), (0, 0.5)],  # [0.55, 0.35, 0.2, 0.6]],
    'lw': 1.,
    'smooth': True,
    'shift': 0.04,  # positive shifts down
    'vertical': True,
    'plotTotal': True,
    'plotOrbs': True,
    'ldos': [0.18148, 0.44062, False],  # bulk true or false
    # 'ldos': [.61, .8, True],  # bulk true or false
    'legend': [True, 14, 4]
}

bandopt = {
    'ylim': (-11, 7),
    'xlim': (0, 0),  # numbers represent branches (0 is first)
    'spin': [Spin.up],  # comment out for both spins
    'showVBM': True,
    'showCBM': True,
    'zero_to_fermi': True,
    'bandRange': [569, 573],
    'bandMarks': True,
    'projections': True,
    'lw': 1.,
    'alpha': 0.25,
    'smooth': [True, 100]
}

clickColor = 'k'
clickWidth = 2
clickAlpha = 1.0
segWidth = 0.2
bandMarkSize = 25
bandMarkColor = 'k'

colorChoice = [False, 2]  # linewidth
no_shift = False
plt.close('all')

plots = [
    # 'bands',
    'dos',
    # 'BZ'
]

if 'zoom' in dosopt:
    if 'bands' in plots or 'vertical' in dosopt:
        dosopt['zoom'].reverse()

locs = []

plot(locs, parse=True, save=True, figsize=(4, 8), dpi=600, fmt='tiff', lw=1)

# plt.subplots_adjust(
#     top=0.9,
#     bottom=0.2,
#     left=0.22,
#     right=0.96,
# )

plt.tight_layout()
