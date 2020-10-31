"""
TO-DO:
    - Remove free variable requirement
    - Allow for fractional atomic contribution
    - Fix checkbox layout
"""

import re
import numpy as np

from sympy import symbols, Lt, lambdify
from sympy.solvers import solve

from shapely.geometry import LineString, Polygon, Point
from descartes import PolygonPatch

import matplotlib.pyplot as plt

from widgets import Slider, PremiumCheckButtons

loc = 'CPL Energies/'  # location of energies

Cu, Bi, W, O = symbols(['Cu', 'Bi', 'W', 'O'])
free, x, y, z = Cu, Bi, W, O

# Globals
SP = None


def getEnthalpy(_dict, f, m=None):

    # collect number of each species
    frac = dict()
    if not m:
        for b in bulks:
            if b in f:
                n = re.search(re.escape(b) + '[0-9]+', f)
                if n:
                    frac[b] = int(
                        re.search(r"[0-9][0-9]*", n.group(0)).group(0))
                else:
                    frac[b] = 1
    elif re.match('Vac', f):
        frac.update((k, v - 1) if k == f.split('_')[1] else (k, v)
                    for k, v in Super['fracs'].items())
    else:
        s = f.split('_')
        frac.update((k, v + 1) if k == s[0] else (k, v)
                    for k, v in Super['fracs'].items())
        frac.update(
            (k, v - 1) if k == s[1] else (k, v) for k, v in frac.items())
    _dict[f]['fracs'] = frac
    # calculate enthalpy
    return _dict[f]['E'] - sum(_dict[f]['fracs'][i] * bulks[i]
                               for i in _dict[f]['fracs'])


# Collect DFT energies and calculate enthalpies
frags = dict()
bulks = dict()
try:
    with open(loc + 'bulkE.txt', 'r') as b, open(loc + 'totalE.txt', 'r') as t:
        for line in b.readlines():
            i = line.split()
            if i[0] == 'O':
                bulks[i[0]] = float(i[1]) / 2
            else:
                bulks[i[0]] = float(i[1])
        for line in t.readlines():
            i = line.split()
            frags[i[0]] = {}
            frags[i[0]]['E'] = float(i[1])
            frags[i[0]]['H'] = round(getEnthalpy(frags, i[0]), 3)
except IOError as e:
    print('\nOperation failed: %s' % e.strerror)

prim = next(iter(frags.keys()))

defects = dict()
try:
    with open(loc + 'defectE.txt', 'r') as d:
        for line in d.readlines():
            i = line.split()
            m = 1
            if i[0] == 'Super':
                Super = dict()
                Super['E'] = float(i[1])
                m = round(Super['E'] / frags[prim]['E'])
                frac = dict()
                frac.update(
                    (k, m * v) for k, v in frags[prim]['fracs'].items())
                Super['fracs'] = frac
                Super['H'] = round(
                    Super['E'] - sum(Super['fracs'][i] * bulks[i]
                                     for i in Super['fracs']), 3)
            else:
                defects[i[0]] = dict()
                defects[i[0]]['E'] = float(i[1])
                defects[i[0]]['H'] = round(
                    getEnthalpy(defects, i[0], m) / m, 3)
                defects[i[0]]['diff'] = {
                    k: v
                    for k, v in zip(bulks, [
                        s - d for s, d in zip(Super['fracs'].values(), defects[
                            i[0]]['fracs'].values())
                    ]) if k in i[0].split('_')
                }
except IOError as e:
    print('\nOperation failed: %s' % e.strerror)

# %% Define equations, inequalities, lines, and shading

for j, f in enumerate(frags):

    # EQUATIONS | INEQUALITIES #
    s = x
    s -= s  # define empty variable
    for l, v in frags[f]['fracs'].items():
        s += v * symbols(l)
    frags[f]['eq'] = s - frags[f]['H']  # generate equations
    if f == prim:  # define primary line
        frags[f]['line'] = solve(frags[f]['eq'], y)[0].subs(z, 0)
        p = solve(frags[f]['eq'], z)[0]  # solve in terms of z variable
    else:
        frags[f]['ineq'] = Lt(s, frags[f]['H'])  # generate inequalities
        sol = frags[f]['eq'].subs(z, p)  # solve by z back-substitution

        # LINES | SHADES #
        if str(y) in str(sol):  # non-vertical lines
            frags[f]['line'] = solve(sol, y)[0]
            frags[f]['shade'] = solve(
                frags[f]['ineq'].subs(z,
                                      solve(frags[prim]['eq'], z)[0]), y)
        else:  # vertical lines
            frags[f]['line'] = solve(sol, x)[0]
            frags[f]['shade'] = solve(
                frags[f]['ineq'].subs(z,
                                      solve(frags[prim]['eq'], z)[0]), x)

# %% Define CPL

fig, ax = plt.subplots(figsize=(10, 8))
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.85, top=0.85)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.yaxis.tick_right()
ax.yaxis.set_label_position('right')
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
plt.minorticks_on()

font = {'family': 'times new roman', 'size': 16}
axes = {'fontsize': 20, 'labelpad': 20}
ticks = {'labelsize': 16}
math = {'rm': 'Roman'}
opt = {'font': font, 'xtick': ticks, 'ytick': ticks, 'mathtext': math}
for i in opt:
    plt.rc(i, **opt[i])

ax.tick_params(axis='both', **ticks)

a = .75  # line opacity
a_s = 0.2  # shade opacity

colors = plt.cm.get_cmap('Dark2', len(frags) - 1)


def draw_CPL(c):
    ax.set_xlabel(r'$\Delta\mu_{%s}$' % str(x), **axes)
    ax.set_ylabel(r'$\Delta\mu_{%s}$' % str(y), **axes)
    lines = []
    shades = []
    poly_x, poly_y, xy, v = np.zeros(4)
    for i, f in enumerate(frags):
        if i == 0:  # primary line
            q = frags[f]['line'].subs(free, c)
            minX = round(float(solve(q, x)[0]), 3)
            minY = round(float(q.subs(x, 0)), 3)
            ax.axis([minX, 0.0, minY, 0.0])
            line = lambdify(x, q)
            t = np.array([minX, 0])
            lines.append(
                ax.plot(t, line(t), color='k', zorder=len(frags), label=f))
            primPoly = Polygon([[minX, 0], [0, 0], [0, minY]])
        else:
            l = frags[f]['line'].subs(free, c)
            s = str(frags[f]['shade'])
            if not str(y) in s:  # vertical
                A = (l, 1)
                B = (l, minY - 1)
                seg = LineString([A, B])
                pts = seg.intersection(primPoly.boundary)
                t = np.array([p[0][0] for p in [p.coords.xy for p in pts]])
                v = np.array([p[1][0] for p in [p.coords.xy for p in pts]])
                if str(">") in s:
                    xVal = 1
                else:
                    xVal = minX - 1
                poly_x = np.array([A[0], xVal, xVal, B[0]])
                poly_y = np.array([A[1], A[1], B[1], B[1]])
            else:  # non-vertical
                line = lambdify(x, l)
                A = (minX - 1000, line(minX - 1000))
                B = (1000, line(1000))
                seg = LineString([A, B])
                pts = seg.intersection(primPoly.boundary)
                t = np.array([p[0][0] for p in [p.coords.xy for p in pts]])
                v = np.array([p[1][0] for p in [p.coords.xy for p in pts]])
                if not str(x) in s:  # horizontal
                    if str(">") in s:
                        poly_x = np.array([A[0], A[0], B[0], B[0]])
                        poly_y = np.array([A[1], 1, 1, B[1]])
                    else:
                        poly_x = np.array([A[0], A[0], B[0], B[0]])
                        poly_y = np.array([A[1], minY - 1, minY - 1, B[1]])
                else:  # slopped
                    poly_x = [A[0], B[0]]
                    poly_y = [A[1], B[1]]
                    if line(0) < line(1):  # positive slope
                        if str(">") in s:
                            poly_x.append(A[0])
                            poly_y.append(B[1])
                        else:
                            poly_x.append(B[0])
                            poly_y.append(A[1])
                    else:  # negative slope
                        if str(">") in s:
                            poly_x.append(B[0])
                            poly_y.append(A[1])
                        else:
                            poly_x.append(A[0])
                            poly_y.append(B[1])

            lines.append(
                ax.plot(t, v, color=colors(i - 1), alpha=a, label=f, lw=2))
            poly = Polygon(np.array([[i, j] for i, j in zip(poly_x, poly_y)]))
            overlap = poly.intersection(primPoly)
            if overlap:
                shades.append(
                    PolygonPatch(overlap,
                                 color=colors(i - 1),
                                 alpha=a_s,
                                 label=f))
            else:
                shades.append([])

    # mark single phase region
    # global SP
    # SP = Polygon(shades[0].get_path().vertices)
    # for act in shades[1:2]:
    #     SP = SP.intersection(Polygon(act.get_path().vertices))
    #     if SP and SP.geom_type != 'LineString':
    #         if SP.geom_type != 'Polygon': SP = SP.geoms[1]
    #         ax.add_patch(PolygonPatch(SP,color='k',alpha=0.75,label='SP')).set_hatch('x')

    return shades


shades = draw_CPL(0.0)

# %% Define widgets


def singlePhase():
    visibility = check.get_status()
    labels = [patch.get_label() for patch in ax.patches]
    global SP
    if 'SP' in labels: ax.patches.remove(ax.patches[labels.index('SP')])
    if visibility.count(True) > 1:
        if len(ax.patches) > 1 and len(ax.patches) == visibility.count(True):
            SP = Polygon(ax.patches[0].get_path().vertices)
            for act in ax.patches:
                SP = SP.intersection(Polygon(act.get_path().vertices))
            if SP and SP.geom_type != 'LineString':
                if SP.geom_type != 'Polygon': SP = SP.geoms[1]
                ax.add_patch(
                    PolygonPatch(SP, color='k', alpha=0.75,
                                 label='SP')).set_hatch('x')


def update(val):
    plt.close(2)
    visibility = check.get_status()
    ax.lines.clear()
    shades.clear()
    ax.patches.clear()
    for vis, shade in zip(visibility, draw_CPL(sfree.val)):
        shades.append(shade)
        if shade:
            if vis:
                ax.add_patch(shade)
                shade.set_visible(True)
            else:
                shade.set_visible(False)
    singlePhase()


def shade(label):
    plt.close(2)
    i = [l.get_label() for l in ax.lines[1:]
         ].index(re.sub(r'\$_\{([0-9]*)\}\$', r'\1', label))
    if shades[i]:
        shades[i].set_visible(not shades[i].get_visible())
        if not shades[i] in ax.patches:
            ax.add_patch(shades[i])
        else:
            ax.patches.remove(shades[i])
    singlePhase()
    plt.draw()


def genDefect(event):
    plt.close(2)
    if event.xdata and 'SP' in [patch.get_label()
                                for patch in ax.patches] and SP.contains(
                                    Point(event.xdata, event.ydata)):

        def delta(i, k):
            if symbols(k) == x:
                val = pts[:, 0][i]
            elif symbols(k) == y:
                val = pts[:, 1][i]
            elif symbols(k) == free:
                val = sfree.val
            else:
                val = solve(frags[prim]['eq'], z)[0].subs({
                    x: pts[:, 0][i],
                    y: pts[:, 1][i],
                    free: sfree.val
                })
            return val

        figDef = plt.figure(figsize=(12, 8))
        gs = figDef.add_gridspec(2, 3)
        axDef = plt.subplot(gs[:, :2])
        axZoom = plt.subplot(gs[0, 2])

        bounds = SP.boundary.coords.xy
        pts = np.array([[x, y] for x, y in zip(bounds[0], bounds[1])])
        xDat = np.arange(1, len(pts))
        yDat = np.zeros([len(xDat), len(defects)])
        for i in range(len(xDat)):
            for j, d in zip(range(len(defects)), defects):
                yDat[i, j] = defects[d]['E'] - Super['E'] + sum([
                    v * (bulks[k] + delta(i, k))
                    for k, v in defects[d]['diff'].items()
                ])
        legDef = []
        for d in defects:
            if re.match('Vac', d):
                label = re.sub('Vac_(.*)', r'V$_{\1}$', d)
            else:
                label = re.sub('_(.*)', r'$_{\1}$', d)
            legDef.append(label)

        axZoom.plot(pts[:, 0], pts[:, 1], 'k--')
        axZoom.set_axis_off()
        for i, p in zip(xDat, pts):
            xs = 0.01
            ys = 0.01
            axZoom.text(p[0] + xs, p[1] + ys, i)

        axDef.plot(xDat, yDat)
        axDef.set_xlim([xDat[0], xDat[-1]])
        axDef.set_ylim(-0.5, 3.0)
        axDef.set_xticks(xDat)
        axDef.set_xticklabels(xDat)
        for i in xDat:
            axDef.axvline(i, color='k', lw=0.5)
        axDef.set_ylabel(r'$\Delta H_f$ (eV)')
        axDef.axhline(y=0, ls='--', lw=0.5, color='k')
        figDef.legend(legDef,
                      ncol=1,
                      prop={'size': 14},
                      loc='center',
                      bbox_to_anchor=(0.78, 0.25))

        plt.show()


if len(frags[prim]['fracs'].keys()) > 3:
    axSlider = plt.axes([0.1, 0.05, 0.75, 0.02])  # [left,bottom,width,height]
    sfree = Slider(axSlider,
                   r'$\Delta\mu_{%s}$' % str(free),
                   round(float(solve(frags[prim]['line'].subs(x, 0))[0]), 0),
                   0.0,
                   valinit=0.0,
                   valstep=0.001,
                   valfmt="%1.3f")
    sfree.on_changed(update)
    sfree.label.set_size(20)

axChecks = plt.axes([0.1, 0.1, 0.3, 0.4], frameon=False)
visibility = [shade.set_visible(False) if shade else [] for shade in shades]
check = PremiumCheckButtons(axChecks,
                            ax.lines[1:],
                            visibility,
                            loc=3,
                            borderaxespad=0)

check.on_clicked(shade)

cid = fig.canvas.mpl_connect('button_press_event', genDefect)

plt.show()
