import re
from typing import Dict, List, Tuple
from matplotlib.backend_bases import Event
import numpy as np

from sympy import symbols, Eq, Lt, lambdify
from sympy.core.symbol import Symbol
from sympy.solvers import solve

from shapely.geometry import LineString, Polygon, Point
from descartes import PolygonPatch

import matplotlib.pyplot as plt
from matplotlib import cm

from widgets import Slider, PremiumCheckButtons

###############################################################################
# CLASSES
###############################################################################


class Fragment:
    def __init__(self, name: str, total_energy: float,
                 color: Tuple[float]) -> None:
        self.name = name
        self.components = self.get_components()
        self.E = total_energy
        self.H = None
        self.color = color

    def get_components(self) -> Dict[str, int]:
        """
        Get fragment componets as {element: fraction} dictionary.

        Example: (CuBiW2O8)

        components = [Cu, Bi, W2, O8]
        parts = [Cu], [Bi], [W, 2], [O, 8]
        compound_parts = {'Cu': 1, 'Bi': 1, 'W': 2, 'O': 8}
        """

        compound_parts = {}
        components = re.findall(r'[A-Z][a-z]?[0-9]?', self.name)

        for component in components:
            parts = re.split(r'(\d+)', component)
            element, n = parts[0], 1
            if len(parts) > 1:
                n = int(parts[1])
            compound_parts[element] = n

        return compound_parts

    def compute_enthalpy(self, bulk_energies: Dict[str, float]) -> None:
        """Compute the enthalpy of the compound."""

        enthalpy = self.E

        for element, n in self.components.items():
            enthalpy -= n * bulk_energies[element]

        self.H = round(enthalpy, 3)

    def generate_ineq(self,
                      axes: Tuple[Symbol],
                      primary: bool = False) -> None:

        x, y, z = axes
        expr = x - x

        for element, fraction in self.components.items():
            expr += fraction * symbols(element)

        if primary:
            Fragment.prim_eq = Eq(expr, self.H)
            self.ineq = solve(Lt(expr, self.H), y).subs(z, 0)
        else:
            ineq = Lt(expr, self.H).subs(z, solve(Fragment.prim_eq, z)[0])

            if ineq.has(y):
                self.ineq = solve(ineq, y)
            else:
                self.ineq = solve(ineq, x)


class CPL:
    def __init__(self,
                 data_location: str,
                 cpl_axes: Tuple[Symbol],
                 line_opacity: float = 0.75,
                 shade_opacity: float = 0.2,
                 free: Symbol = None,
                 SP: Polygon = None,
                 colormap: str = 'Dark2') -> None:

        # define CPL axes
        self.axes = cpl_axes
        self.free = free

        # define color map
        self.colormap = colormap

        # load fragments
        self.frags = self.load_fragments(data_location)
        self.primary = self.frags[0]

        # TODO: Implement defect analysis in OOP

        # compute CPL lines
        self.generate_inequalities()
        self.shades = []

        # plotting parameters
        self.line_opacity = line_opacity
        self.shade_opacity = shade_opacity
        self.SP = SP
        self.check = None
        self.s_free = None

        self.initiate_plot()  # basic plot setup

    def load_fragments(self, loc: str) -> List[Fragment]:

        bulk_energies = self.get_bulk_energies(loc)  # {element: energy}
        total_energies = self.get_total_energies(loc)  # {compound: energy}

        colors = cm.get_cmap(self.colormap, len(total_energies) - 1)

        fragments = []

        for i, (name, energy) in enumerate(total_energies.items()):

            # set color
            if i == 0:
                color = (0., 0., 0., 1.)  # black
            else:
                color = colors(i - 1)

            fragment = Fragment(name, energy, color)
            fragment.compute_enthalpy(bulk_energies)
            fragments.append(fragment)

        return fragments

    def get_bulk_energies(self, loc: str) -> Dict[str, float]:
        """Get bulk energies from file."""

        bulk_energies = {}

        try:
            with open(loc + 'bulkE.txt') as bulks_file:
                for line in bulks_file.readlines():
                    element, energy = line.split()
                    if element == 'O':
                        energy = float(energy) / 2
                    else:
                        energy = float(energy)
                    bulk_energies[element] = energy
        except Exception as err:
            print('\nOperation failed: %s' % err.strerror)
            return {}
        else:
            return bulk_energies

    def get_total_energies(self, loc: str) -> Dict[str, float]:
        """Get total energies from file."""

        total_energies = {}

        try:
            with open(loc + 'totalE.txt') as totals_file:

                for line in totals_file.readlines():
                    compound, energy = line.split()
                    total_energies[compound] = float(energy)

        except Exception as err:
            print('\nOperation failed: %s' % err.strerror)
            return {}
        else:
            return total_energies

    def generate_inequalities(self) -> None:
        """Compute the line equation of the compound w.r.t the CPL parameters."""

        self.primary.generate_ineq(self.axes, True)

        for frag in self.frags[1:]:
            frag.generate_ineq(self.axes)

    def initiate_plot(self) -> None:

        self.fig, self.ax = plt.subplots(figsize=(10, 8))
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.85, top=0.85)
        self.ax.xaxis.tick_top()
        self.ax.xaxis.set_label_position('top')
        self.ax.yaxis.tick_right()
        self.ax.yaxis.set_label_position('right')
        self.ax.spines['left'].set_visible(False)
        self.ax.spines['bottom'].set_visible(False)
        self.ax.set_xlabel(r'$\Delta\mu_{%s}$' % str(self.axes[0]), **axes)
        self.ax.set_ylabel(r'$\Delta\mu_{%s}$' % str(self.axes[1]), **axes)
        plt.minorticks_on()

    def load_widgets(self) -> None:

        if len(self.primary.components) > 3:

            # [left, bottom, width, height]
            ax_slider = plt.axes([0.1, 0.05, 0.75, 0.02])

            self.s_free = Slider(ax_slider,
                                 r'$\Delta\mu_{%s}$' % str(self.free),
                                 -2.0,
                                 0.0,
                                 valinit=0.0,
                                 valstep=0.001,
                                 valfmt="%1.3f")

            self.s_free.on_changed(self.update)
            self.s_free.label.set_size(20)

        ax_checks = plt.axes([0.1, 0.1, 0.3, 0.4], frameon=False)
        visibility = [
            shade.set_visible(False) if shade else [] for shade in self.shades
        ]
        self.check = PremiumCheckButtons(ax_checks,
                                         self.frags[1:],
                                         visibility,
                                         loc=3,
                                         borderaxespad=0)

        self.check.on_clicked(self.apply_shade)

    def draw(self, free_val: float) -> None:

        x, y = self.axes[:2]  # assign axes

        lines = []
        shades = []
        poly_x, poly_y, xy, v = np.zeros(4)

        for i, frag in enumerate(self.frags):

            if i == 0:  # primary phase

                prim = frag.ineq.rhs.subs(self.free, free_val)  # primary line

                # set CPL axes limits
                x_min = round(float(solve(prim, x)[0]), 3)
                y_min = round(float(prim.subs(x, 0)), 3)
                self.ax.set_xlim(x_min, 0.0)
                self.ax.set_ylim(y_min, 0.0)

                f = lambdify(x, prim)  # line function
                t = np.array([x_min, 0])  # x-axis end points

                # plot line and add to lines register
                lines.append(
                    self.ax.plot(t,
                                 f(t),
                                 color='k',
                                 zorder=len(self.frags),
                                 label=frag.name))

                # define thermodynamically stable region for primary
                prim_poly = Polygon([[x_min, 0], [0, 0], [0, y_min]])

            else:  # secondary phases

                sec = frag.ineq.rhs.subs(self.free, free_val)  # secondary line

                if not frag.ineq.has(y):  # vertical

                    A = (sec, 1)
                    B = (sec, y_min - 1)
                    seg = LineString([A, B])
                    pts = seg.intersection(prim_poly.boundary)

                    if isinstance(pts, LineString):
                        t = [x_min, x_min]
                        v = [0, 0]
                    else:
                        t = np.array(
                            [p[0][0] for p in [p.coords.xy for p in pts]])
                        v = np.array(
                            [p[1][0] for p in [p.coords.xy for p in pts]])

                    if frag.ineq.rel_op == '>':
                        xVal = 1
                    else:
                        xVal = x_min - 1

                    poly_x = np.array([A[0], xVal, xVal, B[0]])
                    poly_y = np.array([A[1], A[1], B[1], B[1]])

                else:  # non-vertical

                    f = lambdify(x, sec)
                    A = (x_min - 1000, f(x_min - 1000))
                    B = (1000, f(1000))
                    seg = LineString([A, B])
                    pts = seg.intersection(prim_poly.boundary)

                    if isinstance(pts, LineString):
                        t = [x_min, 0]
                        v = [0, y_min]
                    else:
                        t = np.array(
                            [p[0][0] for p in [p.coords.xy for p in pts]])
                        v = np.array(
                            [p[1][0] for p in [p.coords.xy for p in pts]])

                    if not frag.ineq.has(x):  # horizontal

                        if frag.ineq.rel_op == '>':
                            poly_x = np.array([A[0], A[0], B[0], B[0]])
                            poly_y = np.array([A[1], 1, 1, B[1]])
                        else:
                            poly_x = np.array([A[0], A[0], B[0], B[0]])
                            poly_y = np.array(
                                [A[1], y_min - 1, y_min - 1, B[1]])

                    else:  # slopped

                        poly_x = [A[0], B[0]]
                        poly_y = [A[1], B[1]]

                        if f(0) < f(1):  # positive slope

                            if frag.ineq.rel_op == '>':
                                poly_x.append(A[0])
                                poly_y.append(B[1])
                            else:
                                poly_x.append(B[0])
                                poly_y.append(A[1])

                        else:  # negative slope

                            if frag.ineq.rel_op == '>':
                                poly_x.append(B[0])
                                poly_y.append(A[1])
                            else:
                                poly_x.append(A[0])
                                poly_y.append(B[1])

                lines.append(
                    self.ax.plot(t,
                                 v,
                                 color=frag.color,
                                 alpha=self.line_opacity,
                                 label=frag.name,
                                 lw=2))

                poly = Polygon(
                    np.array([[i, j] for i, j in zip(poly_x, poly_y)]))

                overlap = poly.intersection(prim_poly)
                if overlap:
                    shades.append(
                        PolygonPatch(overlap,
                                     color=frag.color,
                                     alpha=self.shade_opacity,
                                     label=frag.name))
                else:
                    shades.append([])

        return shades

    def single_phase(self) -> None:
        visibility = self.check.get_status()
        labels = [patch.get_label() for patch in self.ax.patches]
        if 'SP' in labels:
            self.ax.patches.remove(self.ax.patches[labels.index('SP')])
        if visibility.count(True) > 1:
            if len(self.ax.patches) > 1 and len(
                    self.ax.patches) == visibility.count(True):
                self.SP = Polygon(self.ax.patches[0].get_path().vertices)
                for act in self.ax.patches:
                    self.SP = self.SP.intersection(
                        Polygon(act.get_path().vertices))
                if self.SP and self.SP.geom_type != 'LineString':
                    if self.SP.geom_type != 'Polygon':
                        self.SP = self.SP.geoms[1]
                    self.ax.add_patch(
                        PolygonPatch(self.SP,
                                     color='k',
                                     alpha=0.75,
                                     label='SP')).set_hatch('x')

    def update(self, val: float) -> None:
        plt.close(2)
        visibility = self.check.get_status()
        self.ax.lines.clear()
        self.shades.clear()
        self.ax.patches.clear()
        for vis, shade in zip(visibility, self.draw(self.s_free.val)):
            self.shades.append(shade)
            if shade:
                if vis:
                    self.ax.add_patch(shade)
                    shade.set_visible(True)
                else:
                    shade.set_visible(False)
        self.single_phase()

    def apply_shade(self, label: str) -> None:
        plt.close(2)
        i = [l.get_label() for l in self.ax.lines[1:]
             ].index(re.sub(r'\$_\{([0-9]*)\}\$', r'\1', label))
        if self.shades[i]:
            self.shades[i].set_visible(not self.shades[i].get_visible())
            if not self.shades[i] in self.ax.patches:
                self.ax.add_patch(self.shades[i])
            else:
                self.ax.patches.remove(self.shades[i])
        self.single_phase()
        plt.draw()


###############################################################################
# PLOTTING PARAMETERS
###############################################################################

font = {'family': 'times new roman', 'size': 16}
axes = {'fontsize': 20, 'labelpad': 20}
ticks = {'labelsize': 16}
math = {'rm': 'Roman'}
params = {'font': font, 'xtick': ticks, 'ytick': ticks, 'mathtext': math}
for p in params:
    plt.rc(p, **params[p])

###############################################################################
# MAIN PROGRAM
###############################################################################

if __name__ == '__main__':

    # define CPL axes
    free, x, y, z = symbols(['Cu', 'Bi', 'W', 'O'])

    loc = 'tests/CBTO/'  # location of energies

    cpl_setup = {
        'data_location': loc,
        'cpl_axes': (x, y, z),
        'free': free,
        'line_opacity': 0.75,
        'shade_opacity': 0.2,
        'SP': None
    }

    cpl = CPL(**cpl_setup)
    cpl.shades = cpl.draw(0.0)
    cpl.load_widgets()

    plt.show()
