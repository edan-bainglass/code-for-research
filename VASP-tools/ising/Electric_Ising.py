import inspect
from typing import Iterable, List, Tuple

import numpy as np
from numpy import exp, dot, ndarray, sin, cos, arccos, arctan2
from numpy.random import random
from numpy.linalg import norm

from pymatgen.core.structure import Structure

import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# =============================================================================
#                                   Classes
# =============================================================================


class System:
    """A representation of the ferroelectric system."""
    def __init__(self, loc: str, supercell_dim: List[int] = []) -> None:
        self.struct: Structure = Structure.from_file(loc + '/POSCAR')

        self.add_charges(loc)  # add valence charges using Bader analysis data

        if supercell_dim:
            self.supercell(supercell_dim)

        self.no_anions = self.strip_anions('O')
        self.cations, self.anions = self.build_lattice()

        self.init_cations = self.cations.copy()
        self.init_polar = sum(self.init_cations[:, 1]) / self.get_volume()

        # average registers
        self.avg_cell = []
        self.avg_p = []

    def add_charges(self, loc: str) -> None:
        """
        Extract Bader charges from Bader analysis output.
        Bader analysis obtained via the Henkelman group code.

        valence charge = pseudopotential charge - Bader analysis charge
        """

        z = {'Sb': 5.0, 'Ta': 11.0, 'O': 6.0}  # pseudopotential charges
        bader_charges = []
        valence_charges = []
        with open(loc + '/ACF.dat') as bader_file:
            for line in bader_file:
                if line.lstrip()[0].isnumeric():
                    bader_charges.append(line.split()[4])

            for i in self.struct.symbol_set:
                for j in self.struct.indices_from_symbol(i):
                    valence_charges.append(z[i] - float(bader_charges[j]))

        self.struct.add_site_property('charge', valence_charges)

    def supercell(self, dim: List[int]) -> None:
        """Supercell transformation."""
        self.struct.make_supercell(dim)

    def strip_anions(self, anion: str) -> Structure:
        """Returns the system as a structure stripped of its anion species."""
        cations: Structure = self.struct.copy()
        cations.remove_species(anion)
        return cations

    def get_volume(self) -> float:
        """Returns the volume of the system."""
        return self.struct.volume

    def get_cations(self) -> ndarray:
        """Return the cationic sites of the system."""
        return self.cations

    def build_lattice(self) -> Tuple[Iterable[float]]:
        """Initializes the lattice with calculated dipole moments."""

        # calculate dipole moments
        def calc_dipoles():

            # calculate center of mass of cation polyhedron
            def center_of_mass(sites):
                center = np.zeros(3)
                total_weight = 0
                for site in sites:
                    wt = site.species.weight  # species_and_occu has been deprecated
                    center += site.coords * wt
                    total_weight += wt
                return center / total_weight

            dips = []
            c = self.struct.site_properties['charge']

            for atom in enumerate(self.no_anions):

                # define local molecule
                nn = self.struct.get_neighbors(atom[1],
                                               2.9,
                                               include_index=True)
                sites = [atom[1]] + [nn[i][0] for i in range(len(nn))]
                cm = center_of_mass(sites)

                # calculate dipole
                dip = np.zeros(3)
                for i in enumerate(sites):
                    if i[0] == 0:
                        # cation
                        dip += (i[1].coords - cm) * c[atom[0]]
                    else:
                        # oxygens
                        dip += (i[1].coords - cm) * c[nn[i[0] - 1][2]]
                dips.append(dip)

            for i in range(len(oxy)):
                dips.append([])  # add empty dipoles on oxygen sites

            return dips

        coords = list(self.struct.cart_coords)
        labels = [
            str(x[1]) + str(x[0] + 1) for x in enumerate(self.struct.species)
        ]
        oxy = self.struct.indices_from_symbol('O')
        dipoles = calc_dipoles()

        # combine coordinates, dipoles, and labels for each site
        lat = np.array([(c, d, l) for c, d, l in zip(coords, dipoles, labels)],
                       dtype=object)

        # split cations and anions
        cations = lat[:min(oxy)]
        anions = lat[min(oxy):]

        return (cations, anions)


class ElectricIsing:
    """An Ising model simulation object."""
    def __init__(self, system: System) -> None:

        self.system = system
        self.lat = system.cations

    def run_simulation(self, temperatures: List[int], eql_swps: int,
                       pol_swps: int, d_ang: int) -> None:
        """
        Executes Ising algorithm to obtain average structure polarization.
        """
        for t in temperatures:

            # thermalize lattice
            for _ in range(eql_swps):
                self.sweep(t, d_ang)

            # calculate polarization
            P = np.zeros(3)
            dips = 0.0
            for n in range(pol_swps):
                self.sweep(t, d_ang)
                P += sum(self.lat[:, 1])
                dips += self.lat[:, 1]

            # normalize
            self.system.avg_p.append(P / pol_swps / self.system.get_volume())
            dips = dips / pol_swps

            # save dipole moment configuration at current temperature
            cell = self.lat.copy()
            for i in range(len(self.lat)):
                cell[i][1] = dips[i]
            self.system.avg_cell.append(cell)

    def sweep(self, t: int, d_ang: float = 0.5) -> None:
        """
        Sweep through lattice, perturbing dipole moments and computing the
        energetic cost. If favorable, new dipole moment is retained. Total
        dipole-dipole interaction potential energy and polarizations are
        obtained in the process.

        Units:

        [r] = Å
        [p] = eÅ, or Debye
        [E] = eV

        E = 1/4πε0 * 1/rn^3 * [p1 • p2 - 3(p1 • rh)(p2 • rh)]
        """

        #######################################################################
        # Helper functions
        #######################################################################

        def perturb(p1: ndarray) -> ndarray:
            """Return the perturbed dipole moment."""
            p = p1.copy()
            theta = arccos(p[2] / norm(p)) + d_ang * (2 * random() - 1)
            phi = arctan2(p[1], p[0]) + d_ang * (2 * random() - 1)
            p[0] = norm(p1) * sin(theta) * cos(phi)
            p[1] = norm(p1) * sin(theta) * sin(phi)
            p[2] = norm(p1) * cos(theta)

            return p

        def calcE(p: ndarray) -> float:
            """Calculate dipole-dipole interaction potential."""
            E = 0.0
            for j in range(len(self.lat)):
                if j != i:
                    r2, p2 = self.lat[j][0], self.lat[j][1]
                    r = r2 - r1
                    rn = norm(r)
                    rh = r / rn
                    E += 14.4 / rn**3 * (dot(p, p2) -
                                         3 * dot(p, rh) * dot(p2, rh))
            return E

        #######################################################################

        # perform sweep
        for i in range(len(self.lat)):

            r1, p1 = self.lat[i][:2]
            e1 = calcE(p1)
            p1p = perturb(p1)
            e2 = calcE(p1p)
            dE = round(e2 - e1, 10)

            # if perturbation is energetically favorable
            # or thermally accessible...
            if dE < 0 or random() <= exp(-dE / (8.617e-5 * t)):
                self.lat[i][1] = p1p  # keep perturbation


class Plotter:
    """A plotting tool for an Electric Ising model."""
    def __init__(self, system: System) -> None:
        self.system: System = system

    def draw_lattice(self,
                     fig_size: Tuple[int] = (5, 5),
                     ax=None,
                     tilt=0,
                     turn=0,
                     oxy=False,
                     title=None,
                     labels=False) -> None:
        """Draws lattice with superimposed dipole moments."""

        A, B, C = self.system.struct.lattice.abc

        # if inspect.stack()[1][3] == '<module>': plt.close()

        if not ax:
            fig = plt.figure(figsize=fig_size)
            ax = fig.add_subplot(111, projection='3d')

        ax.set_axis_off()
        ax.set_xlim3d(0, A)
        ax.set_ylim3d(0, B)
        ax.set_zlim3d(0, C)

        # plot cell boundaries
        self.plot_cell(ax, [(0, 0, 0), (0, B, 0), (A, 0, 0), (0, 0, C)])

        # add oxygens if requested
        if oxy:
            lat = np.array(
                list(self.system.cations) + list(self.system.anions))
        else:
            lat = self.system.cations

        for i in enumerate(lat):
            if i[0] in self.system.struct.indices_from_symbol('Sb'):
                c, r = 'sandybrown', .4
            elif i[0] in self.system.struct.indices_from_symbol('Ta'):
                c, r = 'darkkhaki', .6
            else:
                c, r = 'darkred', .15
            x = i[1][0][0]
            y = i[1][0][1]
            z = i[1][0][2]

            self.draw_sphere(ax, x, y, z, c, r)  # draw atom
            if labels:
                ax.text(x + .25, y + .25, z + .25, i[1][2],
                        fontsize=10)  # add labels

            if i[0] not in self.system.struct.indices_from_symbol('O'):
                u = i[1][1][0]
                v = i[1][1][1]
                w = i[1][1][2]
                mag = norm(i[1][1])

                ax.quiver(x,
                          y,
                          z,
                          u,
                          v,
                          w,
                          length=mag * 0.025,
                          color='k',
                          arrow_length_ratio=0.1,
                          pivot='middle')  # draw dipole moment

        if self.system.init_polar.any():
            P = self.system.init_polar / norm(self.system.init_polar)

            # draw polarization vector
            ax.quiver(A / 2,
                      B / 2,
                      C / 2,
                      P[0],
                      P[1],
                      P[2],
                      length=4.,
                      color='g',
                      arrow_length_ratio=0.1,
                      pivot='middle')

        ax.set_aspect('auto')  # scale plot
        ax.view_init(tilt, turn)  # set viewing direction
        if title: ax.set_title(title, {'fontsize': 20})

    def plot_cell(self, ax, bounds) -> None:
        """Draw cell boundaries."""

        bounds_array = [np.array(list(i)) for i in bounds]

        points = []
        points += bounds_array
        vectors = [
            bounds_array[1] - bounds_array[0],
            bounds_array[2] - bounds_array[0],
            bounds_array[3] - bounds_array[0]
        ]

        points += [bounds_array[0] + vectors[0] + vectors[1]]
        points += [bounds_array[0] + vectors[0] + vectors[2]]
        points += [bounds_array[0] + vectors[1] + vectors[2]]
        points += [bounds_array[0] + vectors[0] + vectors[1] + vectors[2]]

        points = np.array(points)

        edges = [[points[0], points[3], points[5], points[1]],
                 [points[1], points[5], points[7], points[4]],
                 [points[4], points[2], points[6], points[7]],
                 [points[2], points[6], points[3], points[0]],
                 [points[0], points[2], points[4], points[1]],
                 [points[3], points[6], points[7], points[5]]]

        faces = Poly3DCollection(edges, linewidths=1, edgecolors='k')
        faces.set_facecolor((0, 0, 0, 0))

        ax.add_collection3d(faces)

    def draw_sphere(self, ax, xs, ys, zs, c, r) -> None:
        """Draw spheres on atomic sites."""

        phi = np.linspace(0, 2 * np.pi, n_phi)
        theta = np.linspace(0, np.pi, n_theta)

        x = r * np.outer(cos(phi), sin(theta)) + xs
        y = r * np.outer(sin(phi), sin(theta)) + ys
        z = r * np.outer(np.ones(np.size(phi)), cos(theta)) + zs

        # ax.plot_surface(x, y, z, rstride=1, cstride=1, color=c, antialiased=False)
        ax.plot_wireframe(x, y, z, rstride=1, cstride=1, color=c)

    def plot_stats(self, temperatures: List[int], figsize=(8, 6)) -> None:
        """Plot statistical results."""
        P = [norm(i) * 1e6 for i in self.system.avg_p]
        plt.figure(figsize=figsize)
        plt.title("Ferro-Electric Ising Model", fontsize=30, pad=20)
        plt.scatter(temperatures, P, color='IndianRed', marker='.')
        plt.xlabel("Temperature (K)", fontsize=20, labelpad=20)
        plt.ylabel(r"Polarization ($\mathrm{\mu D\cdot\AA^{-3}}$)",
                   fontsize=20,
                   labelpad=20)
        plt.xlim(0, max(temperatures))
        plt.ylim(min(P), max(P))
        plt.tick_params(axis='both', labelsize=10)
        plt.tight_layout()


# =============================================================================
#                          General Plotting Options
# =============================================================================

font = {'family': 'arial', 'weight': 'normal', 'size': 22}
rc('font', **font)

# =============================================================================
#                                Main Program
# =============================================================================

n_phi = 15
n_theta = 15

loc = "tests/SbTaO4"

system = System(loc, supercell_dim=[2, 2, 1])
plotter = Plotter(system)
plotter.draw_lattice(title='Initial Lattice')

# define an Ising Simulator
ising = ElectricIsing(system)

# define temperature (broken into several ranges for resolution)
TR1 = np.linspace(100, 1000, 10)
TR2 = np.linspace(1025, 2000, 40)
TR3 = np.linspace(2200, 4000, 10)
temperatures = np.concatenate((TR1, TR2, TR3))
# T = TR1 # tester

# run simulation
delta_angle = np.radians(0.5)  # perturbation cutoff angle
ising.run_simulation(temperatures, 1, 1, delta_angle)
plotter.draw_lattice(title='Final Lattice')

# plot statistical results
plotter.plot_stats(temperatures)

plt.show()
