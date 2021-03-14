import re
import inspect

import numpy as np
from numpy import exp, dot, sin, cos, arccos, arctan2
from numpy.random import random
from numpy.linalg import norm

from pymatgen.io.vasp import Potcar
from pymatgen.core import structure
import pymatgen.analysis.ferroelectricity.polarization as ferro

import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# =============================================================================
#                                   Functions
# =============================================================================


def draw_sphere(ax, xs, ys, zs, c, r):
    ''' Draw spheres on atomic sites '''

    phi = np.linspace(0, 2 * np.pi, n_phi)
    theta = np.linspace(0, np.pi, n_theta)

    x = r * np.outer(cos(phi), sin(theta)) + xs
    y = r * np.outer(sin(phi), sin(theta)) + ys
    z = r * np.outer(np.ones(np.size(phi)), cos(theta)) + zs

    # ax.plot_surface(x, y, z, rstride=1, cstride=1, color=c, antialiased=False)
    ax.plot_wireframe(x, y, z, rstride=1, cstride=1, color=c)


def plot_cell(ax, bounds):
    ''' Draw cell boundaries '''

    bounds_array = [np.array(list(i)) for i in bounds]

    points = []
    points += bounds_array
    vectors = [
        bounds_array[1] - bounds_array[0], bounds_array[2] - bounds_array[0],
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


def draw_lattice(lat,
                 P=np.array([]),
                 ax=None,
                 tilt=0,
                 turn=0,
                 aspect=1.85,
                 oxy=False,
                 title=None,
                 labels=False):
    ''' Draws lattice with superimposed dipole moments '''

    A, B, C = struct.lattice.a, struct.lattice.b, struct.lattice.c

    if inspect.stack()[1][3] == '<module>': plt.close()

    if not ax:
        fig = plt.figure(figsize=(5, 10))
        ax = fig.add_subplot(111, projection='3d')

    ax.set_axis_off()
    ax.set_xlim3d(0, A)
    ax.set_ylim3d(0, B)
    ax.set_zlim3d(0, C)

    # plot cell bounderies
    bounds = [(0, 0, 0), (0, B, 0), (A, 0, 0), (0, 0, C)]
    plot_cell(ax, bounds)

    # add oxygens if requested
    if oxy:
        lat = np.array(list(lat) + list(oxygens))

    for i in enumerate(lat):
        if i[0] in struct.indices_from_symbol('Sb'):
            #            continue
            c, r = 'sandybrown', .4
        elif i[0] in struct.indices_from_symbol('Ta'):
            #            continue
            c, r = 'darkkhaki', .6
        else:
            #            continue
            c, r = 'darkred', .15
        x = i[1][0][0]
        y = i[1][0][1]
        z = i[1][0][2]

        draw_sphere(ax, x, y, z, c, r)  # draw atom
        if labels:
            ax.text(x + .25, y + .25, z + .25, i[1][2],
                    fontsize=10)  # add labels

        if i[0] not in struct.indices_from_symbol('O'):
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
                      length=mag * 2.5,
                      color='k',
                      arrow_length_ratio=0.1,
                      pivot='middle')  # draw dipole moment

    if P.any():
        P = P / norm(P)
        ax.quiver(A / 2,
                  B / 2,
                  C / 2,
                  P[0],
                  P[1],
                  P[2],
                  length=4.,
                  color='g',
                  arrow_length_ratio=0.1,
                  pivot='middle')  # draw polarization vector

    ax.set_aspect('auto')  # scale plot
    ax.view_init(tilt, turn)  # set viewing direction
    if title: ax.set_title(title, {'fontsize': 20})


def draw_combo(s=None, tilt=0, turn=0, oxy=False, aspect=1.85):
    ''' Plot initial vs. selected lattice '''

    if s is None:
        t = T[-1]
        lat = cations
        P = avgP[-1]
    else:
        t = T[s]
        lat = avgCell[s]
        P = avgP[s]

    plt.close()
    fig = plt.figure(figsize=(10, 10))
    ax0 = fig.add_subplot(121, projection='3d', title='Initial Lattice')
    draw_lattice(init_cations, init_P, ax0, tilt, turn, aspect, oxy)
    ax1 = fig.add_subplot(122, projection='3d', title=f'T = {t} K')
    draw_lattice(lat, P, ax1, tilt, turn, aspect, oxy)


def getCharges():
    ''' Extract Bader charges from Bader analysis output '''

    POTCAR = Potcar.from_file(loc + '/POTCAR')
    z = ferro.zval_dict_from_potcar(POTCAR)
    bader = open(loc + '/ACF.dat')
    c = []
    s = [line.split() for line in bader if re.match('[0-9]+', line.split()[0])]
    for i in unit.symbol_set:
        for j in unit.indices_from_symbol(i):
            c.append(z[i] - float(s[j][4]))
    bader.close()

    return c


def build_lattice():
    ''' Initialize lattice with calculated dipole moments. '''

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
        c = struct.site_properties['charge']

        for atom in enumerate(struct_no_O):
            # define local molecule
            nn = struct.get_neighbors(atom[1], 2.9, include_index=True)
            sites = [atom[1]] + [nn[i][0] for i in range(len(nn))]
            cm = center_of_mass(sites)
            # calculate dipole
            dip = np.zeros(3)
            for i in enumerate(sites):
                if i[0] == 0:
                    dip += (i[1].coords - cm) * c[atom[0]]  # cation
                else:
                    dip += (i[1].coords - cm) * c[nn[i[0] - 1][2]]  # oxygens
            dips.append(dip)

        for i in range(len(oxy)):
            dips.append([])  # add empty dipoles on oxygen sites

        return dips

    coords = list(struct.cart_coords)
    labels = [str(x[1]) + str(x[0] + 1) for x in enumerate(struct.species)]
    oxy = struct.indices_from_symbol('O')
    dipoles = calc_dipoles()

    # combine coordinates, labels, and dipoles
    lat = np.array([i for p in zip(coords, dipoles, labels)
                    for i in p]).reshape(-1, 3)

    # split cations and oxygens
    cations = lat[:min(oxy)]
    oxygens = lat[min(oxy):]

    return cations, oxygens


def sweep(func, lat):
    ''' Sweep through lattice, perturbing dipole moments and computing the
        energetic cost. If favorable, new dipole moment is retained. Total
        dipole-dipole interaction potential energy and polarizations are
        obtained in the process. [r] = A, [p] = eA, or Debye, [E] = eV '''
    def equilibrate():
        def perturb():
            a = p1.copy()
            # perturb angles
            theta = arccos(p1[2] / norm(p1)) + dAng * (2 * random() - 1)
            phi = arctan2(p1[1], p1[0]) + dAng * (2 * random() - 1)
            # redefine dipole moment
            p1[0] = norm(a) * sin(theta) * cos(phi)
            p1[1] = norm(a) * sin(theta) * sin(phi)
            p1[2] = norm(a) * cos(theta)

        def calcE():
            # calculate dipole-dipole interaction potential
            E = 0.
            for j in range(L):
                if j != i:
                    r2, p2 = lat[j][0], lat[j][1]
                    r = r2 - r1
                    rn = norm(r)
                    rh = r / rn
                    E += 14.4 / rn**3 * (dot(p1, p2) -
                                         3 * dot(p1, rh) * dot(p2, rh))
            return E

        # equilibriate lattice @ current temperature
        r1, p1 = lat[i][0], lat[i][1].copy()
        e1 = calcE()
        perturb()
        e2 = calcE()
        eDif = round(e2 - e1, 10)
        if eDif < 0:
            lat[i][1] = p1
        elif random() <= exp(-eDif / (8.617e-5 * t)):
            lat[i][1] = p1
        else:
            pass

    # perform sweep
    L = len(lat)
    for i in range(L):
        eval(func)


def Ising(t, lat):
    ''' Executes Ising algorithm to obtain average structure polarization '''

    # thermalize lattice
    for n in range(equiSwp):
        sweep('equilibrate()', lat)

    # calculate polarization
    P = np.zeros(3)
    dips = 0.
    for n in range(polSwp):
        sweep('equilibrate()', lat)
        P += sum(lat[:, 1])
        dips += lat[:, 1]

    # normalize
    avgP.append(P / polSwp / struct.volume)
    dips = dips / polSwp

    # save dipole moment configuration at current temperature
    cell = lat.copy()
    for i in range(len(lat)):
        cell[i][1] = dips[i]
    avgCell.append(cell)


# =============================================================================
#                               General Plotting Options
# =============================================================================

font = {'family': 'arial', 'weight': 'normal', 'size': 22}
rc('font', **font)

# =============================================================================
#                                   Main Program
# =============================================================================

n_phi = 20
n_theta = 20

loc = r'tests/SbTaO4'

cell_type = 'super'

# get unit cell
unit = structure.Structure.from_file(loc + '/POSCAR')
unit.add_site_property('charge', getCharges())

if cell_type == 'super':
    # construct super cell
    sup = unit.copy()
    sup.make_supercell([2, 2, 1])
    struct = sup
    aspect = 1
else:
    struct = unit
    aspect = 1.85

# construct cation-only structure
struct_no_O = struct.copy()
struct_no_O.remove_species('O')

# initialize lattice with calculated dipole moments
cations, oxygens = build_lattice()
init_cations = cations.copy()
init_P = sum(init_cations[:, 1]) / struct.volume
draw_lattice(init_cations, title='Initial Lattice', P=init_P, aspect=aspect)

# define equilibriation and polarization sweeps
equiSwp = 1
polSwp = 2

dAng = np.radians(0.5)  # perturbation cutoff angle

# define temperature (broken into several ranges for resolution)
TR1 = np.linspace(100, 1000, 10)
TR2 = np.linspace(1025, 2000, 40)
TR3 = np.linspace(2200, 4000, 10)
T = np.concatenate((TR1, TR2, TR3))
# T = TR1 # tester

# define cell register to store final cell per temperature
avgCell = []

# define polarization register to store average polarizations per temperature
avgP = []

# sweep through the lattice for each temperature step
for t in T:
    Ising(t, cations)

# draw_combo()  # start and end polarization states

# plot statistical results
# plt.close()
# P = [norm(i) * 1e6 for i in avgP]
# fig = plt.figure(figsize=(10, 7))
# plt.title("Ferro-Electric Ising Model", fontsize=30, pad=20)
# plt.scatter(T, P, color='IndianRed', marker='.')
# plt.xlabel("Temperature (K)", fontsize=20, labelpad=20)
# plt.ylabel(u"Polarization ($\mathrm{\mu D\cdot\AA^{-3}}$)",
#            fontsize=20,
#            labelpad=20)
# plt.axis([0, max(T), 0, max(P)], 'tight')
# plt.tick_params(axis='both', labelsize=10)

plt.show()
