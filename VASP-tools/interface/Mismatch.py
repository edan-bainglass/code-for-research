import numpy as np


def calcErr(params):

    args = list(zip(params[0], params[1]))

    # calculate interface lattice parameter averages
    avg = [np.average(pair) for pair in args]

    # calculate percent errors
    err = [np.abs(l - avg[i]) / l * 100 for i, l in enumerate(args)]

    return [avg, err]


with open(r'Interface\unitCells.txt', 'r') as f:
    latParams = [[float(i) for i in line.split()] for line in f.readlines()]

scale = [[2, 2, 1], [1, 2, 3]]  # CAREFUL!!! Not the same for all orientations

for i, c1 in enumerate(latParams[0]):  # material 1 lattice parameters

    # set interface parameters for material 1
    a1 = latParams[0][(i + 1) % 3] * scale[0][(i + 1) % 3]
    b1 = latParams[0][(i + 2) % 3] * scale[0][(i + 2) % 3]

    for j, c2 in enumerate(latParams[1]):  # material 1 lattice parameters

        # set interface parameters for material 2
        a2 = latParams[1][(j + 1) % 3] * scale[1][(j + 1) % 3]
        b2 = latParams[1][(j + 2) % 3] * scale[1][(j + 2) % 3]

        # create strings for interface planes
        plane1 = ['1' if z == i else '0' for z in range(3)]
        plane2 = ['1' if z == j else '0' for z in range(3)]

        for k in range(2):  # two orientations

            # alternate between and mark plane orientation for material 1
            p1 = [a1, b1] if k == 0 else [b1, a1]
            s = '' if k == 0 else '*'

            # calculate percent errors
            [avg, err] = calcErr([p1, [a2, b2]])
            r = 'GOOD!!!' if all(e < 5 for l in err for e in l) else ''

            if r == 'GOOD!!!':

                print('\n(%s)%s - (%s): %4s\n' %
                      (' '.join(plane1), s, ' '.join(plane2), r))

                print('a -> {0:7.3f} {1:7.3f} {2:7.3f} {3:7.2f}% {4:6.2f}%'.
                      format(p1[0], a2, avg[0], *err[0]) + '\n'
                      'b -> {0:7.3f} {1:7.3f} {2:7.3f} {3:7.2f}% {4:6.2f}%'.
                      format(p1[1], b2, avg[1], *err[1]) + '\n' +
                      'c -> {0:7.3f} {1:7.3f}'.format(c1, c2))
