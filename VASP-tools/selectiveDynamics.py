# -*- coding: utf-8 -*-

loc = 'C:/Users/edanb/Desktop'
filename = 'tempPOS.vasp'

axis = 1  # 0, 1, 2 = x, y, z
l = 43.2661589986
coordType = 'Cartesian'

with open(loc + '/' + filename, 'r') as f, open(loc + '/SD.vasp', 'w') as g:

    for line in range(5):
        g.write(f.readline())

    atoms = f.readline()
    g.write(atoms)
    atoms = atoms.split()
    nums = f.readline()
    g.write(nums)
    nums = nums.split()
    g.write('SD\n')
    g.write(f.readline())

    for s in range(len(nums)):

        lines = [[float(value) for value in f.readline().split()]
                 for line in range(int(nums[s]))]

        for line in sorted(lines, key=lambda line: line[axis]):
            SD = ' F F F' if .17 * l < line[axis] < .36 * l else ' T T T'
            g.write('%18.9f %18.9f %18.9f' % tuple(line) + SD + ' \n')
