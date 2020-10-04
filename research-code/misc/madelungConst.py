import numpy as np

"Calculate Madelung constant by cubic summation"


def BCC(L):
    M = 0
    step = 0.5
    for n in np.arange(step, L + step, step):
        for i in np.arange(-n, n + step, step):
            for j in np.arange(-n, n + step, step):
                for k in np.arange(-n, n + step, step):
                    # print(n,i,j,k)
                    if abs(i) != n and abs(j) != n and abs(k) != n: continue
                    if any([int(l) != l for l in [i, j, k]
                            ]) and not all([int(l) != l for l in [i, j, k]]):
                        continue
                    print(n, i, j, k, 2 * (i + j + k),
                          (i**2 + j**2 + k**2)**-0.5)
                    print()
                    M += (-1)**abs(2 * (i + j + k)) / (i**2 + j**2 + k**2)**0.5
    print(M)


def FCC(L):
    M = 0
    step = 0.5
    for n in np.arange(step, L + step, step):
        for i in np.arange(-n, n + step, step):
            for j in np.arange(-n, n + step, step):
                for k in np.arange(-n, n + step, step):
                    # print(n,i,j,k)
                    if abs(i) != n and abs(j) != n and abs(k) != n: continue
                    if any([int(l) != l for l in [i, j, k]
                            ]) and not all([int(l) != l for l in [i, j, k]]):
                        continue
                    # print(n,i,j,k,2*(i+j+k),(i**2+j**2+k**2)**0.5)
                    # print()
                    M += (-1)**abs(2 * (i + j + k)) / (i**2 + j**2 + k**2)**0.5
    print(M)


def SC(L):
    M = 0
    step = 1
    for n in np.arange(step, L + step, step):
        for i in np.arange(-n, n + step, step):
            for j in np.arange(-n, n + step, step):
                for k in np.arange(-n, n + step, step):
                    if abs(i) != n and abs(j) != n and abs(k) != n: continue
                    M += (-1)**abs(i + j + k) / (i**2 + j**2 + k**2)**0.5
    print(M)
