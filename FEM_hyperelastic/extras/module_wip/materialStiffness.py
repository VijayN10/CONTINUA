import numpy as np

def materialStiffness(ndof, ncoord, B, J, materialprops):
    mu1 = materialprops[0]
    K1 = materialprops[1]

    dl = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    if ncoord == 2:
        Bqq = B[0, 0] + B[1, 1] + 1
        C = np.zeros((ndof,ncoord,ndof,ncoord))
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        C[i, j, k, l] = mu1 * (dl[i, k] * B[j, l] + B[i, l] * dl[j, k] - (2/3) * (B[i, j] * dl[k, l] + dl[i, j] * B[k, l]) + (2/3) * Bqq * dl[i, j] * dl[k, l] / 3) / J**(2/3) + K1 * (2 * J - 1) * J * dl[i, j] * dl[k, l]


    return C