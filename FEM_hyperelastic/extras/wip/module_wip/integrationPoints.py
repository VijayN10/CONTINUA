import numpy as np

def integrationPoints(ncoord, nelnodes, npoints, elident):
    # xi = np.zeros((npoints,ncoord))
    xi=np.zeros((ncoord,npoints))

#   1D elements
    if ncoord == 1:
        if npoints == 1:
            xi[0][0] = 0.
        elif npoints == 2:
            xi[0][0] = -0.5773502692
            xi[0][1] = -xi[0][0]
        elif npoints == 3:
            xi[0][0] = -0.7745966692
            xi[0][1] = 0.0
            xi[0][2] = -xi[0][0]

#   2D elements
    elif ncoord == 2:

    #   Triangular element
        if nelnodes == 3 or nelnodes == 6:
            if npoints == 1:
                xi[0][0] = 1. / 3.
                xi[1][0] = 1. / 3.
            elif npoints == 3:
                xi[0][0] = 0.6
                xi[1][0] = 0.2
                xi[0][1] = 0.2
                xi[1][1] = 0.6
                xi[0][2] = 0.2
                xi[1][2] = 0.2
            elif npoints == 4:
                xi[0][0] = 1. / 3.
                xi[1][0] = 1. / 3.
                xi[0][1] = 0.6
                xi[1][1] = 0.2
                xi[0][2] = 0.2
                xi[1][2] = 0.6
                xi[0][3] = 0.2
                xi[1][3] = 0.2

    #   Rectangular element
        elif nelnodes == 4 or nelnodes == 8:
            if npoints == 1:
                xi[0][0] = 0.
                xi[1][0] = 0.
            elif npoints == 4:
                xi[0][0] = -0.5773502692
                xi[1][0] = xi[0][0]
                xi[0][1] = -xi[0][0]
                xi[1][1] = xi[0][0]
                xi[0][2] = xi[0][0]
                xi[1][2] = -xi[0][0]
                xi[0][3] = -xi[0][0]
                xi[1][3] = -xi[0][0]
            elif npoints == 9:
                xi[0][0] = -0.7745966692
                xi[1][0] = xi[0][0]
                xi[0][1] = 0.0
                xi[1][1] = xi[0][0]
                xi[0][2] = -xi[0][0]
                xi[1][2] = xi[0][0]
                xi[0][3] = xi[0][0]
                xi[0][4] = 0.0
                xi[1][4] = 0.0
                xi[0][5] = -xi[0][0]
                xi[1][5] = 0.0
                xi[0][6] = xi[0][0]
                xi[1][6] = -xi[0][0]
                xi[0][7] = 0
                xi[1][7] = -xi[0][0]
                xi[0][8] = -xi[0][0]
                xi[1][8] = -xi[0][0]

    #   3D elements not added

    return xi