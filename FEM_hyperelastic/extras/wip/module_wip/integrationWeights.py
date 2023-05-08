def integrationWeights(ncoord, nelnodes, npoints, elident):

    w = [0 for i in range(npoints)]

    # 1D elements
    if ncoord == 1:
        if npoints == 1:
            w[0] = 2
        elif npoints == 2:
            w = [1, 1]
        elif npoints == 3:
            w = [0.555555555, 0.888888888, 0.555555555]

    # 2D elements
    elif ncoord == 2:

        # Triangular element
        if (nelnodes == 3 or nelnodes == 6):
            if npoints == 1:
                w[0] = 0.5
            elif npoints == 3:
                w[0] = 1/6
                w[1] = 1/6
                w[2] = 1/6
            elif npoints == 4:
                w = [-27/96, 25/96, 25/96, 25/96]

        # Rectangular element
        elif (nelnodes == 4 or nelnodes == 8):
            if npoints == 1:
                w[0] = 4
            elif npoints == 4:
                w = [1, 1, 1, 1]
            elif npoints == 9:
                w1D = [0.555555555, 0.888888888, 0.55555555555]
                for j in range(3):
                    for i in range(3):
                        n = 3*(j-1) + i
                        w[n] = w1D[i] * w1D[j]
    
    # 3D elements
    elif ncoord == 3:
        if (nelnodes == 4 or nelnodes == 10):
            if npoints == 1:
                w[0] = 1/6
            elif npoints == 4:
                w = [1/24, 1/24, 1/24, 1/24]
        elif (nelnodes == 8 or nelnodes == 20):
            if npoints == 1:
                w[0] = 8
            elif npoints == 8:
                w = [1, 1, 1, 1, 1, 1, 1, 1]
            elif npoints == 27:
                w1D = [0.555555555, 0.888888888, 0.55555555555]
                for k in range(3):
                    for j in range(3):
                        for i in range(3):
                            n = 9*(k-1) + 3*(j-1) + i
                            w[n] = w1D[i] * w1D[j] * w1D[k]

    return w