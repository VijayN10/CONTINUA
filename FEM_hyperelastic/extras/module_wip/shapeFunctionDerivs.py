import numpy as np

def shapeFunctionDerivs(nelnodes,ncoord,elident,xi):

    dNdxi = np.zeros((nelnodes,ncoord))

    #
    # 1D elements
    #
    if (ncoord == 1):
        if (nelnodes==2):
            dNdxi[0,0] = 0.5
            dNdxi[1,0] = -0.5
        elif (nelnodes == 3):
            dNdxi[0,0] = -0.5+xi[0][0]
            dNdxi[1,0] = 0.5+xi[0][0]
            dNdxi[2,0] = -2.*xi[0][0]

    #
    # 2D elements
    #
    elif (ncoord == 2):

        # Triangular element
        if ( nelnodes == 3 ):
            dNdxi[0,0] = 1.
            dNdxi[1,1] = 1.
            dNdxi[2,0] = -1.
            dNdxi[2,1] = -1.
        elif ( nelnodes == 6 ):
            xi3 = 1.-xi[0][0]-xi[1][0]
            dNdxi[0,0] = 4.*xi[0][0]-1.
            dNdxi[1,1] = 4.*xi[1][0]-1.
            dNdxi[2,0] = -(4.*xi3-1.)
            dNdxi[2,1] = -(4.*xi3-1.)
            dNdxi[3,0] = 4.*xi[1][0]
            dNdxi[3,1] = 4.*xi[0][0]
            dNdxi[4,0] = -4.*xi[1][0]
            dNdxi[4,1] = -4.*xi[0][0]
            dNdxi[5,0] = 4.*xi3 - 4.*xi[0][0]
            dNdxi[5,1] = 4.*xi3 - 4.*xi[1][0]

        # Rectangular element
        elif ( nelnodes == 4 ):
            dNdxi[0,0] = -0.25*(1.-xi[1][0])
            dNdxi[0,1] = -0.25*(1.-xi[0][0])
            dNdxi[1,0] = 0.25*(1.-xi[1][0])
            dNdxi[1,1] = -0.25*(1.+xi[0][0])
            dNdxi[2,0] = 0.25*(1.+xi[1][0])
            dNdxi[2,1] = 0.25*(1.+xi[0][0])
            dNdxi[3,0] = -0.25*(1.+xi[1][0])
            dNdxi[3,1] = 0.25*(1.-xi[0][0])

        # nelnodes == 8 CONDITION NEED TO ADD
    return dNdxi