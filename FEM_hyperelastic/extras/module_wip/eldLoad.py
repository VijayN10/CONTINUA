


def eldload(ncoord, ndof, nfacenodes, elident, coords, traction):
    npoints = numberOfIntegrationPoints(ncoord-1, nfacenodes, elident)
    xi =  np.zeros((ncoord-1,1))
    dxdxi = np.zeros((ncoord,ncoord-1))
    r = np.zeros((ndof*nfacenodes,1))
   
    xilist = integrationPoints(ncoord-1, nfacenodes, npoints, elident)
    w = integrationWeights(ncoord-1, nfacenodes, npoints, elident)

    for intpt in range(npoints):
        for i in range(ncoord-1):
            xi[i] = xilist[i][intpt]
        N = shapeFunctions(nfacenodes, ncoord-1, elident, xi)
        dNdxi = shapeFunctionDerivs(nfacenodes, ncoord-1, elident, xi)

        # Compute the jacobian matrix && its determinant

        for i in range(ncoord):
            for j in range(ncoord-1):
                dxdxi[i][j] = 0
                for a in range(nfacenodes):
                    dxdxi[i][j] += coords[i][a] * dNdxi[a][j]
   
        if ncoord == 2:
            dt = math.sqrt(dxdxi[0][0]**2 + dxdxi[1][0]**2)
        elif ncoord == 3:
            dt = math.sqrt(((dxdxi[1][0]*dxdxi[2][1]) - (dxdxi[1][1]*dxdxi[2][0]))**2 + ((dxdxi[0][0]*dxdxi[2][1]) - (dxdxi[0][1]*dxdxi[2][0]))**2 + ((dxdxi[0][0]*dxdxi[1][1]) - (dxdxi[0][1]*dxdxi[1][0]))**2)
   
        for a in range(nfacenodes):
            for i in range(ndof):
                row = ndof * (a-1) + i
                r[row] += N[a] * traction[i] * w[intpt] * dt
    return r