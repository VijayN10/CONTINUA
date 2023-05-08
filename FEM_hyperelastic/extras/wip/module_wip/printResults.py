def printResults(outfile, nprops, materialprops, ncoord, ndof, nnode, coords, nelem, maxnodes, connect, nelnodes, elident, nfix, fixnodes, ndload, dloads, dofs):

    outfile.write("Nodal Displacements: \n")
    if ndof == 2:
        outfile.write(" Node    Coords          u1          u2 \n")
        for i in range(nnode):
            outfile.write("{:3d} {:8.4f} {:8.4f} {:8.4f} {:8.4f}\n".format(
                i, coords[0, i], coords[1, i], dofs[2*i][0], dofs[2*i+1][0]
            ))

    outfile.write('\n\n Strains and Stresses and Deformation Gradient \n')

    lmncoord = np.zeros((ncoord,maxnodes))
    displacements = np.zeros((ndof,maxnodes))

    for lmn in range(nelem):

        outfile.write(' \n Element; {} '.format(lmn))
        if ncoord == 2:
            outfile.write('  \n int pt    Coords          B_11      B_22     B_12      s_11       s_22      s_12      F_11      F_12      F_21      F_22 \n')
        
    
        for a in range(1, nelnodes + 1):
            for i in range(1, ncoord + 1):
                lmncoord[i - 1][a - 1] = coords[i - 1][connect[a - 1][lmn - 1] - 1]
            for i in range(1, ndof + 1):
                displacements[i - 1][a - 1] = dofs[ndof * (connect[a - 1][lmn - 1] - 1) + i - 1]

        n = nelnodes
        ident = elident

        npoints = numberOfIntegrationPoints(ncoord,n, elident)
        dNdx = np.zeros((n,ncoord))
        dxdxi = np.zeros((ncoord,ncoord))
        xi = np.zeros((ncoord,1))
        x = np.zeros((ncoord,1))
        F = np.zeros((ncoord,ncoord))

        # Set up integration points
        xilist = integrationPoints(ncoord,n,npoints, elident)

        for intpt in range(npoints):

            for i in range(ncoord):
                xi[i] = xilist[i][intpt]

            N = shapeFunctions(nelnodes, ncoord, elident, xi)
            dNdxi = shapeFunctionDerivs(nelnodes, ncoord, elident, xi)

            for i in range(ncoord):
                x[i][0] = 0
                for a in range(n):
                    x[i][0] += lmncoord[i][a]*N[a][0]

            for i in range(ncoord):
                for j in range(ncoord):
                    dxdxi[i][j] = 0
                    for a in range(nelnodes):
                        dxdxi[i][j] += lmncoord[i][a] * dNdxi[a][j]
            
            dxidx = np.linalg.inv(dxdxi)
            dt = np.linalg.det(dxdxi)

  
            for a in range(nelnodes):
                for i in range(ncoord):
                    dNdx[a][i] = 0
                    for j in range(ncoord):
                        dNdx[a][i] += dNdxi[a][j] * dxidx[j][i]
            

            for i in range(ncoord):
                for j in range(ncoord):
                    F[i][j] = 0
                    if i == j:
                        F[i][i] = 1
                    for a in range(nelnodes):
                        F[i][j] += displacements[i][a] * dNdx[a][j]

            J = np.linalg.det(F)
            B = np.matmul(F, np.transpose(F))


            Finv = np.linalg.inv(F)
            dNdxs = np.zeros((nelnodes, ncoord))
            for a in range(nelnodes):
                for i in range(ncoord):
                    for j in range(ncoord):
                        dNdxs[a, i] += dNdx[a, j] * Finv[j, i]

            stress = kirchhoffStress(ndof,ncoord,B,J,materialprops)
            stress = stress/J

            if ncoord == 2:
                # print("{:5d} {:7.4f} {:7.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f}".format(intpt, x[0][0], x[1][0], B[0][0], B[1][1], B[0][1], stress[0][0], stress[1][1], stress[0][1]), file=outfile)
                print("{:5d} {:7.4f} {:7.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f}".format(intpt, x[0][0], x[1][0], B[0][0], B[1][1], B[0][1], stress[0][0], stress[1][1], stress[0][1], F[0][0], F[1][0], F[1][0], F[1][1] ), file=outfile)
