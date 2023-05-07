def globalResidual(ncoord, ndof, nnode, coords, nelem, maxnodes, elident, nelnodes, connect, materialprops, dofs):
    
    # Assemble the global stiffness matrix

    resid = np.zeros((ndof*nnode,1))
    lmncoord = np.zeros((ncoord,maxnodes))
    lmndof = np.zeros((ndof,maxnodes))
    rel = np.zeros((ndof*maxnodes,ndof*maxnodes))

    # Loop over all the elements

    for lmn in range(1, nelem + 1):

        # Extract coords of nodes, DOF for the current element

        for a in range(1, nelnodes + 1):
            for i in range(1, ncoord + 1):
                lmncoord[i - 1][a - 1] = coords[i - 1][connect[a - 1][lmn - 1] - 1]
            for i in range(1, ndof + 1):
                lmndof[i - 1][a - 1] = dofs[ndof * (connect[a - 1][lmn - 1] - 1) + i - 1]

        n = nelnodes
        ident = elident[lmn - 1]
        rel = elresid(ncoord, ndof, n, ident, lmncoord, materialprops, lmndof)

        # Add the current element residual to the global residual

        for a in range(1, nelnodes + 1):
            for i in range(1, ndof + 1):
                rw = ndof * (connect[a - 1][lmn - 1] - 1) + i - 1
                resid[rw] += rel[ndof * (a - 1) + i - 1]
    return resid