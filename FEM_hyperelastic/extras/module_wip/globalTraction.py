def globalTraction(ncoord, ndof, nnodes, ndload, coords, nelnodes, elident, connect, dloads, dofs):
    r = np.zeros((ndof*nnodes, 1))
    traction = np.zeros((ndof, 1))

    for load in range(ndload):

        # Extract the coords of the nodes on the appropriate element face
        
        lmn = dloads[0, load]
        face = dloads[1, load]
        n = nelnodes
        ident = elident
        nfnodes = nfacenodes(ncoord, n)
        nodelist = facenodes(ncoord, n, ident, face)
        lmncoord = np.zeros((ncoord, nfnodes))
        for a in range(nfnodes):
            for i in range(ncoord):
                lmncoord[i, a] = coords[i, connect[int(nodelist[a]-1), int(dloads[0, load]-1)]-1]
            # for i in range(ndof):
            #     lmndof[i, a] = dofs[ndof*(connect[int(nodelist[a]), int(dloads[0, load])]-1)+i]
        
        # Compute the element load vector

        for i in range(ndof):
            traction[i] = dloads[i+2, load]
        rel = eldload(ncoord, ndof, nfnodes, ident, lmncoord, traction)

        # Assemble the element load vector into global vector
        
        for a in range(nfnodes):
            for i in range(ndof):
                rw = (connect[int(nodelist[a])-2, int(dloads[0, load])-1])*ndof+i
                r[rw] = r[rw] + rel[(a-1)*ndof+i]
    return r