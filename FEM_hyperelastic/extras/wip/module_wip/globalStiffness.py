import numpy as np
from elStif import elStif


def globalStiffness(ncoord, ndof, nnode, coords, nelem, maxnodes, elident, nelnodes, connect, materialprops, dofs):
    
    # Assemble the global stiffness matrix

    stif = np.zeros((ndof*nnode,ndof*nnode))
    lmncoord = np.zeros((ncoord,maxnodes))
    lmndof = np.zeros((ndof,maxnodes))
    
    # Loop over all the elements

    for lmn in range(nelem):
        
        # Extract coords of nodes, DOF for the current element

        for a in range(nelnodes):
            for i in range(ncoord):
                lmncoord[i][a] = coords[i][connect[a][lmn]-1]
            for i in range(ndof):
                lmndof[i][a] = dofs[ndof*(connect[a][lmn]-1)+i]
        n = nelnodes
        ident = elident[lmn]
        kel = elStif(ncoord, ndof, n, ident, lmncoord, materialprops, lmndof)
        
        # Add the current element stiffness:the global stiffness

        for a in range(nelnodes):
            for i in range(ndof):
                for b in range(nelnodes):
                    for k in range(ndof):
                        rw = ndof*(connect[a][lmn]-1)+i
                        cl = ndof*(connect[b][lmn]-1)+k
                        stif[rw][cl] += kel[ndof*(a-1)+i][ndof*(b-1)+k]
    return stif