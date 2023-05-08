#======================= Lists of nodes on element faces =============
#
#    This procedure returns the list of nodes on an element face
#    The nodes are ordered so that the element face forms either
#    a 1D line element or a 2D surface element for 2D or 3D problems
#@jit(nopython=False, parallel=True)
def facenodes(ncoord,nelnodes,elident,face):
    i3 = [2,3,1]
    i4 = [2,3,4,1]
    fnlist = np.zeros((nfacenodes(ncoord,nelnodes),1))

    if (ncoord == 2): 
        if (nelnodes == 3): 
            fnlist[0] = face
            fnlist[1] = i3[face-1]
        elif (nelnodes == 6): 
            fnlist[0] = face
            fnlist[1] = i3[face-1]
            fnlist[2] = face+3
        elif (nelnodes==4): 
            fnlist[0] = face
            fnlist[1] = i4[face-1]
        elif (nelnodes==8): 
            fnlist[0] = face
            fnlist[1] = i4[face-1]
            fnlist[2] = face+4
    return fnlist