#====================== No. nodes on element faces ================
#
#   This procedure returns the number of nodes on each element face
#   for various element types.  This info is needed for computing
#   the surface integrals associated with the element traction vector
#@jit(nopython=False, parallel=True)
def nfacenodes(ncoord,nelnodes):
    if (ncoord == 2): 
        if (nelnodes == 3 or nelnodes == 4):
            n = 2
        elif (nelnodes == 6 or nelnodes == 8):
            n = 3
    elif (ncoord == 3): 
        if (nelnodes == 4):
            n = 3
        elif (nelnodes == 10):
            n = 6
        elif (nelnodes == 8):
            n = 4
        elif (nelnodes == 20):
            n = 8
    return n