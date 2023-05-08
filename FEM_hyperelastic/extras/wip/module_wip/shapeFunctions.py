import numpy as np

def shapeFunctions(nelnodes, ncoord, elident, xi):
    
    N = np.zeros((nelnodes,1))
    # 
    #   1D elements
    # 
    if (ncoord == 1):
        if (nelnodes==2):
            N[0] = 0.5*(1+xi[0][0])
            N[1] = 0.5*(1-xi[0][0])
        elif (nelnodes == 3):
            N[0] = -0.5*xi[0]*(1-xi[0][0])
            N[1] =  0.5*xi[0]*(1+xi[0][0])
            N[2] = (1-xi[0])*(1+xi[0][0])
#
#  2D elements
#
    elif (ncoord == 2):
  # 
  #     Triangular element
  # 
        if ( nelnodes == 3 ):
            N[0] = xi[0][0]
            N[1] = xi[0][0]
            N[2] = 1-xi[0][0]-xi[1][0]              
        elif ( nelnodes == 6 ):
            xi3 = 1-xi[0][0]-xi[1][0]
            N[0] = (2*xi[0][0]-1)*xi[0][0]
            N[1] = (2*xi[1][0]-1)*xi[1][0]
            N[2] = (2*xi3-1)*xi3
            N[3] = 4*xi[0][0]*xi[1][0]
            N[4] = 4*xi[1][0]*xi3
            N[5] = 4*xi3*xi[0][0]
  # 
  #     Rectangular element
  #                   
        elif ( nelnodes == 4 ):
            N[0] = 0.25*(1-xi[0][0])*(1-xi[1][0])
            N[1] = 0.25*(1+xi[0][0])*(1-xi[1][0])
            N[2] = 0.25*(1+xi[0][0])*(1+xi[1][0])
            N[3] = 0.25*(1-xi[0][0])*(1+xi[1][0])
        elif (nelnodes == 8):
            N[0] = -0.25*(1-xi[0][0])*(1-xi[1][0])*(1+xi[0][0]+xi[1][0])
            N[1] = 0.25*(1+xi[0][0])*(1-xi[1][0])*(xi[0][0]-xi[1][0]-1)
            N[2] = 0.25*(1+xi[0][0])*(1+xi[1][0])*(xi[0][0]+xi[1][0]-1)
            N[3] = 0.25*(1-xi[0][0])*(1+xi[1][0])*(xi[1][0]-xi[0][0]-1)
            N[4] = 0.5*(1-xi[0][0]*xi[0][0])*(1-xi[1][0])
            N[5] = 0.5*(1+xi[0][0])*(1-xi[1][0]*xi[1][0])
            N[6] = 0.5*(1-xi[0][0]*xi[0][0])*(1+xi[1][0])
            N[7] = 0.5*(1-xi[0][0])*(1-xi[1][0]*xi[1][0])

    return N