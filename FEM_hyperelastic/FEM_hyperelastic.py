import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import math

# INPUT

#          Example 2D and 3D hyperelastic FEM code
#
#        Variables read from input file;
#        nprops              No. material parameters
#        materialprops(i)    List of material parameters
#        ncoord              No. spatial coords (2 for 2D, 3 for 3D)
#        ndof                No. degrees of freedom per node (2 for 2D, 3 for 3D)
#                            (here ndof=ncoord, but the program allows them to be different
#                            to allow extension to plate & beam elements with C^1 continuity)
#        nnode               No. nodes
#        coords(i,j)         ith coord of jth node, for i=1..ncoord; j=1..nnode
#        nelem               No. elements
#        maxnodes            Max no. nodes on any one element (used for array dimensioning)
#        nelnodes(i)         No. nodes on the ith element
#        elident(i)          An integer identifier for the ith element.  Not used
#                            in this code but could be used to switch on reduced integration,
#                            etc.
#        connect(i,j)        List of nodes on the jth element
#        nfix                Total no. prescribed displacements
#        fixnodes(i,j)       List of prescribed displacements at nodes
#                            fixnodes(1,j) Node number
#                            fixnodes(2,j) Displacement component number (1, 2 or 3)
#                            fixnodes(3,j) Value of the displacement
#        ndload              Total no. element faces subjected to tractions
#        dloads(i,j)         List of element tractions
#                            dloads(1,j) Element number
#                            dloads(2,j) face number
#                            dloads(3,j), dloads(4,j), dloads(5,j) Components of traction
#                            (assumed uniform)


# Total no. material parameters and list of parameters

nprops = 2
materialprops = np.array([[1], [10]])     # materialprops[mu1,K1]      # Python first position is 0 --> Need to change all code?
# mu1 = 1
# k1 = 10
mul = materialprops[0][0]
k1 = materialprops[1][0]

# no. coords (1:3), no. DOF, no. nodes and nodal coordinates

ncoord = 2
ndof = 2
nnode = 4

coords = np.array([[0, 1, 1, 0],[0, 0, 1, 1]])

# No. elements and connectivity

nelem = 1
maxnodes = 4
nelnodes = 4
elident = np.array([[1]])                 # elident = np.array([[1]])
connect = np.array([[1],[2],[3],[4]])     # connect = np.array([[1],[2],[3],[4]])

# No. nodes with prescribed displacements, with the prescribed displacements

nfix = 4
fixnodes = np.array([[1, 1, 2, 4],[1, 2, 2, 1],[0, 0, 0, 0]])

# No. loaded element faces, with the loads

ndload = 1
dloads = np.array([[1],[2],[3],[0]])      # dloads = np.array([[1],[2],[3],[0]])
 

# def read_input_file(infile):
#     with open(infile) as f:
#         cellarray = f.read().split()

#     nprops = int(cellarray[2])
#     materialprops = [0]*nprops
#     cellno = 2
#     for i in range(nprops):
#         cellno += 2
#         materialprops[i] = float(cellarray[cellno])

#     cellno += 2
#     ncoord = int(cellarray[cellno])
#     cellno += 2
#     ndof = float(cellarray[cellno])
#     cellno += 2
#     nnode = int(cellarray[cellno])

#     coords = [[0 for i in range(ncoord)] for j in range(nnode)]
#     cellno += 1
#     for i in range(nnode):
#         for j in range(ncoord):
#             cellno += 1
#             coords[i][j] = float(cellarray[cellno])

#     cellno += 2
#     nelem = int(cellarray[cellno])
#     cellno += 2
#     maxnodes = int(cellarray[cellno])
#     connect = [[0 for i in range(nelem)] for j in range(maxnodes)]
#     nelnodes = [0]*nelem
#     elident = [0]*nelem
#     cellno += 3
#     for i in range(nelem):
#         cellno += 1
#         elident[i] = int(cellarray[cellno])
#         cellno += 1
#         nelnodes[i] = int(cellarray[cellno])
#         for j in range(nelnodes[i]):
#             cellno += 1
#             connect[j][i] = int(cellarray[cellno])

#     cellno += 2
#     nfix = int(cellarray[cellno])
#     cellno += 3
#     fixnodes = [[0 for i in range(3)] for j in range(nfix)]
#     for i in range(nfix):
#         cellno += 1
#         fixnodes[i][0] = int(cellarray[cellno])
#         cellno += 1
#         fixnodes[i][1] = int(cellarray[cellno])
#         cellno += 1
#         fixnodes[i][2] = float(cellarray[cellno])

#     cellno += 2
#     ndload = int(cellarray[cellno])
#     cellno += 3
#     dloads = [[0 for i in range(2+ndof)] for j in range(ndload)]
#     for i in range(ndload):
#         cellno += 1
#         dloads[i][0] = int(cellarray[cellno])
#         cellno += 1
#         dloads[i][1] = int(cellarray[cellno])
#         for j in range(ndof):
#             cellno += 1
#             dloads[i][j+2] = float(cellarray[cellno])
    
#     return nprops,materialprops,ncoord,ndof,nnode,coords,nelem,maxnodes,connect,nelnodes,elident,nfix,fixnodes,ndload,dloads



######################################################################################################################################################


#================= Material Stiffness ==================================
#
#    Computes material stiffness tensor C_{ijkl} 
#    Currently coded either for plane strain or general 3D.

def materialstiffness(ndof, ncoord, B, J, materialprops):
    mu1 = materialprops[0]
    K1 = materialprops[1]

    dl = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    if ncoord == 2:
        Bqq = B[0, 0] + B[1, 1] + 1
        C = np.zeros((ndof,ncoord,ndof,ncoord))
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        C[i, j, k, l] = mu1 * (dl[i, k] * B[j, l] + B[i, l] * dl[j, k] - (2/3) * (B[i, j] * dl[k, l] + dl[i, j] * B[k, l]) + (2/3) * Bqq * dl[i, j] * dl[k, l] / 3) / J**(2/3) + K1 * (2 * J - 1) * J * dl[i, j] * dl[k, l]


    return C

######################################################################################################################################################

#================= Stress ==================================
#
#   Computes stress sigma_{ij} given B_{ij}

def Kirchhoffstress(ndof,ncoord,B,J,materialprops):
  stress = np.zeros((ndof,ncoord))
  dl = [[1,0,0],[0,1,0],[0,0,1]]

  mu1 = materialprops[0][0]
  K1 = materialprops[1][0]
  Bkk = sum([B[i][i] for i in range(ndof)])
  if (ndof==2):
      Bkk = Bkk + 1
  for i in range(ndof):
      for j in range(ncoord):
          stress[i][j] = mu1*(B[i][j] - Bkk*dl[i][j]/3.)/J**(2/3) + K1*J*(J-1)*dl[i][j]

  return stress

######################################################################################################################################################

#====================== No. integration points =============================
#
#   Defines the number of integration points:be used for
#   each element type

def numberofintegrationpoints(ncoord, nelnodes, elident):
    if (ncoord == 1):
        n = nelnodes
    elif (ncoord == 2):
        if (nelnodes == 3):
            n = 1
        elif (nelnodes == 6):
            n = 3
        elif (nelnodes == 4):
            n = 4
        elif (nelnodes == 8):
            n = 9
    elif (ncoord == 3):
        if (nelnodes == 4):
            n = 1
        elif (nelnodes == 10):
            n = 4
        elif (nelnodes == 8):
            n = 8
        elif (nelnodes == 20):
            n = 27
    return n

######################################################################################################################################################

#====================== INTEGRATION POINTS ==================================
#
#   Defines positions of integration points

def integrationpoints(ncoord, nelnodes, npoints, elident):
    # xi = np.zeros((npoints,ncoord))
    xi=np.zeros((ncoord,npoints))

#   1D elements
    if ncoord == 1:
        if npoints == 1:
            xi[0][0] = 0.
        elif npoints == 2:
            xi[0][0] = -0.5773502692
            xi[0][1] = -xi[0][0]
        elif npoints == 3:
            xi[0][0] = -0.7745966692
            xi[0][1] = 0.0
            xi[0][2] = -xi[0][0]

#   2D elements
    elif ncoord == 2:

    #   Triangular element
        if nelnodes == 3 or nelnodes == 6:
            if npoints == 1:
                xi[0][0] = 1. / 3.
                xi[1][0] = 1. / 3.
            elif npoints == 3:
                xi[0][0] = 0.6
                xi[1][0] = 0.2
                xi[0][1] = 0.2
                xi[1][1] = 0.6
                xi[0][2] = 0.2
                xi[1][2] = 0.2
            elif npoints == 4:
                xi[0][0] = 1. / 3.
                xi[1][0] = 1. / 3.
                xi[0][1] = 0.6
                xi[1][1] = 0.2
                xi[0][2] = 0.2
                xi[1][2] = 0.6
                xi[0][3] = 0.2
                xi[1][3] = 0.2

    #   Rectangular element
        elif nelnodes == 4 or nelnodes == 8:
            if npoints == 1:
                xi[0][0] = 0.
                xi[1][0] = 0.
            elif npoints == 4:
                xi[0][0] = -0.5773502692
                xi[1][0] = xi[0][0]
                xi[0][1] = -xi[0][0]
                xi[1][1] = xi[0][0]
                xi[0][2] = xi[0][0]
                xi[1][2] = -xi[0][0]
                xi[0][3] = -xi[0][0]
                xi[1][3] = -xi[0][0]
            elif npoints == 9:
                xi[0][0] = -0.7745966692
                xi[1][0] = xi[0][0]
                xi[0][1] = 0.0
                xi[1][1] = xi[0][0]
                xi[0][2] = -xi[0][0]
                xi[1][2] = xi[0][0]
                xi[0][3] = xi[0][0]
                xi[0][4] = 0.0
                xi[1][4] = 0.0
                xi[0][5] = -xi[0][0]
                xi[1][5] = 0.0
                xi[0][6] = xi[0][0]
                xi[1][6] = -xi[0][0]
                xi[0][7] = 0
                xi[1][7] = -xi[0][0]
                xi[0][8] = -xi[0][0]
                xi[1][8] = -xi[0][0]

    #   3D elements not added

    return xi

######################################################################################################################################################

#================= INTEGRATION WEIGHTS ==================================
#
#   Defines integration weights w_i

def integrationweights(ncoord, nelnodes, npoints, elident):

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

######################################################################################################################################################

#================= SHAPE FUNCTIONS ==================================
#
#        Calculates shape functions for various element types

def shapefunctions(nelnodes, ncoord, elident, xi):
    
    N = np.zeros((nelnodes,1))
    # 
    #   1D elements
    # 
    if (ncoord == 1):
        if (nelnodes==2):
            N[0] = 0.5*(1+xi[0])
            N[1] = 0.5*(1-xi[0])
        elif (nelnodes == 3):
            N[0] = -0.5*xi[0]*(1-xi[0])
            N[1] =  0.5*xi[0]*(1+xi[0])
            N[2] = (1-xi[0])*(1+xi[0])
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
            N[2] = 1-xi[0][0]-xi[0][0]              
        elif ( nelnodes == 6 ):
            xi3 = 1-xi[0][0]-xi[0][0]
            N[0] = (2*xi[0][0]-1)*xi[0][0]
            N[1] = (2*xi[0][0]-1)*xi[0][0]
            N[2] = (2*xi3-1)*xi3
            N[3] = 4*xi[0][0]*xi[0][0]
            N[4] = 4*xi[0][0]*xi3
            N[5] = 4*xi3*xi[0][0]
  # 
  #     Rectangular element
  #                   
        elif ( nelnodes == 4 ):
            N[0] = 0.25*(1-xi[0][0])*(1-xi[0][0])
            N[1] = 0.25*(1+xi[0][0])*(1-xi[0][0])
            N[2] = 0.25*(1+xi[0][0])*(1+xi[0][0])
            N[3] = 0.25*(1-xi[0][0])*(1+xi[0][0])
        elif (nelnodes == 8):
            N[0] = -0.25*(1-xi[0][0])*(1-xi[0][0])*(1+xi[0][0]+xi[0][0])
            N[1] = 0.25*(1+xi[0][0])*(1-xi[0][0])*(xi[0][0]-xi[0][0]-1)
            N[2] = 0.25*(1+xi[0][0])*(1+xi[0][0])*(xi[0][0]+xi[0][0]-1)
            N[3] = 0.25*(1-xi[0][0])*(1+xi[0][0])*(xi[0][0]-xi[0][0]-1)
            N[4] = 0.5*(1-xi[0][0]*xi[0][0])*(1-xi[0][0])
            N[5] = 0.5*(1+xi[0][0])*(1-xi[0][0]*xi[0][0])
            N[6] = 0.5*(1-xi[0][0]*xi[0][0])*(1+xi[0][0])
            N[7] = 0.5*(1-xi[0][0])*(1-xi[0][0]*xi[0][0])

    return N

######################################################################################################################################################

#
#================= SHAPE FUNCTION DERIVATIVES ======================
#

def shapefunctionderivs(nelnodes,ncoord,elident,xi):

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

######################################################################################################################################################
#
#================= ELEMENT RESIDUAL VECTOR ================================
#

def elresid(ncoord, ndof, nelnodes, elident, coord, materialprops, displacement):
   
#  Assemble the element residual force
#
#    Arguments:
#
#      ncoord             No. coordinates (2 or 3 for 2D or 3D problem)
#      ndof               No. degrees of freedom per node (often ndof = ncoord)
#      nelnodes           No. nodes on the element
#      elident            Element identifier (not used here - for future enhancements!)
#      coords[i,a]        ith coord of ath node
#      materialprops      Material properties passed on to constitutive procedures
#      displacement[i,a]  ith displacement component at ath node
#
#   Local variables
#      npoints            No. integration points
#      xi[i,inpt]         ith local coord of integration point no. intpt
#      w[intpt]           weight for integration point no. intpt
#      N[a]               Shape function associated with ath node on element
#      dNdxi[a,i]         Derivative of ath shape function wrt ith local coord
#      dNdx[a,i]          Derivative of ath shape function wrt ith global coord
#      dxdxi[i,j]         Derivative of ith global coord wrt jth local coord
#      dxidx[i,j]         Derivative of ith local coord wrt jth global coord
#      det                Determinant of jacobian
#      strain[i,j]        strain_ij components
#      stress[i,j]        stress_ij components
#      r[row]             Residual vector
   
    npoints = numberofintegrationpoints(ncoord, nelnodes, elident)
    dxdxi = np.zeros((ncoord, ncoord))
    dxidx = np.zeros((ncoord, ncoord))
    dNdxs = np.zeros((nelnodes, ncoord))
    rel = np.zeros((ndof * nelnodes, 1))
    xi = np.zeros((ncoord,npoints))
    dNdx = np.zeros((nelnodes, ncoord))
    F = np.zeros((ncoord,ncoord))

# Set up integration points and weights
    
    xilist = integrationpoints(ncoord, nelnodes, npoints, elident)
    w = integrationweights(ncoord, nelnodes, npoints, elident)

    # Loop over the integration points
    for intpt in range(npoints):

        # Compute shape functions and derivatives wrt local coords

        for i in range(ncoord):
            xi[i] = xilist[i][intpt]
        N = shapefunctions(nelnodes, ncoord, elident, xi)
        dNdxi = shapefunctionderivs(nelnodes, ncoord, elident, xi)

        # Compute the jacobian matrix and its determinant

        for i in range(ncoord):
            for j in range(ncoord):
                dxdxi[i][j] = 0
                for a in range(nelnodes):
                    dxdxi[i][j] += coord[i][a] * dNdxi[a][j]
        dxidx = np.linalg.inv(dxdxi)
        dt = np.linalg.det(dxdxi)

        # Convert shape function derivatives to derivatives wrt global coords

        for a in range(nelnodes):
            for i in range(ncoord):
                dNdx[a][i] = 0
                for j in range(ncoord):
                    dNdx[a][i] += dNdxi[a][j] * dxidx[j][i]

        # Compute the deformation gradients by differentiating displacements

        for i in range(ncoord):
            for j in range(ncoord):
                F[i][j] = 0
                if i == j:
                    F[i][i] = 1
                for a in range(nelnodes):
                    F[i][j] += displacement[i][a] * dNdx[a][j]

        # Compute Bbar and J

        J = np.linalg.det(F)
        B = np.dot(F, np.transpose(F))

        # Convert shape function derivatives to derivatives wrt spatial coords

        Finv = np.linalg.inv(F)
        for a in range(nelnodes):
            for i in range(ncoord):
                dNdxs[a][i] = 0
                for j in range(ncoord):
                    dNdxs[a][i] += dNdx[a][j] * Finv[j][i]
        
        # Compute the stress

        stress = Kirchhoffstress(ndof,ncoord,B,J,materialprops)

        # Compute the element residual

        for a in range(nelnodes):
            for i in range(ncoord):
                for row in range(ndof):
                    rel[row + a * ndof] += stress[row][i] * dNdxs[a][i] * w[intpt] * dt
    return rel
       
######################################################################################################################################################
#
#================= ELEMENT STIFFNESS MATRIX ================================
#

def elstif(ncoord, ndof, nelnodes, elident, coord, materialprops, displacement):
    
#  Assemble the element stiffness
#
#    Arguments;
#
#      ncoord             No. coordinates (2 or 3 for 2D or 3D problem)
#      ndof               No. degrees of freedom per node (often ndof = ncoord)
#      nelnodes           No. nodes on the element
#      elident            Element identifier (not used here - for future enhancements!)
#      coords(i,a)        ith coord of ath node
#      materialprops      Material properties passed on:constitutive procedures
#      displacement(i,a)  ith displacement component at ath node
#
#   Local variables
#      npoints            No. integration points
#      xi(i,inpt)         ith local coord of integration point no. intpt
#      w(intpt)           weight for integration point no. intpt
#      N(a)               Shape function associated with ath node on element
#      dNdxi(a,i)         Derivative of ath shape function wrt ith local coord
#      dNdx(a,i)          Derivative of ath shape function wrt ith global coord
#      dxdxi(i,j)         Derivative of ith global coord wrt jth local coord
#      dxidx(i,j)         Derivative of ith local coord wrt jth global coord
#      det                Determinant of jacobian
#      strain(i,j)        strain_ij components
#      dsde(i,j,k,l)      Derivative of stress_ij with respect:strain_kl
#      kel(row,col)       Rows && cols of element stiffness

    npoints = numberofintegrationpoints(ncoord, nelnodes, elident)
    dNdx = np.zeros((nelnodes, ncoord))
    dxdxi = np.zeros((ncoord, ncoord))
    strain = np.zeros((ndof, ncoord))
    kel = np.zeros((ndof*nelnodes, ndof*nelnodes))

    # Set up integration points && weights 

    xilist = integrationpoints(ncoord,nelnodes,npoints,elident)
    w = integrationweights(ncoord,nelnodes,npoints,elident)

    # Loop over the integration points

    for intpt in range(npoints):

        # Compute shape functions && derivatives wrt local coords

        xi = np.zeros((ncoord,1))
        for i in range(0,ncoord):
          xi[i] = xilist[i,intpt] 

        N = shapefunctions(nelnodes, ncoord, elident, xi)
        dNdxi = shapefunctionderivs(nelnodes, ncoord, elident, xi)

        # Compute the jacobian matrix && its determinant

        for i in range(ncoord):
            for j in range(ncoord):
                dxdxi[i][j] = 0
                for a in range(nelnodes):
                    dxdxi[i][j] += coord[i][a] * dNdxi[a][j]
        
        dxidx = np.linalg.inv(dxdxi)
        dt = np.linalg.det(dxdxi)

        # Convert shape function derivatives:derivatives wrt global coords
        
        for a in range(nelnodes):
            for i in range(ncoord):
                dNdx[a, i] = 0
                for j in range(ncoord):
                    dNdx[a, i] += dNdxi[a, j] * dxidx[j, i]

        # Compute the deformation gradients by differentiating displacements

        F = np.eye(ncoord)
        for i in range(ncoord):
            for j in range(ncoord):
                for a in range(nelnodes):
                    F[i][j] += displacement[i][a] * dNdx[a][j]
        
        # Compute Bbar and J

        J = np.linalg.det(F)
        B = np.dot(F, F.T)

        # Convert shape function derivatives to derivatives wrt spatial coords

        Finv = np.linalg.inv(F)
        dNdxs = np.zeros((nelnodes, ncoord))
        for a in range(nelnodes):
            for i in range(ncoord):
                for j in range(ncoord):
                    dNdxs[a, i] += dNdx[a, j] * Finv[j, i]

        # Compute the stress

        stress = Kirchhoffstress(ndof, ncoord, B, J, materialprops)

        # Compute the material tangent stiffness (d stress/d strain)
        # ds/de is just C_ijkl for linear elasticity - this notation is used
        # to allow extension to nonlinear problems

        dsde = materialstiffness(ndof, ncoord, B, J, materialprops)
        
        # Compute the element stiffness

        for a in range(nelnodes):
            for i in range(ndof):
                for b in range(nelnodes):
                    for k in range(ndof):
                        row = ndof*(a-1)+i
                        col = ndof*(b-1)+k
                        for j in range(ncoord):
                            for l in range(ncoord):
                                kel[row][col] += dsde[i][j][k][l] * dNdxs[b][l] * dNdxs[a][j] * w[intpt] * dt
                        kel[row][col] -= stress[i][j] * dNdxs[a][k] * dNdxs[b][j] * w[intpt] * dt
      
    return kel

######################################################################################################################################################

#====================== No. nodes on element faces ================
#
#   This procedure returns the number of nodes on each element face
#   for various element types.  This info is needed for computing
#   the surface integrals associated with the element traction vector

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

######################################################################################################################################################

#======================= Lists of nodes on element faces =============
#
#    This procedure returns the list of nodes on an element face
#    The nodes are ordered so that the element face forms either
#    a 1D line element or a 2D surface element for 2D or 3D problems
#
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

######################################################################################################################################################

#
#====================== ELEMENT DISTRIBUTED LOAD VECTOR ==============
#

def eldload(ncoord, ndof, nfacenodes, elident, coords, traction):
    npoints = numberofintegrationpoints(ncoord-1, nfacenodes, elident)
    xi =  np.zeros((ncoord-1,1))
    dxdxi = np.zeros((ncoord,ncoord-1))
    r = np.zeros((ndof*nfacenodes,1))
   
    xilist = integrationpoints(ncoord-1, nfacenodes, npoints, elident)
    w = integrationweights(ncoord-1, nfacenodes, npoints, elident)

    for intpt in range(npoints):
        for i in range(ncoord-1):
            xi[i] = xilist[i][intpt]
        N = shapefunctions(nfacenodes, ncoord-1, elident, xi)
        dNdxi = shapefunctionderivs(nfacenodes, ncoord-1, elident, xi)

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

######################################################################################################################################################
#
#====================== Assemble the global residual vector =================
#

def globalresidual(ncoord, ndof, nnode, coords, nelem, maxnodes, elident, nelnodes, connect, materialprops, dofs):
    
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

######################################################################################################################################################

#
#====================== Assemble the global stiffness matrix =================
#

def globalstiffness(ncoord, ndof, nnode, coords, nelem, maxnodes, elident, nelnodes, connect, materialprops, dofs):
    
    # Assemble the global stiffness matrix

    Stif = np.zeros((ndof*nnode,ndof*nnode))
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
        kel = elstif(ncoord, ndof, n, ident, lmncoord, materialprops, lmndof)
        
        # Add the current element stiffness:the global stiffness

        for a in range(nelnodes):
            for i in range(ndof):
                for b in range(nelnodes):
                    for k in range(ndof):
                        rw = ndof*(connect[a][lmn]-1)+i
                        cl = ndof*(connect[b][lmn]-1)+k
                        Stif[rw][cl] += kel[ndof*(a-1)+i][ndof*(b-1)+k]
    return Stif

######################################################################################################################################################

#
#===================== Assemble the global traction vector =============
#

def globaltraction(ncoord, ndof, nnodes, ndload, coords, nelnodes, elident, connect, dloads, dofs):
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

######################################################################################################################################################


def print_results(outfile, nprops, materialprops, ncoord, ndof, nnode, coords, nelem, maxnodes, connect, nelnodes, elident, nfix, fixnodes, ndload, dloads, dofs):

    outfile.write("Nodal Displacements: \n")
    if ndof == 2:
        outfile.write(" Node Coords u1 u2 \n")
        for i in range(nnode):
            outfile.write("{:3d} {:8.4f} {:8.4f} {:8.4f} {:8.4f}\n".format(
                i, coords[0, i], coords[1, i], dofs[2*i][0], dofs[2*i+1][0]
            ))

    outfile.write('\n\n Strains and Stresses \n')

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

        npoints = numberofintegrationpoints(ncoord,n, elident)
        dNdx = np.zeros((n,ncoord))
        dxdxi = np.zeros((ncoord,ncoord))
        xi = np.zeros((ncoord,1))
        x = np.zeros((ncoord,1))
        F = np.zeros((ncoord,ncoord))

        # Set up integration points
        xilist = integrationpoints(ncoord,n,npoints, elident)

        for intpt in range(npoints):

            for i in range(ncoord):
                xi[i] = xilist[i][intpt]

            N = shapefunctions(nelnodes, ncoord, elident, xi)
            dNdxi = shapefunctionderivs(nelnodes, ncoord, elident, xi)

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

            stress = Kirchhoffstress(ndof,ncoord,B,J,materialprops)
            stress = stress/J

            if ncoord == 2:
                # print("{:5d} {:7.4f} {:7.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f}".format(intpt, x[0][0], x[1][0], B[0][0], B[1][1], B[0][1], stress[0][0], stress[1][1], stress[0][1]), file=outfile)
                print("{:5d} {:7.4f} {:7.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f}".format(intpt, x[0][0], x[1][0], B[0][0], B[1][1], B[0][1], stress[0][0], stress[1][1], stress[0][1], F[0][0], F[1][0], F[1][0], F[1][1] ), file=outfile)



# def plotmesh(coords, ncoord, nnode, connect, nelem, elident, nelnodes, color):
#     f2D_3 = [1, 2, 3]
#     f2D_4 = [1, 2, 3, 4]
#     f2D_6 = [1, 4, 2, 5, 3, 6]
#     f2D_8 = [1, 5, 2, 6, 3, 7, 4, 8]
#     plt.hold(True)
#     if ncoord == 2:  # Plot a 2D mesh
#         for lmn in range(nelem):
#             x = []
#             for i in range(nelnodes[lmn]):
#                 x.append([coords[0][connect[i][lmn]], coords[1][connect[i][lmn]]])
#             plt.scatter([point[0] for point in x], [point[1] for point in x], marker='o', color='r')
#             if nelnodes[lmn] == 3:
#                 poly = Polygon(x, closed=True, fill=False, edgecolor=color)
#                 plt.gca().add_patch(poly)
#             elif nelnodes[lmn] == 4:
#                 poly = Polygon(x, closed=True, fill=False, edgecolor=color)
#                 plt.gca().add_patch(poly)
#             elif nelnodes[lmn] == 6:
#                 poly = Polygon(x, closed=True, fill=False, edgecolor=color)
#                 plt.gca().add_patch(poly)
#             elif nelnodes[lmn] == 8 or nelnodes[lmn] == 9:
#                 poly = Polygon(x, closed=True, fill=False, edgecolor=color)
#                 plt.gca().add_patch(poly)
#     plt.axis('equal')
#     plt.hold(False)

#######


# infile = "hyperelastic_quad4.txt"
# infile = open('hyperelastic_quad4.txt', 'r')
outfile = open('FEM_results.txt', 'w')

# nprops, materialprops, ncoord, ndof, nnode, coords, nelem, maxnodes, connect, nelnodes, elident, nfix, fixnodes, ndload, dloads = read_input_file("hyperelastic_quad4.txt")

# infile.close()
# plt.figure()
# plt.triplot(coords[:,0], coords[:,1], connect)
# plt.show()


##########################################################################################################################################


#============================ MAIN FEM ANALYSIS PROCEDURE ========================
#
#   w           Nodal displacements.  Let w_i^a be ith displacement component
#               at jth node.  Then dofs contain (w_1^1, w_2^1, w_1^2, w_2^2....) for 2D
#               and (w_1^1, w_2^1, w_3^1, w_1^2, w_2^2, w_3^2....) for 3D
#   dw          Correction to nodal displacements.  Let w_i^a be ith displacement component
#               at jth node.  Then dofs contain (w_1^1, w_2^1, w_1^2, w_2^2....) for 2D
#               and (w_1^1, w_2^1, w_3^1, w_1^2, w_2^2, w_3^2....) for 3D
#   K           Global stiffness matrix.  Stored as [K_1111 K_1112  K_1121  K_1122...
#                                                    K_1211 K_1212  K_1221  K_1222...
#                                                    K_2111 K_2112  K_2121  K_2122...]
#               for 2D problem and similarly for 3D problem
#   F           Force vector.  Currently only includes contribution from tractions
#               acting on element faces (i.e. body forces are neglected)
#   R           Volume contribution to residual
#   b           RHS of equation system


w = np.zeros((nnode*ndof,1))

#  Here we specify how the Newton Raphson iteration should run
#  Load is applied in nsteps increments;
#  tol is the tolerance used in checking Newton-Raphson convergence
#  maxit is the max no. Newton-Raphson iterations
#  relax is the relaxation factor (Set to 1 unless big convergence problems)

nsteps = 5
tol = 0.0001
maxit = 30
relax = 1.0

forcevdisp = np.zeros((2,nsteps+1))
forcevdisp[0,0] = 0
forcevdisp[1,0] = 0

for step in range(1, nsteps+1):
    loadfactor = step/nsteps

    err1 = 1.0
    nit = 0

    print(f'\n Step {step} Load {loadfactor}\n')

    while ((err1>tol) and (nit<maxit)): # Newton Raphson loop
        nit = nit + 1
   
        K = globalstiffness(ncoord,ndof,nnode,coords, 
                nelem,maxnodes,elident,nelnodes,connect,materialprops,w)
        F = globaltraction(ncoord,ndof,nnode,ndload,coords, 
                    nelnodes,elident,connect,dloads,w)
        R = globalresidual(ncoord,ndof,nnode,coords, 
                nelem,maxnodes,elident,nelnodes, 
                connect,materialprops,w)

        b = loadfactor*F - R
         
        # Fix constrained nodes.
        for n in range(nfix):
            rw = int((ndof*(fixnodes[0][n]-1) + fixnodes[1][n])-1)
            for cl in range(ndof*nnode):
                K[rw][cl] = 0
            K[rw][rw] = 1.
            b[rw] = loadfactor*fixnodes[2][n] - w[rw]
        
        # Solve for the correction
        dw = np.linalg.solve(K,b)
        
        # Check convergence
        w = w + relax*dw
        wnorm = float(np.dot(w.T,w))
        err1 = float(np.dot(dw.T,dw))
        err2 = float(np.dot(b.T,b))
        err1 = math.sqrt(err1/wnorm)
        err2 = math.sqrt(err2)/(ndof*nnode)

        print(f'Iteration number {nit} Correction {err1} Residual {err2} tolerance {tol}')
    
    outfile.write(f'\n\n Step {step} Load {loadfactor}\n')

    print_results(outfile, 
        nprops,materialprops,ncoord,ndof,nnode,coords, 
        nelem,maxnodes,connect,nelnodes,elident, 
        nfix,fixnodes,ndload,dloads,w)
    
    # Store traction and displacement for plotting later
    forcevdisp[1,step] = loadfactor*dloads[2,0]
    forcevdisp[0,step] = w[2][0]

plt.plot(forcevdisp[0,:], forcevdisp[1,:], 'r', linewidth=3)
plt.xlabel('Displacement', fontsize=16)
plt.ylabel('Force', fontsize=16)
plt.show()


#================================= POST-PROCESSING =================================
#
# Create a plot of the deformed mesh


defcoords = np.zeros((ndof,nnode))
scalefactor = 1.0
for i in range(nnode):
    for j in range(ndof):
        defcoords[j,i] = coords[j,i] + scalefactor*w[ndof*(i-1)+j]

# title_string ="Hyperelastic"

# plt.scatter(coords[0],coords[1],s=20,color='m',marker='X',label='Undeformed coordinates')
# plt.scatter(defcoords[0],defcoords[1],s=10,color='r',label='Deformed coordinates')
# plt.title(title_string)
# plt.legend()
# plt.show()

                    
# plt.figure()
# # plotmesh(coords,ncoord,nnode,connect,nelem,elident,nelnodes,'g')
# plt.hold(True)
# # plotmesh(defcoords,ncoord,nnode,connect,nelem,elident,nelnodes,'r')

# plt.close(outfile)
