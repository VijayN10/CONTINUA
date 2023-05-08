import numpy as np


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
materialprops = np.array([[1], [10]])                                        # materialprops = np.array([[2], [1]])    # materialprops[mu1,K1]      # Python first position is 0 --> change in all code
                                        # first row implies multiscale, Second row implies type of RVE
                                    
# mu1 = 1
# k1 = 10
model = materialprops[0][0]

if model == 1 :
    mul = materialprops[0][0]
    k1 = materialprops[1][0]
elif model == 2:
    rvetype = materialprops[1][0]

# no. coords (1:3), no. DOF, no. nodes and nodal coordinates

ncoord = 2
ndof = 2
nnode = 9

coords = np.array([[0, 1, 2, 0, 1, 2, 0, 1, 2],
                   [0, 0, 0, 1, 1, 1, 2, 2, 2]])

# No. elements and connectivity

nelem = 4
maxnodes = 4 
nelnodes = 4   
elident = np.array([[1],[2],[3],[4]])                 # elident = np.array([[1]])

connect = np.array([[1, 2, 5, 4],
                    [2, 3, 6, 5],
                    [5, 6, 9, 8],
                    [4, 5, 8, 7]])     # connect = np.array([[1],[2],[3],[4]])

# No. nodes with prescribed displacements, with the prescribed displacements

nfix = 6
fixnodes = np.array([[1, 1, 2, 3, 4, 7],
                     [1, 2, 2, 2, 1, 1],
                     [0, 0, 0, 0, 0, 0]])

# No. loaded element faces, with the loads

ndload = 2
dloads = np.array([[2, 3],
                   [2, 2],
                   [3, 3],
                   [0, 0]])      # dloads = np.array([[1],[2],[3],[0]])


# Name for the Giraffe input file (without identification number)

name = "tex"