#%%
from numba import jit
import numpy as np
from numpy.linalg import eig, inv
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import math
import subprocess
import time
import os
import shutil
import sys


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


def meshGen(nelemx, nelemy, width, height, maxnodes, nelnodes):

    """
    Generate a structured mesh of quadrilateral elements.

    Args:
        nelemx: number of elements along the x-axis
        nelemy: number of elements along the y-axis
        width: total width of the mesh
        height: total height of the mesh
        maxnodes: maximum number of nodes per element
        nelnodes: number of nodes per element

    Returns:
        coords: node coordinates (2D array, shape (2, nnodes))
        connect: element connectivity (2D array, shape (nelem, maxnodes))
        elident: element identifiers (1D array, shape (nelem,))
        nnode: number of nodes
        nelem: number of elements
        maxnodes: maximum number of nodes per element
        nelnodes: number of nodes per element
    """

    # Compute node coordinates
    dx = width / nelemx
    dy = height / nelemy
    x = np.linspace(0, width, nelemx+1)   # x-coordinates of nodes
    y = np.linspace(0, height, nelemy+1)  # y-coordinates of nodes
    X, Y = np.meshgrid(x, y) 
    coords = np.vstack([X.ravel(), Y.ravel()])

    # Compute element connectivity
    
    nnode = (nelemx+1) * (nelemy+1)
    nelem = nelemx * nelemy
    
    # Check if only one element is present
    if nelem == 1:
        proceed = input("Meshing not done: only one element present. Do you want to proceed? (y/n)")
        if proceed != "y":
            sys.exit()
            
    elident = np.arange(1, nelem+1).reshape(-1, 1)
    connect = np.zeros((nelem, maxnodes), dtype=int)
    for j in range(nelemy):
        for i in range(nelemx):
            n1 = j*(nelemx+1) + i  # node index of lower left corner of element
            n2 = n1 + 1    # node index of lower right corner of element
            n3 = n2 + nelemx+1  # node index of upper right corner of element
            n4 = n1 + nelemx+1  # node index of upper left corner of element
            if j % 2 == 0:
                e = j*nelemx + i  # for even-numbered rows, elements are numbered left to right
                connect[e] = [n1+1, n2+1, n3+1, n4+1]
            else:
                e = (j+1)*nelemx - i - 1   # for odd-numbered rows, elements are numbered right to left
                connect[e] = [n1+1, n2+1, n3+1, n4+1]
    
    connect = connect.T

    # Return mesh information
    return coords, connect, elident, nnode, nelem, maxnodes, nelnodes


nelemx = 2
nelemy = 2
width = 1
height = 1
maxnodes = 4
nelnodes = 4
coords, connect, elident, nnode, nelem, maxnodes, nelnodes = meshGen(nelemx, nelemy, width, height, maxnodes, nelnodes)



left_nodes = []  # List to store the node numbers on the extreme left side of the mesh
for i in range(nelemy + 1):
    left_nodes.append(i * (nelemx + 1) + 1)
# print(left_nodes)

bottom_nodes = []  # List to store the node numbers on the extreme bottom of the mesh
for i in range(nelemx + 1):
    bottom_nodes.append(i + 1)
# print(bottom_nodes)


nfix = len(bottom_nodes) + len(left_nodes)
fixnodes = np.zeros((3, nfix), dtype=int)

# Adding bottom nodes to the first row
fixnodes[0, :len(bottom_nodes)] = bottom_nodes

# Adding left nodes to the first row
fixnodes[0, len(bottom_nodes):] = left_nodes

# Sorting the first row in ascending order
fixnodes[0, :] = np.sort(fixnodes[0, :])


# Determine the values in the second row based on the position of the nodes
for i in range(nfix):
    node = fixnodes[0, i]
    if node in bottom_nodes and node in left_nodes:
        # Node is at the corner of the grid
        fixnodes[1, i] = 0
    elif node in bottom_nodes:
        # Node is on the bottom row
        fixnodes[1, i] = 2
    elif node in left_nodes:
        # Node is on the left column
        fixnodes[1, i] = 1

# Set the first two elements of the second row to 1 and 2 respectively
fixnodes[1, 0] = 1
fixnodes[1, 1] = 2

# Third row will be all zeros
fixnodes[2, :] = 0

# print(fixnodes)


# Find the elements on the rightmost side
right_elems = []
for j in range(nelemy):
    if j % 2 == 0:
        e = (j + 1) * nelemx - 1
        right_elems.append(e + 1)
    else:
        e = j * nelemx
        right_elems.append(e + 1)

# print(right_elems)


# print(right_elems)
ndload = len(right_elems)
# print(ndload)
dloads = np.zeros((4, ndload))
dloads[0] = right_elems
dloads[1] = 2
dloads[2] = 3
dloads[3] = 0

# print(dloads)


# Plot nodes with labels
plt.scatter(coords[0], coords[1], s=10, color='black')
for i in range(nnode):
    plt.text(coords[0,i]+0.05, coords[1,i]+0.05, i+1, color='red', fontsize=12)

    
plt.xlabel('x')
plt.ylabel('y')
plt.title('Structured Mesh')
plt.axis('equal')
plt.show()




# nfix = 8
# fixnodes = np.array([[1, 1, 2, 3, 4, 5, 9, 13],
#                      [1, 2, 2, 2, 2, 1, 1, 1],
#                      [0, 0, 0, 0, 0, 0, 0, 0]])

# # No. loaded element faces, with the loads

# ndload = 3
# dloads = np.array([[3, 4, 9],
#                    [2, 2, 2],
#                    [3, 3, 3],
#                    [0, 0, 0]])      # dloads = np.array([[1],[2],[3],[0]])


# Name for the Giraffe input file (without identification number)

name = "tex"

######################################################################################################################################################


if model ==2:
    with open("Giraffe/RVE_data/grfTop" + str(rvetype) + ".txt", "r") as file:
        top = file.read()

    with open("Giraffe/RVE_data/grfBottom" + str(rvetype) + ".txt", "r") as file:
        bottom = file.read()

    # Read RVE data file
    with open("Giraffe/RVE_data/rve{}.txt".format(rvetype), "r") as file:
        data = [float(x) for x in file.read().split()]
        Lxx, Lyy, t, nx, ny = data[:5]
        nx = int(nx)
        ny = int(ny)

    # Get the nodal coordinates of the ends of the beams

    X = np.zeros(((nx + ny) * 2, 3))

    # Read the nodedata text files from data folder
    with open("Giraffe/RVE_data/nodedata" + str(rvetype)  + ".txt", "r") as file:
        for i, line in enumerate(file):
            X[i] = [float(x) for x in line.strip().split()]



######################################################################################################################################################


# To create Giraffe input file 
# It creates new folder (with the same name of Giaffe input file) at location of executables. 
# The folder contains the created file.

def giraffeInputGenerator(rvetype, name, F, nx, ny, X, step, loadfactor, nit, intpt, aa):

    # Consider rvetype                 
    rr = rvetype  


    # Create new folder with name of .inp file
    inp = name + '_' + str(rr) + '_step' + str(step) + '_load' + str(loadfactor) + '_iteration' + str(nit) +'_intgPt' + str(intpt) + '_aa'  + str(aa)               # inp = tex_0                     rvetype_step_load_iteration_integrationPoint_aa
    folder = "Giraffe/" + inp                              # folder = E:/Softwares/01/Giraffe/tex_0
    os.makedirs("Giraffe/" + inp, exist_ok=True)               # Creates tex_0 folder

    
    # Deformation gradient F
    # F = np.array([[1.2, 0.2],   # Assumption used for testing this function
    #             [0.2, 1.2]])


    # Get the strain from the deformation gradient
    strain = F - np.eye(2)

    # Get alpha values from strain tensor
    alpha1 = strain[0,0]
    alpha2 = strain[0,1]
    alpha3 = strain[1,0]
    alpha4 = strain[1,1]
    
    # Convert alpha into a matrix
    alpha = np.array([[alpha1, alpha2],[alpha3, alpha4]])
    # print(f'alpha = {alpha}')


    # Initialize disp and calculate the displacement
    # of the ends of the beams
    disp = np.zeros(((nx + ny) * 2, 2))
    disp[:,0] = alpha1*X[:,0] + alpha2*X[:,1]
    disp[:,1] = alpha3*X[:,0] + alpha4*X[:,1]

    

    displacement_block = []

    for i in range(1, ((nx + ny) * 2) + 1):
        block = ("""\n
NodalDisplacement {} NodeSet {} CS 125 NTimes 2
//Time UX UY UZ ROTX ROTY ROTZ
0\t0\t0\t0\t0\t0\t0
1\t{}\t{}\t0\t0\t0\t0\n
""").format(i, i + 4, disp[i-1][0], disp[i-1][1])
        displacement_block.append(block)

  
    # Open the file for writing

    filepath = 'Giraffe/' + inp + '/' + inp + '.inp'         # filepath - E:/Softwares/01/Giraffe/tex_0/tex_0.inp
    with open(filepath, 'w') as file:
        
    # Top part of Giraffe input file

        file.write(top) 

    # Displacement block

        file.write("\n\nDisplacements {}\n\n".format((nx + ny) * 2))
        for i in range((nx + ny) * 2):
            file.write(displacement_block[i])


    # Bottom part of Giraffe input file

        file.write(bottom)


    return inp, folder

######################################################################################################################################################

# Run giraffe
#@jit(nopython=False, parallel=True)
def runGiraffe(inp):

    # path : E:/Softwares/01/Giraffe/
    # folder : E:/Softwares/01/Giraffe/tex_0
    # inp  : tex_0

    # path of Giraffe.exe
    # giraffe = path + "Giraffe.exe"
    
    # cur_folder = os.path.dirname(os.getcwd())
    # print(cur_folder)
    
    cur_folder = os.getcwd()
    # print(cur_folder)    
   
    # Changing directory
    os.chdir('Giraffe')

    p = subprocess.Popen(["Giraffe.exe"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    # p = subprocess.Popen([r"E:\Softwares\01\Giraffe\Giraffe.exe"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    input_file_name = inp                                                            # input("Enter the name of the input file: ")
    # p.communicate(input=input_file_name.encode())
    
    output, error = p.communicate(input=input_file_name.encode())
    if error:
        print("Error: ", error.decode())
    else:
        print("Output: ", output.decode())


    # # report = r"E:\Softwares\01\Giraffe\tex_0\simulation_report.txt"
    # report = 'Giraffe/' + inp + '/simulation_report.txt'

    # while True:
    #     with open(report, "r") as f:
    #         lines = f.readlines()
    #         if lines:
    #             last_line = lines[-1]
    #             if "Total solution time:" in last_line:
    #                 print("Simulation completed.\n", last_line)
    #                 break
    #             else:
    #                 print(last_line, end="")
    #         time.sleep(1)   
    
    # Going back to current directory or folder
    os.chdir(cur_folder)
    
    opFilePath = 'Giraffe/' + inp + '/monitors/monitor_nodeset_1.txt'

    return opFilePath


######################################################################################################################################################

# Check the size of output file

# The function returns a boolean value indicating whether the file is empty or not. 
# If the file is empty, os.path.getsize('.txt') will return 0 and checkGiraffeOutputs will return True. 
# If the file is non-empty, os.path.getsize('.txt') will return a value greater than 0 and checkGiraffeOutputs will return False.
#@jit(nopython=False, parallel=True)
def checkGiraffeOutputs(opFilePath):

    if os.path.getsize(opFilePath) == 0:
        print("File is empty! We are running the Giraffe again...")
        return True
    
    return False


######################################################################################################################################################

# def giraffeStress(opFilePath):

#@jit(nopython=False, parallel=True)
def giraffeStress(Lxx, Lyy, t, inp):
      
    # # Changing t temporarily    # Ask
    # t = 1
    
    # File for left side
    # monfilepath = "4x4Job/monitors/monitor_nodeset_2.txt"
    monfilepath = 'Giraffe/' + inp + '/monitors/monitor_nodeset_2.txt'
    data_left = np.genfromtxt(monfilepath, skip_header=1)
    posx = np.abs(data_left[:,1])
    posx = posx - posx[0]
    posy = np.abs(data_left[:,2])
    posy = posy - posy[0]
    forcex = np.abs(data_left[:,7])
    forcey = np.abs(data_left[:,8])
    data_left = np.column_stack((posx, posy, forcex, forcey))

    # File for top side
    # monfilepath = "4x4Job/monitors/monitor_nodeset_4.txt"
    monfilepath = 'Giraffe/' + inp + '/monitors/monitor_nodeset_4.txt'
    data_top = np.genfromtxt(monfilepath, skip_header=1)
    posx = np.abs(data_top[:,1])
    posx = posx - posx[0]
    posy = np.abs(data_top[:,2])
    posy = posy - posy[0]
    forcex = np.abs(data_top[:,7])
    forcey = np.abs(data_top[:,8])
    data_top = np.column_stack((posx, posy, forcex, forcey))

    # Get the forces
    Fxx = data_left[-1,2]
    Fyy = data_top[-1,3]
    Fxy = data_left[-1,3]
    Fyx = data_top[-1,2]

    # Stress
    stress = np.zeros((2,2))
    stress[0,0] = Fxx/(1000*Lyy*t)
    stress[1,1] = Fyy/(1000*Lxx*t)
    stress[0,1] = Fxy/(1000*Lyy*t)
    stress[1,0] = Fyx/(1000*Lxx*t)

    return stress

######################################################################################################################################################

# def giraffeStiffness(opFilePath)

######################################################################################################################################################

#@jit(nopython=False, parallel=True)
def Voigt2normal(Dnumerical):
    
    dsde = np.zeros((2,2,2,2))
    
    dsde[0,0,0,0] = Dnumerical[0,0]
    dsde[0,0,1,1] = Dnumerical[0,1]
    dsde[0,0,0,1] = Dnumerical[0,2]
    dsde[0,0,1,0] = Dnumerical[0,2]

    dsde[1,1,0,0] = Dnumerical[1,0]
    dsde[1,1,1,1] = Dnumerical[1,1]
    dsde[1,1,0,1] = Dnumerical[1,2]
    dsde[1,1,1,0] = Dnumerical[1,2]

    dsde[0,1,0,0] = Dnumerical[2,0]
    dsde[0,1,1,1] = Dnumerical[2,1]
    dsde[0,1,0,1] = Dnumerical[2,2]
    dsde[0,1,1,0] = Dnumerical[2,2]

    dsde[1,0,0,0] = Dnumerical[2,0]
    dsde[1,0,1,1] = Dnumerical[2,1]
    dsde[1,0,0,1] = Dnumerical[2,2]
    dsde[1,0,1,0] = Dnumerical[2,2]

    return dsde


#%% 
######################################################################################################################################################
# AREA TO TEST CODE 21/02/2023

#%% 
######################################################################################################################################################

#================= Material Stiffness ==================================
#
#    Computes material stiffness tensor C_{ijkl} 
#    Currently coded either for plane strain or general 3D.
#@jit(nopython=False, parallel=True)
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
#@jit(nopython=False, parallel=True)
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
#@jit(nopython=False, parallel=True)
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
    return n

######################################################################################################################################################

#====================== INTEGRATION POINTS ==================================
#
#   Defines positions of integration points
#@jit(nopython=False, parallel=True)
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
#@jit(nopython=False, parallel=True)
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

    return w

######################################################################################################################################################

#================= SHAPE FUNCTIONS ==================================
#
#        Calculates shape functions for various element types
#@jit(nopython=False, parallel=True)
def shapefunctions(nelnodes, ncoord, elident, xi):
    
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

######################################################################################################################################################

#
#================= SHAPE FUNCTION DERIVATIVES ======================
#
#@jit(nopython=False, parallel=True)
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
#@jit(nopython=False, parallel=True)
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
#@jit(nopython=False, parallel=True)
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
        
                    
        with open('F_values.txt', 'a') as file:
            for row in F:
                for val in row:
                    file.write(f'{val:.6f} ')
            file.write('\n')
        
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

        if model == 1:

            # Compute the stress

            stress = Kirchhoffstress(ndof, ncoord, B, J, materialprops)

            # Compute the material tangent stiffness (d stress/d strain)
            # ds/de is just C_ijkl for linear elasticity - this notation is used
            # to allow extension to nonlinear problems

            dsde = materialstiffness(ndof, ncoord, B, J, materialprops)
        

        elif model == 2:
            
            aa = 10
            # Create a new Giraffe file
            inp, folder = giraffeInputGenerator(rvetype, name, F, nx, ny, X, step, loadfactor, nit, intpt, aa)

            # Run Giraffe
            print(f"Giraffe is running to get stess for intpt = {intpt}...")
            opFilePath = runGiraffe(inp)

            # Check if Giraffe has run

            # Check if the output file is empty using `checkGiraffeOutputs` function
            flag_giraffe = checkGiraffeOutputs(opFilePath)  # size of text file

            # Set a counter to limit the number of tries
            counterGiraffe = 0

            # Keep looping while the output file is empty and counter is less than 10
            while flag_giraffe and counterGiraffe < 10:

                # Run the `runGiraffe` function and update the `opFilePath' or are are we running the same file?
                opFilePath = runGiraffe(inp)

                # Check the output file again
                flag_giraffe = checkGiraffeOutputs(opFilePath)

                # Increment the counter
                counterGiraffe += 1

            # Check if flag_giraffe is False
            if not flag_giraffe:
                print("Calculating stress and dsde")
                            
            else:
                # If flag_giraffe is True, print message and exit program
                print("File is still empty after 10 tries.")
                exit()   

            # Get the stress from Giraffe outputs
            stress = giraffeStress(Lxx, Lyy, t, inp)


            # # print("Task completed and now deleting folder...")
            # shutil.rmtree(folder)        # Deleting the folder         
            # # print("folder deleted...")
            
            
            ## STEP 01

            # # Convert to the right stress measures          # WIP
            # stress = np.dot(F, stress.T)                 

            # Get the Right Cauchy Green deformation tensor
            C = np.dot(F.T, F)

            # Extract the eigenvalues and eigenvectors of C
            # These are the principle stretches and the principle directions
            val, vec = eig(C)     # check val, vec or vec, val

            # Initialize and get the stretch tensor by decomposition
            U = np.zeros((2,2))
            for i in range(2):
                U += np.sqrt(val[i])*np.outer(vec[:,i], vec[:,i])

            # Get rotation matrix using F = RU decomposition
            Rot = np.dot(F, inv(U))

            ## STEP 02

            # Machine precision
            sqeps = np.sqrt(np.finfo(float).eps)
                
            # Initialize a matrix for the Numerical tangent
            Dnumerical = np.zeros((3,3))
            io = [0, 1, 0, 1]
            it = [0, 1, 1, 0]

            ## STEP 03

            # Start the loop to calculate the numerical tangent

            # Calculate the numerical tangent
            for aa in range(3):
                    
                # Initialize to zero
                dC = np.zeros((2,2))
                        
                # If clause
                if(aa < 2):
                    if(C[aa,aa] == 0):
                        e = sqeps
                    else:
                        e = sqeps*C[aa,aa]
                    dC[aa,aa] = 1
                else:
                    if(C[io[aa],it[aa]] == 0):
                        e = sqeps
                    else:
                        e = sqeps*C[io[aa],it[aa]]
                    dC[io[aa],it[aa]] = 0.5
                    dC[it[aa],io[aa]] = 0.5

                # Calculate the pertured Left Cauchy green deformation tensor (LCGDT)
                Cp = C + e*dC
                
                # Calcualte the perturbed Related perturbed kinematic measures
                # Get eigenvalues and vector for perturbed LCGDT
                valp, vecp = eig(Cp)

                # Get perturbed stretch tensor
                Up = np.zeros((2,2))
                for kk in range(2):
                    Up += np.sqrt(valp[kk])*np.outer(vecp[:,kk], vecp[:,kk])

                # Get perturbed deformation gradient
                Fp = np.dot(Rot, Up)
                
                # F = Fp  # TBC

                # Create a new Giraffe file
                inp, folder = giraffeInputGenerator(rvetype, name, Fp, nx, ny, X, step, loadfactor, nit, intpt, aa)
                # inp, Lxx, Lyy, t, folder = giraffeInputGenerator(rvetype, name, F, Lxx, Lyy, t, nx, ny, X)  # Send Fp here (Need to add other function?)

                # Run Giraffe
                print(f"Giraffe is running to get stresspm for intpt = {intpt} and aa = {aa}...")
                opFilePath = runGiraffe(inp)

                # Check if Giraffe has run
                
                # Check if the output file is empty using `checkGiraffeOutputs` function
                flag_giraffe = checkGiraffeOutputs(opFilePath)  # size of text file

                # Set a counter to limit the number of tries
                counterGiraffe = 0

                # Keep looping while the output file is empty and counter is less than 10
                while flag_giraffe and counterGiraffe < 10:

                    # Run the `runGiraffe` function and update the `opFilePath' or are are we running the same file?
                    opFilePath = runGiraffe(inp)

                    # Check the output file again
                    flag_giraffe = checkGiraffeOutputs(opFilePath)

                    # Increment the counter
                    counterGiraffe += 1

                # Check if flag_giraffe is False
                if not flag_giraffe:
                    print("Calculating stress and dsde")
                            
                else:
                    # If flag_giraffe is True, print message and exit program
                    print("File is still empty after 10 tries.")
                    exit()   

                # Get the stress from Giraffe outputs
                stresspm = giraffeStress(Lxx, Lyy, t, inp)

                # Convert the stress measure
                stresspm = np.dot(Fp, stresspm.T)                           
                    
                # Get the numerical tangent
                Dnumerical[0,aa] = (2/e)*(stresspm[0,0] - stress[0,0])
                Dnumerical[1,aa] = (2/e)*(stresspm[1,1] - stress[1,1])
                Dnumerical[2,aa] = (2/e)*(stresspm[0,1] - stress[0,1])

                # # Delete the Giraffe case folder
                # # print("Task completed and now deleting folder...")
                # shutil.rmtree(folder)                                                # WIP Not removing folder
                # # print("folder deleted...")

            # Convert the obtained numerical tangent to the form
            # usable by this code 
            dsde = Voigt2normal(Dnumerical) 

            # print(dsde)

   
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

######################################################################################################################################################

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
            fnlist[1] = i4[int(face)-1]
        elif (nelnodes==8): 
            fnlist[0] = face
            fnlist[1] = i4[face-1]
            fnlist[2] = face+4
    return fnlist

######################################################################################################################################################

#
#====================== ELEMENT DISTRIBUTED LOAD VECTOR ==============
#
#@jit(nopython=False, parallel=True)
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
#@jit(nopython=False, parallel=True)
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
#@jit(nopython=False, parallel=True)
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
#@jit(nopython=False, parallel=True)
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

#@jit(nopython=False, parallel=True)
def print_results(outfile, nprops, materialprops, ncoord, ndof, nnode, coords, nelem, maxnodes, connect, nelnodes, elident, nfix, fixnodes, ndload, dloads, dofs):

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


# Plot Mesh
def plotmesh(coords, ncoord, nnode, connect, nelem, elident, nelnodes, color):
    f2D_3 = np.array([1,2,3])
    f2D_4 = np.array([1,2,3,4])
    f2D_6 = np.array([1,4,2,5,3,6])
    f2D_8 = np.array([1,5,2,6,3,7,4,8])

    fig, ax = plt.subplots()
    for lmn in range(nelem):
        x = []
        for i in range(nelnodes):
            x.append(coords[:, connect[i, lmn]-1])
        x = np.array(x)
        ax.scatter(x[:, 0], x[:, 1], color='r', edgecolors='none')
        if nelnodes == 3:
            ax.add_patch(Polygon(x[f2D_3-1], facecolor='none', edgecolor=color))
        elif nelnodes == 4:
            ax.add_patch(Polygon(x[f2D_4-1], facecolor='none', edgecolor=color))
        elif nelnodes == 6:
            ax.add_patch(Polygon(x[f2D_6-1], facecolor='none', edgecolor=color))
        elif nelnodes == 8 or nelnodes == 9:
            ax.add_patch(Polygon(x[f2D_8-1], facecolor='none', edgecolor=color))
    ax.autoscale(True)
    ax.set_aspect('equal')
    plt.show()


#######

# # PlotMesh untility has not added
# infile = "hyperelastic_quad4.txt"
# infile = open('hyperelastic_quad4.txt', 'r')

# Commented outfile for boosting 
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
        T = globaltraction(ncoord,ndof,nnode,ndload,coords, 
                    nelnodes,elident,connect,dloads,w)
        R = globalresidual(ncoord,ndof,nnode,coords, 
                nelem,maxnodes,elident,nelnodes, 
                connect,materialprops,w)

        b = loadfactor*T - R
         
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
    
    # Commented to boost the speed
    outfile.write(f'\n\n Step {step} Load {loadfactor}\n')

    # Commented to boost the speed
    print_results(outfile, 
        nprops,materialprops,ncoord,ndof,nnode,coords, 
        nelem,maxnodes,connect,nelnodes,elident, 
        nfix,fixnodes,ndload,dloads,w)
    
    # Store traction and displacement for plotting later
    forcevdisp[1,step] = loadfactor*dloads[2,0]
    forcevdisp[0,step] = w[2][0]
    
    # calculate deformed coordinates
    defcoords = np.zeros((ndof, nnode))
    scalefactor = 1.0
    for i in range(nnode):
        for j in range(ndof):
            defcoords[j, i] = coords[j, i] + scalefactor * w[ndof * (i - 1) + j]


#================================= POST-PROCESSING =================================



plotmesh(coords, ncoord, nnode, connect, nelem, elident, nelnodes, 'g')
# plotmesh(defcoords, ncoord, nnode, connect, nelem, elident, nelnodes, 'r')

plt.plot(forcevdisp[0,:], forcevdisp[1,:], 'r', linewidth=3)
plt.xlabel('Displacement', fontsize=16)
plt.ylabel('Force', fontsize=16)
plt.show()