import numpy as np
import sys
import matplotlib.pyplot as plt


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
width = 2
height = 2
maxnodes = 4
nelnodes = 4
coords, connect, elident, nnode, nelem, maxnodes, nelnodes = meshGen(nelemx, nelemy, width, height, maxnodes, nelnodes)

# # Print mesh information
# print(f"Mesh generated with {nelemx*nelemy} elements and {(nelemx+1)*(nelemy+1)} nodes.")
# print(f"Element connectivity:\n{connect}")
# print(f"Node coordinates:\n{coords}")
# print(f"elident:\n{elident}")
# print(f"nnode:\n{nnode}")
# print(f"nelem:\n{nelem}")
# print(f"maxnodes:\n{maxnodes}")
# print(f"nelnodes:\n{nelnodes}")

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



