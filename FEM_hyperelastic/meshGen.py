import numpy as np
import sys


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
    

    # Check if only one element is present
    if nelem == 1:
        proceed = input("Meshing not done: only one element present. Do you want to proceed? (y/n)")
        if proceed != "y":
            sys.exit()

    # Return mesh information
    return coords, connect, elident, nnode, nelem, maxnodes, nelnodes


nelemx = 3
nelemy = 3
width = 3
height = 3
maxnodes = 4
nelnodes = 4
coords, connect, elident, nnode, nelem, maxnodes, nelnodes = meshGen(nelemx, nelemy, width, height, maxnodes, nelnodes)

# Print mesh information
print(f"Mesh generated with {nelemx*nelemy} elements and {(nelemx+1)*(nelemy+1)} nodes.")
print(f"Element connectivity:\n{connect}")
print(f"Node coordinates:\n{coords}")
print(f"elident:\n{elident}")
print(f"nnode:\n{nnode}")
print(f"nelem:\n{nelem}")
print(f"maxnodes:\n{maxnodes}")
print(f"nelnodes:\n{nelnodes}")