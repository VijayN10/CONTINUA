import numpy as np

def loadRveData(rvetype):
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

    return top, bottom, Lxx, Lyy, t, nx, ny, X
