def create_displacement_block(nx, ny, disp):
    displacement_block = []

    for i in range(1, nx * ny + 1):
        block = "NodalDisplacement {} NodeSet {} CS 505 AngularParameters EulerVector NTimes 2\n//Time UX UY UZ ROTX ROTY ROTZ\n0\t0\t0\t0\t0\t0\t0\n1\t{}\t{}\t0\t0\t0\t0".format(i, i + 4, disp[i-1][0], disp[i-1][1])
        displacement_block.append(block)

    return displacement_block

disp = [(-0.55714285714285600, -0.35714285714285700),
        (-0.41428571428571300, -0.27142857142857100),
        (-0.27142857142857100, -0.18571428571428500),
        (-0.12857142857142800, -0.09999999999999940),
        (0.12857142857142800, 0.09999999999999950),
        (0.27142857142857100, 0.18571428571428500),
        (0.41428571428571300, 0.27142857142857100),
        (0.55714285714285600, 0.35714285714285700),
        (-0.54285714285714200, -0.34285714285714200),
        (-0.37142857142857100, -0.22857142857142800),
        (-0.19999999999999900, -0.11428571428571400),
        (-0.02857142857142790, 0.00000000000000047),
        (0.02857142857142790, -0.00000000000000044),
        (0.19999999999999900, 0.11428571428571400),
        (0.37142857142857100, 0.22857142857142800),
        (0.54285714285714200, 0.34285714285714200)]

nx = 4
ny = 4

displacement_block = create_displacement_block(nx, ny, disp)

print("Displacements {}".format(nx * ny))
for i in range(nx * ny):
    print(displacement_block[i])