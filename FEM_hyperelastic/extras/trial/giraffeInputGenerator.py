import numpy as np
import os


rvetype = 0
name = 'tex'
F = np.array([[1.2, 0.2],   # Assumption used for testing this function
            [0.2, 1.2]])


def giraffeInputGenerator(rvetype, name, F):

    # Consider rvetype                 
    rr = rvetype  

    # Create new folder with name of .inp file
    inp = name + '_' + str(rr)                       # inp = tex_0
    folder = "Giraffe/" + inp                              # folder = E:/Softwares/01/Giraffe/tex_0
    os.makedirs("Giraffe/" + inp, exist_ok=True)               # Creates tex_0 folder

    # Take inputs from text file and assign following variables 
    with open("Giraffe/RVE_data/rve"+ str(rr)  +".txt", "r") as file:
        Lxx, Lyy, t, nx, ny = [float(x) for x in file.read().split()]
        nx = int(nx)
        ny = int(ny)

    # Get the nodal coordinates of the ends of the beams

    X = np.zeros(((nx + ny) * 2, 3))

    # Read the nodedata text files from data folder
    with open("Giraffe/RVE_data/nodedata" + str(rr)  + ".txt", "r") as file:
        for i, line in enumerate(file):
            X[i] = [float(x) for x in line.strip().split()]
    
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
    disp = np.zeros(((nx + ny) * 2,2))
    disp[:,0] = alpha1*X[:,0] + alpha2*X[:,1]
    disp[:,1] = alpha3*X[:,0] + alpha4*X[:,1]
    # print(f'disp = {disp}')

    # Read the top portion from grfTop{i}.txt
    with open("Giraffe/RVE_data/grfTop" + str(rr)  + ".txt", "r") as top_file:
        top = top_file.read()

    # Read the bottom portion from grfBottom{i}.txt
    with open("Giraffe/RVE_data/grfBottom" + str(rr) + ".txt", "r") as bottom_file:
        bottom = bottom_file.read()


    # Open the file for writing

    filepath = 'Giraffe/' + inp + '/' + inp + '.inp'         # filepath - E:/Softwares/01/Giraffe/tex_0/tex_0.inp
    with open(filepath, 'w') as file:
        
    # Top part of Giraffe input file

        file.write(top) 

    # Displacement block

    # Left

        file.write("""
        
    Displacements 8
    ///////////////////////////////
    ///////////left////////////////
    ///////////////////////////////
    NodalDisplacement 1 NodeSet 5 CS 125 NTimes 2
    //Time UX UY UZ ROTX ROTY ROTZ
    0	0	0	0	0	0	0     
        """)   
        file.write(f'2 {disp[0,0]:.12f} {disp[0,1]:.12f} 0 0 0 0\n')
        file.write("""
    NodalDisplacement 2 NodeSet 6 CS 125 NTimes 2
    //Time UX UY UZ ROTX ROTY ROTZ
    0	0	0	0	0	0	0
        """)
        file.write(f'2 {disp[1,0]:.12f} {disp[1,1]:.12f} 0 0 0 0\n\n')

    # Right

        file.write("""
    ///////////////////////////////
    ///////////right///////////////
    ///////////////////////////////
    NodalDisplacement 3 NodeSet 7 CS 125 NTimes 2
    //Time UX UY UZ ROTX ROTY ROTZ
    0	0	0	0	0	0	0
        """)
        file.write(f'2 {disp[2,0]:.12f} {disp[2,1]:.12f} 0 0 0 0\n')
        file.write("""
    NodalDisplacement 4 NodeSet 8 CS 125 NTimes 2
    //Time UX UY UZ ROTX ROTY ROTZ
    0	0	0	0	0	0	0
        """)
        file.write(f'2 {disp[3,0]:.12f} {disp[3,1]:.12f} 0 0 0 0\n\n')

    # Bottom

        file.write("""
    ///////////////////////////////
    ///////////bottom//////////////
    ///////////////////////////////
    NodalDisplacement 5 NodeSet 9 CS 125 NTimes 2
    //Time UX UY UZ ROTX ROTY ROTZ
    0	0	0	0	0	0	0
        """)
        file.write(f'2 {disp[4,0]:.12f} {disp[4,1]:.12f} 0 0 0 0\n')
        file.write("""
    NodalDisplacement 6 NodeSet 10 CS 125 NTimes 2
    //Time UX UY UZ ROTX ROTY ROTZ
    0	0	0	0	0	0	0
        """)
        file.write(f'2 {disp[5,0]:.12f} {disp[5,1]:.12f} 0 0 0 0\n\n')

    # Top

        file.write("""
    ///////////////////////////////
    ///////////top/////////////////
    ///////////////////////////////
    NodalDisplacement 7 NodeSet 11 CS 125 NTimes 2
    //Time UX UY UZ ROTX ROTY ROTZ
    0	0	0	0	0	0	0
        """)
        file.write(f'2 {disp[6,0]:.12f} {disp[6,1]:.12f} 0 0 0 0\n')
        file.write("""
    NodalDisplacement 8 NodeSet 12 CS 125 NTimes 2
    //Time UX UY UZ ROTX ROTY ROTZ
    0	0	0	0	0	0	0
        """)
        file.write(f'2 {disp[7,0]:.12f} {disp[7,1]:.12f} 0 0 0 0\n\n')


    # Bottom part of Giraffe input file

        file.write(bottom)


    return inp, Lxx, Lyy, t, folder


inp, Lxx, Lyy, t, folder = giraffeInputGenerator(rvetype, name, F)