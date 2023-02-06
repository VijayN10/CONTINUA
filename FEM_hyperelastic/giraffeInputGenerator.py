import numpy as np
import os
import subprocess

####################################################################################################

nrve = 3    # number of rve

# Loop to create giraffe input files for each rve

for rr in range(nrve):

    # Path for Giraffe executables 
    path = "E:/Softwares/01/Giraffe/"
    giraffe = path + "Giraffe.exe"
    inp = 'tex_' + str(rr)    # Name of the Giraffe input file


    # Create new folder with name of .inp file
    folder = path + inp
    os.makedirs(folder, exist_ok=True)

    ######################################################################################################

    # Take inputs from text file and assign following variables 
    with open("data/rve"+ str(rr)  +".txt", "r") as file:
        Lxx, Lyy, t, nx, ny = [float(x) for x in file.read().split()]
        nx = int(nx)
        ny = int(ny)

    # Lxx = 6.4e-3
    # Lyy = 6.4e-3
    # t = 0.41e-3
    # nx = 2
    # ny = 2

    ######################################## INPUT OF F AND X ############################################
    # The F is used to calculate strain which then calulates alpha values 

    # Get the nodal coordinates of the ends of the beams
    X = np.zeros(((nx + ny) * 2, 3))

    with open("data/nodedata" + str(rr)  + ".txt", "r") as file:
        for i, line in enumerate(file):
            X[i] = [float(x) for x in line.strip().split()]
    print(f'X = {X}')                   

    # Deformation gradient
    F = np.array([[1.2, 0.2],
                [0.2, 1.2]])            # Assumption

    # Get the strain from the deformation gradient
    strain = F - np.eye(2)

    # Get alpha values from strain tensor
    alpha1 = strain[0,0]
    alpha2 = strain[0,1]
    alpha3 = strain[1,0]
    alpha4 = strain[1,1]

    #######################################################################################################

    # Convert alpha into a matrix
    alpha = np.array([[alpha1, alpha2],[alpha3, alpha4]])
    print(f'alpha = {alpha}')


    # Initialize disp and calculate the displacement
    # of the ends of the beams
    disp = np.zeros(((nx + ny) * 2,2))
    disp[:,0] = alpha1*X[:,0] + alpha2*X[:,1]
    disp[:,1] = alpha3*X[:,0] + alpha4*X[:,1]
    print(f'disp = {disp}')

    ############################################## Generate the Giraffe input file #######################################

    # Read the top portion from grfTop{i}.txt
    with open("data/grfTop" + str(rr)  + ".txt", "r") as top_file:
        top = top_file.read()

    # Read the bottom portion from grfBottom{i}.txt
    with open("data/grfBottom" + str(rr) + ".txt", "r") as bottom_file:
        bottom = bottom_file.read()


    # Open the file for writing

    with open(folder + '/' + inp + '.inp', 'w') as file:
        
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


    ############################################################### Run Giraffe.exe ################################################################

    # p = subprocess.Popen([giraffe], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # out, err = p.communicate(input=inp.encode() + b"\n")          # To write file name after executing giraffe
    # print(out.decode())

    # try:
    #     subprocess.run([giraffe], timeout=60)
    # except subprocess.TimeoutExpired:
    #     print("The .exe file took too long to run.")
