import numpy as np
import os

from input import rvetype
from rveDataLoader import loadRveData



top, bottom, Lxx, Lyy, t, nx, ny, X = loadRveData(rvetype)

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