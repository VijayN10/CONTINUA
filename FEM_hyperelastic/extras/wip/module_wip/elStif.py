
import numpy as np

from input import *

from numberOfIntegrationPoints import numberOfIntegrationPoints
from integrationPoints import integrationPoints
from integrationWeights import integrationWeights
from shapeFunctions import shapeFunctions
from shapeFunctionDerivs import shapeFunctionDerivs
from kirchhoffStress import kirchhoffStress
from materialStiffness import materialStiffness
from giraffeInputGenerator import giraffeInputGenerator
from runGiraffe import runGiraffe


from rveDataLoader import loadRveData



top, bottom, Lxx, Lyy, t, nx, ny, X = loadRveData(rvetype)




def elStif(ncoord, ndof, nelnodes, elident, coord, materialprops, displacement):
    
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

    npoints = numberOfIntegrationPoints(ncoord, nelnodes, elident)
    dNdx = np.zeros((nelnodes, ncoord))
    dxdxi = np.zeros((ncoord, ncoord))
    strain = np.zeros((ndof, ncoord))
    kel = np.zeros((ndof*nelnodes, ndof*nelnodes))

    # Set up integration points && weights 

    xilist = integrationPoints(ncoord,nelnodes,npoints,elident)
    w = integrationWeights(ncoord,nelnodes,npoints,elident)

    # Loop over the integration points

    for intpt in range(npoints):

        # Compute shape functions && derivatives wrt local coords

        xi = np.zeros((ncoord,1))
        for i in range(0,ncoord):
          xi[i] = xilist[i,intpt] 

        N = shapeFunctions(nelnodes, ncoord, elident, xi)
        dNdxi = shapeFunctionDerivs(nelnodes, ncoord, elident, xi)

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

            stress = kirchhoffStress(ndof, ncoord, B, J, materialprops)

            # Compute the material tangent stiffness (d stress/d strain)
            # ds/de is just C_ijkl for linear elasticity - this notation is used
            # to allow extension to nonlinear problems

            dsde = materialStiffness(ndof, ncoord, B, J, materialprops)
        

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