import math
import matplotlib.pyplot as plt

from input import *

from rveDataLoader import loadRveData

from globalStiffness import globalStiffness
from globalResidual import globalResidual
from globalTraction import globalTraction

from printResults import printResults
from plotMesh import plotMesh



if model == 2:
    top, bottom, Lxx, Lyy, t, nx, ny, X = loadRveData(rvetype)
            




outfile = open('FEM_results.txt', 'w')

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
   
        K = globalStiffness(ncoord,ndof,nnode,coords, 
                nelem,maxnodes,elident,nelnodes,connect,materialprops,w)
        T = globalTraction(ncoord,ndof,nnode,ndload,coords, 
                    nelnodes,elident,connect,dloads,w)
        R = globalResidual(ncoord,ndof,nnode,coords, 
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
    printResults(outfile, 
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



plotMesh(coords, ncoord, nnode, connect, nelem, elident, nelnodes, 'g')
# plotMesh(defcoords, ncoord, nnode, connect, nelem, elident, nelnodes, 'r')

plt.plot(forcevdisp[0,:], forcevdisp[1,:], 'r', linewidth=3)
plt.xlabel('Displacement', fontsize=16)
plt.ylabel('Force', fontsize=16)
plt.show()
