import numpy as np
import matplotlib.pyplot as plt
import math

# INPUT

# Total no. material parameters and list of parameters

nprops = 2
materialprops = np.array([[1], [10]])     # materialprops[mu1,K1]      # Python first position is 0 --> Need to change all code?
mu1 = 1
k1 = 10

# no. coords (1:3), no. DOF, no. nodes and nodal coordinates

ncoord = 2
ndof = 2
nnode = 4

coords = np.array([[0, 1, 1, 0],[0, 0, 1, 1]])

# No. elements and connectivity

nelem = 1
maxnodes = 4
nelnodes = 4
elident = np.array([[0]])                 # elident = np.array([[1]])
connect = np.array([[0],[1],[2],[3]])     # connect = np.array([[1],[2],[3],[4]])

# No. nodes with prescribed displacements, with the prescribed displacements

nfix = 4
fixnodes = np.array([[1, 1, 2, 4],[1, 2, 2, 1],[0, 0, 0, 0]])

# No. loaded element faces, with the loads

ndload = 1
dloads = np.array([[0],[1],[2],[0]])      # dloads = np.array([[1],[2],[3],[0]])
 

##############

def materialstiffness(ndof,ncoord,B,J,materialprops):
    
    mul = materialprops[0]
    K1 = materialprops[1]

    dl = np.array([[1,0,0],[0,1,0],[0,0,1]])

    C = np.zeros((ndof,ncoord,ndof,ncoord))
    
    if ncoord == 2:
        Bqq = B[0,0]+B[1,1] + 1
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    for l in range(0,2):
                        C[i,j,k,l] =  mu1*( dl[i,k]*B[j,l]+B[i,l]*dl[j,k]\
                                    - (2/3)*(B[i,j]*dl[k,l]+dl[i,j]*B[k,l])\
                                    + (2/3)*Bqq*dl[i,j]*dl[k,l]/3 )/J**(2/3)\
                                    + K1*(2*J-1)*J*dl[i,j]*dl[k,l]
                            
    return C 

################

def Kirchhoffstress(ndof,ncoord,B,J,materialprops):

    stress = np.zeros((ndof,ncoord))
    dl = np.array([[1,0,0],[0,1,0],[0,0,1]])

    # mul = materialprops[0]
    K1 = materialprops[1]

    Bkk = B.trace()
    if ndof==2: 
        Bkk = Bkk + 1
    for i in range(0,ndof):
        for j in range(0,ncoord):
            stress[i,j] = mu1*(B[i,j] - Bkk*dl[i,j]/3)/J**(2/3) + K1*J*(J-1)*dl[i,j]

    return stress

#################

def numberofintegrationpoints(ncoord,nelnodes,elident):
    if ncoord == 1:
        n = nelnodes  
    elif ncoord == 2:
        if nelnodes == 3:
            n = 1
        if nelnodes == 6:
            n = 3
        if nelnodes == 4:
            n = 4
        if nelnodes == 8:
            n = 9
    return n

##################

def integrationpoints(ncoord,nelnodes,npoints,elident):
    
  xi = np.zeros((ncoord,npoints))
# 
#   1D elements
# 
  if (ncoord == 1):
    if (npoints==1):
        xi[0,0] = 0
    elif (npoints == 2):
        xi[0,0] = -0.5773502692
        xi[0,1] = - xi[0,0]
    elif (npoints == 3):
        xi[0,0] = -0.7745966692
        xi[0,1] = 0.0
        xi[0,2] = -xi[0,0]

# 
#   2D elements
# 
  elif (ncoord == 2):
# 
#     Triangular element
# 
    if ( nelnodes == 3 or nelnodes == 6 ):
      if (npoints == 1):
        xi[0,0] = 1/3
        xi[1,0] = 1/3
      elif (npoints == 3):
        xi[0,0] = 0.6
        xi[1,0] = 0.2
        xi[0,1] = 0.2
        xi[1,1] = 0.6
        xi[0,2] = 0.2
        xi[1,2] = 0.2
      elif (npoints == 4):
        xi[0,0] = 1/3
        xi[1,0] = 1/3
        xi[0,1] = 0.6
        xi[1,1] = 0.2
        xi[0,2] = 0.2
        xi[1,2] = 0.6
        xi[0,3] = 0.2
        xi[1,3] = 0.2
# 
#     Rectangular element
#                   
    elif ( nelnodes==4 or nelnodes==8 ):

      if (npoints == 1):
        xi[0,0] = 0
        xi[1,0] = 0
      elif (npoints == 4):
        xi[0,0] = -0.5773502692
        xi[1,0] = xi[0,0]
        xi[0,1] = -xi[0,0]
        xi[1,1] = xi[0,0]
        xi[0,2] = xi[0,0]
        xi[1,2] = -xi[0,0]
        xi[0,3] = -xi[0,0]
        xi[1,3] = -xi[0,0]
      elif (npoints == 9): 
        xi[0,0] = -0.7745966692
        xi[1,0] = xi[0,0]
        xi[0,1] = 0.0
        xi[1,1] = xi[0,0]
        xi[0,2] = -xi[0,0]
        xi[1,2] = xi[0,0]
        xi[0,3] = xi[0,0]
        xi[1,3] = 0.0
        xi[0,4] = 0.0
        xi[1,4] = 0.0
        xi[0,5] = -xi[0,0]
        xi[1,5] = 0.0
        xi[0,6] = xi[0,0]
        xi[1,6] = -xi[0,0]
        xi[0,7] = 0
        xi[1,7] = -xi[0,0]
        xi[0,8] = -xi[0,0]
        xi[1,8] = -xi[0,0]

  return xi

###############

def integrationweights(ncoord,nelnodes,npoints,elident):
  w = np.zeros((npoints,1))

# 
#   1D elements
# 
  if (ncoord == 1):
    if (npoints == 1):
      w[0] = 2
    elif (npoints == 2):
      w = np.array([1,1])
    elif (npoints == 3):
      w = np.array([0.555555555,0.888888888,0.555555555])

# 
#   2D elements
# 
  elif (ncoord == 2):
# 
#     Triangular element
# 
    if ( nelnodes == 3 or nelnodes == 6 ):
      if (npoints == 1):
        w[0] = 0.5
      elif (npoints == 3):
        w[0] = 1/6
        w[1] = 1/6
        w[2] = 1/6
      elif (npoints == 4):
        w = np.array([-27/96,25/96,25/96,25/96])

# 
#     Rectangular element
#                   
    elif ( nelnodes==4 or nelnodes==8 ):

      if (npoints == 1):
        w[0] = 4
      elif (npoints == 4):
        w = np.array([1,1,1,1])
      elif (npoints == 9 ):
        w1D = np.array([0.555555555,0.888888888,0.55555555555])
        for j in range(0,3):
          for i in range(0,3):
            n = 3*(j-1)+i
            w[n] = w1D[i]*w1D[j]
  
  return w

##################

def shapefunctions(nelnodes,ncoord,elident,xi):
       
  N = np.zeros((nelnodes,1))

  # 
  #   1D elements
  # 
  if (ncoord == 1):
    if (nelnodes==2):
      N[0] = 0.5*(1+xi[0][0])
      N[1] = 0.5*(1-xi[0][0])
    elif (nelnodes == 3):
      N[0] = -0.5*xi[0][0]*(1-xi[0][0])
      N[1] =  0.5*xi[0][0]*(1+xi[0][0])
      N[2] = (1-xi[0][0])*(1+xi[0][0])

  # %
  # %  2D elements
  # %
  elif (ncoord == 2):
  # 
  #     Triangular element
  # 
    if ( nelnodes == 3 ):
      N[0] = xi[0][0]
      N[1] = xi[0][0]
      N[2] = 1-xi[0][0]-xi[0][0]              
    elif ( nelnodes == 6 ):
      xi3 = 1-xi[0][0]-xi[0][0]
      N[0] = (2*xi[0][0]-1)*xi[0][0]
      N[1] = (2*xi[0][0]-1)*xi[0][0]
      N[2] = (2*xi3-1)*xi3
      N[3] = 4*xi[0][0]*xi[0][0]
      N[4] = 4*xi[0][0]*xi3
      N[5] = 4*xi3*xi[0][0]
  # 
  #     Rectangular element
  #                   
    elif ( nelnodes == 4 ):
      N[0] = 0.25*(1-xi[0][0])*(1-xi[0][0])
      N[1] = 0.25*(1+xi[0][0])*(1-xi[0][0])
      N[2] = 0.25*(1+xi[0][0])*(1+xi[0][0])
      N[3] = 0.25*(1-xi[0][0])*(1+xi[0][0])
    elif (nelnodes == 8):
      N[0] = -0.25*(1-xi[0][0])*(1-xi[0][0])*(1+xi[0][0]+xi[0][0])
      N[1] = 0.25*(1+xi[0][0])*(1-xi[0][0])*(xi[0][0]-xi[0][0]-1)
      N[2] = 0.25*(1+xi[0][0])*(1+xi[0][0])*(xi[0][0]+xi[0][0]-1)
      N[3] = 0.25*(1-xi[0][0])*(1+xi[0][0])*(xi[0][0]-xi[0][0]-1)
      N[4] = 0.5*(1-xi[0][0]*xi[0][0])*(1-xi[0][0])
      N[5] = 0.5*(1+xi[0][0])*(1-xi[0][0]*xi[0][0])
      N[6] = 0.5*(1-xi[0][0]*xi[0][0])*(1+xi[0][0])
      N[7] = 0.5*(1-xi[0][0])*(1-xi[0][0]*xi[0][0])

  return N

############################

def shapefunctionderivs(nelnodes,ncoord,elident,xi):
    
  dNdxi = np.zeros((nelnodes,ncoord))

# 
#  1D elements
# 
  if (ncoord == 1): 
      if (nelnodes==2): 
          dNdxi[0,0] = 0.5
          dNdxi[1,0] = -0.5
      elif (nelnodes == 3):
          dNdxi[0,0] = -0.5+xi[0][0]
          dNdxi[1,0] =  0.5+xi[0][0]
          dNdxi[2,0] = -2.*xi[0][0]

# 
#   2D elements
#

  elif (ncoord == 2):
# 
#     Triangular element
# 
    if ( nelnodes == 3 ):
      dNdxi[0,0] = 1
      dNdxi[1,1] = 1
      dNdxi[2,0] = -1
      dNdxi[2,1] = -1               
    elif ( nelnodes == 6 ):
      xi3 = 1-xi[0][0]-xi[1][0]
      dNdxi[0,0] = 4*xi[0][0]-1
      dNdxi[1,1] = 4*xi[1][0]-1
      dNdxi[2,0] = -(4*xi3-1)
      dNdxi[2,1] = -(4*xi3-1)
      dNdxi[3,0] = 4*xi[1][0]
      dNdxi[3,1] = 4*xi[0][0]
      dNdxi[4,0] = -4*xi[1][0]
      dNdxi[4,1] = -4*xi[0][0]
      dNdxi[5,0] = 4*xi3 - 4*xi[0][0]
      dNdxi[5,1] = 4*xi3 - 4*xi[1][0]
# 
#     Rectangular element
#                   
    elif ( nelnodes == 4 ):
      dNdxi[0,0] = -0.25*(1-xi[1][0])
      dNdxi[0,1] = -0.25*(1-xi[0][0])
      dNdxi[1,0] = 0.25*(1-xi[1][0])
      dNdxi[1,1] = -0.25*(1+xi[0][0])
      dNdxi[2,0] = 0.25*(1+xi[1][0])
      dNdxi[2,1] = 0.25*(1+xi[0][0])
      dNdxi[3,0] = -0.25*(1+xi[1][0])
      dNdxi[3,1] = 0.25*(1-xi[0][0])
    elif (nelnodes == 8):
      dNdxi[0,0] = 0.25*(1-xi[1][0])*(2.*xi[0][0]+xi[1][0])
      dNdxi[1,2] = 0.25*(1-xi[0][0])*(xi[0][0]+2.*xi[1][0])
      dNdxi[1,0] = 0.25*(1-xi[1][0])*(2.*xi[0][0]-xi[1][0])
      dNdxi[1,1] = 0.25*(1+xi[0][0])*(2.*xi[1][0]-xi[0][0])
      dNdxi[2,0] = 0.25*(1+xi[1][0])*(2.*xi[0][0]+xi[1][0])
      dNdxi[2,1] = 0.25*(1+xi[0][0])*(2.*xi[1][0]+xi[0][0])
      dNdxi[3,0] = 0.25*(1+xi[1][0])*(2.*xi[0][0]-xi[1][0])
      dNdxi[3,1] = 0.25*(1-xi[0][0])*(2.*xi[1][0]-xi[0][0])
      dNdxi[4,0] = -xi[0][0]*(1-xi[1][0])
      dNdxi[4,1] = -0.5*(1-xi[0][0]*xi[0][0])
      dNdxi[5,0] = 0.5*(1-xi[1][0]*xi[1][0])
      dNdxi[5,1] = -(1+xi[0][0])*xi[1][0]
      dNdxi[6,0] = -xi[0][0]*(1+xi[1][0])
      dNdxi[6,1] = 0.5*(1-xi[0][0]*xi[0][0])
      dNdxi[7,0] = -0.5*(1-xi[1][0]*xi[1][0])
      dNdxi[7,1] = -(1-xi[0][0])*xi[1][0]
    
  return dNdxi

####################

def elresid(ncoord,ndof,nelnodes,elident,coords,materialprops,displacement):

  npoints = numberofintegrationpoints(ncoord,nelnodes,elident)
  dxdxi = np.zeros((ncoord,ncoord))
  dxidx = np.zeros((ncoord,ncoord))
  dNdxs = np.zeros((nelnodes,ncoord))
  rel = np.zeros((ndof*nelnodes,1))

  xilist = integrationpoints(ncoord,nelnodes,npoints,elident)
  w = integrationweights(ncoord,nelnodes,npoints,elident)
  xi = np.zeros((ncoord,npoints))  # Added TBC
  dNdx = np.zeros((nelnodes,ncoord)) # Added TBC
  coords = np.array([[0, 1, 1, 0],[0, 0, 1, 1]])  # Added here, if not added then it takes default as zero matrix creating problem with calculating dxdxi (giving singular matrix)
  

  F = np.zeros((ncoord,ncoord))  
  displacement = np.zeros((ncoord,nelnodes))   # added (initialize 2 x 4 matrix)

  for intpt in range(0,npoints):

    for i in range(0,ncoord):
      xi[i] = xilist[i,intpt]
   
    N = shapefunctions(nelnodes,ncoord,elident,xi)
    dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi)

    for i in range(0,ncoord):
      for j in range(0,ncoord):
        dxdxi[i,j] = 0.0
        for a in range(0,nelnodes):
          dxdxi[i,j] = dxdxi[i,j] + coords[i,a]*dNdxi[a,j]

    dxidx = np.linalg.inv(dxdxi)
    dt = np.linalg.det(dxdxi)

    for a in range(0,nelnodes):
      for i in range(0,ncoord):
        dNdx[a,i] = 0.0
        for j in range(0,ncoord):
          dNdx[a,i] = dNdx[a,i] + dNdxi[a,j]*dxidx[j,i]


    for i in range(0,ncoord):
      for j in range(0,ncoord):
        F[i,j] = 0.0
        if (i==j):
          F[i,i] = 1.0
        for a in range(0,nelnodes):
          F[i,j] = F[i,j] + (displacement[i,a]*dNdx[a,j])

    J = np.linalg.det(F)
    B = F*(F.transpose())


    Finv = np.linalg.inv(F)
    for a in range(0,nelnodes):
      for i in range(0,ncoord):
        dNdxs[a,i] = 0.0
        for j in range(0,ncoord):
          dNdxs[a,i] = dNdxs[a,i] + dNdx[a,j]*Finv[j,i]  


    stress = Kirchhoffstress(ndof,ncoord,B,J,materialprops)
            
    for a in range(0,nelnodes):
      for i in range(0,ndof):
        row = ndof*(a-1)+i
        for j in range(0,ncoord):
          rel[row] = rel[row] + stress[i,j]*dNdxs[a,j]*w[intpt]*dt

  return rel

#######################

def elstif(ncoord,ndof,nelnodes,elident,coords,materialprops,displacement):
    
    npoints = numberofintegrationpoints(ncoord,nelnodes,elident)
    dNdx = np.zeros((nelnodes,ncoord))
    dxdxi = np.zeros((ncoord,ncoord))
    strain = np.zeros((ndof,ncoord))
    kel = np.zeros((ndof*nelnodes,ndof*nelnodes))
    
    # Set up integration points && weights

    xilist = integrationpoints(ncoord,nelnodes,npoints,elident)
    w = integrationweights(ncoord,nelnodes,npoints,elident)

    # Loop over the integration points

    

    for intpt in range(0,npoints):

        # Compute shape functions && derivatives wrt local coords

        xi = np.zeros((ncoord,1))  #added  ncoord - 1 ---> ncoords (1,1) ---> (2,1)
        
        for i in range(0,ncoord):
          xi[i] = xilist[i,intpt]
        
        N = shapefunctions(nelnodes,ncoord,elident,xi)
        dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi)

        # Compute the jacobian matrix && its determinant
        
        for i in range(0,ncoord):
            for j in range(0,ncoord):
                dxdxi[i,j] = 0.0
                for a in range(0,nelnodes):
                    dxdxi[i,j]= dxdxi[i,j] + coords[i,a]*dNdxi[a,j]
        
        dxidx = np.linalg.inv(dxdxi)
        dt = np.linalg.det(dxdxi)
        
        # Convert shape function derivatives:derivatives wrt global coords

        for a in range(0,nelnodes):
            for i in range(0,ncoord):
                dNdx[a,i] = 0.0
                for j in range(0,ncoord):
                    dNdx[a,i] = dNdx[a,i] + dNdxi[a,j]*dxidx[j,i]
        

        # Compute the deformation gradients by differentiating displacements

        F = np.zeros((ncoord,ncoord))  # TBC Calculated in elresid function?

        displacement = np.zeros((ncoord,nelnodes))   # added (initialize 2 x 4 matrix)

        # print(F)
        # print(F.shape)

        for i in range(0,ncoord):
            for j in range(0,ncoord):
                F[i,j] = 0.0
                if i==j: 
                    F[i,i] = 1
                    for a in range(0,nelnodes):
                        F[i,j] = F[i,j] + (displacement[i,a]*dNdx[a,j])

        # print(F)
        # print(F.shape)
        
        # Compute Bbar and J

        J = np.linalg.det(F)
        B = F*(F.transpose())       


        # Convert shape function derivatives to derivatives wrt spatial coords

        dNdxs = np.zeros((nelnodes,ncoord))   # Is it correct?

        Finv = np.linalg.inv(F)
        for a in range(0,nelnodes):
          for i in range(0,ncoord):
            dNdxs[a,i] = 0.0
            for j in range(0,ncoord):
              dNdxs[a,i] = dNdxs[a,i] + dNdx[a,j]*Finv[j,i]
      
        # Compute the stress

        stress = Kirchhoffstress(ndof,ncoord,B,J,materialprops)

        # Compute the material tangent stiffness (d stress/d strain) ds/de is just C_ijkl for linear elasticity - 
        # this notation is used to allow extension to nonlinear problems
        
        dsde = materialstiffness(ndof,ncoord,B,J,materialprops)

        # Compute the element stiffness

        for a in range(0,nelnodes):
            for i in range(0,ndof):
                for b in range(0,nelnodes):
                    for k in range(0,ndof):
                        row = ndof*(a) + i
                        col = ndof*(b) + k
                        for l in range(0,ncoord):
                            for j in range(0, ncoord):
                                kel[row,col] = kel[row,col] + dsde[i,j,k,l]*dNdxs[b,l]*dNdxs[a,j]*w[intpt]*dt
                                kel[row,col] = kel[row,col] - stress[i,j]*dNdxs[a,k]*dNdxs[b,j]*w[intpt]*dt

    return kel

##############

def nfacenodes(ncoord,nelnodes,elident,face):
    if (ncoord == 2): 
        if (nelnodes == 3 or nelnodes == 4):
            n = 2
        elif (nelnodes == 6 or nelnodes == 8): 
            n = 3
    return n

############

def facenodes(ncoord,nelnodes,elident,face):
    
    i3 = np.array([2,3,1])
    i4 = np.array([2,3,4,1])

    list = np.zeros((nfacenodes(ncoord,nelnodes,elident,face),1))
    
    if ncoord == 2:
        if nelnodes ==3:
            list[0] = face
            list[1] = i3[face]
        elif (nelnodes == 6): 
            list[0] = face
            list[1] = i3[face]
            list[2] = face+3
        elif (nelnodes==4):
            list[0] = face
            list[1] = i4[face]
        elif (nelnodes==8):
            list[0] = face
            list[1] = i4[face]
            list[2] = face+4
            
    return list

###############################

def eldload(ncoord,ndof,nfacenodes,elident,coords,traction):

  npoints = numberofintegrationpoints(ncoord-1,nfacenodes,elident)

  xi =  np.zeros((ncoord-1,1))
  dxdxi = np.zeros((ncoord,ncoord-1))
  r = np.zeros((ndof*nfacenodes,1))
  
  xilist = integrationpoints(ncoord-1,nfacenodes,npoints,elident)
  w = integrationweights(ncoord-1,nfacenodes,npoints,elident)

  for intpt in range(0,npoints):

    for i in range(0,ncoord-1):
      xi[i] = xilist[i,intpt]

    N = shapefunctions(nfacenodes,ncoord-1,elident,xi)
    dNdxi = shapefunctionderivs(nfacenodes,ncoord-1,elident,xi)

# %
# %     Compute the jacobian matrix && its determinant
# %
    for i in range(0,ncoord):
        for j in range(0,ncoord-1):
            dxdxi[i,j] = 0.0
            for a in range(0,nfacenodes):
                dxdxi[i,j] = dxdxi[i,j] + coords[i,a]*dNdxi[a,j]

    if ncoord ==2:
        dt = np.sqrt(dxdxi[0,0]**2+dxdxi[1,0]**2)

    for a in range(0,nfacenodes):
        for i in range(0,ndof):
            row = ndof*(a-1)+i;
            r[row] = r[row] + N[a]*traction[i]*w[intpt]*dt
            
  return r

####################

def globalstiffness(ncoord,ndof,nnode,coords,nelem,
                    maxnodes,elident,nelnodes,connect,materialprops,dofs):

    Stif = np.zeros((ndof*nnode,ndof*nnode))
    lmncoord = np.zeros((ncoord,maxnodes))
    lmndof = np.zeros((ndof,maxnodes))

    # %
    # %   Loop over all the elements
    # %
    for lmn in range(0,nelem): 
    # %
    # %   Extract coords of nodes, DOF for the current element
    # %
        for a in range(0,nelnodes):
            for i in range(0,ncoord):
                lmncoord[i,a] = coords[i,connect[a,lmn]]  # Coords is having dimensions 2x4, thus j should be 0,1,2,3. But, connect[3,0]=4. coords[i,4] gives en error; Solution --> coords[i,4-1]  OR Changed connect to start from 0 instead of 1
            for i in range(0,ndof):
                lmndof[i,a] = dofs[ndof*(connect[a,lmn]-1)+i]

        n = nelnodes
        ident = elident[lmn] 
        kel = elstif(ncoord,ndof,n,ident,lmncoord,materialprops,lmndof)
    # %
    # %   Add the current element stiffness:the global stiffness
    # %
        for a in range(0,nelnodes):
            for i in range(0,ndof):
                for b in range(0,nelnodes):
                    for k in range(0,ndof):
                        rw = ndof*(connect[a,lmn])+i
                        cl = ndof*(connect[b,lmn])+k
                        Stif[rw,cl] = Stif[rw,cl] + kel[ndof*(a)+i,ndof*(b)+k]
    
    return Stif

################

def globaltraction(ncoord,ndof,nnodes,ndload,coords,
                   nelnodes,elident,connect,dloads,dofs):

  r = np.zeros((ndof*nnodes,1))
  traction = np.zeros((ndof,1))

  for load in range(0,ndload):
# %
# %     Extract the coords of the nodes on the appropriate element face
# %
      lmn = int(dloads[0,load])
      face = int(dloads[1,load])
      n = nelnodes
      ident = elident[lmn]
      nfnodes = nfacenodes(ncoord,n,ident,face)
      nodelist = facenodes(ncoord,n,ident,face)   
      lmncoord = np.zeros((ncoord,nfnodes))
      for a in range(0,nfnodes):
          for i in range(0,ncoord):
              lmncoord[i,a] = coords[i,connect[int(nodelist[a]),int(dloads[0,load])]]
        #   for i in range(0,ndof):
        #       lmndof[i,a] = dofs(ndof*(connect[nodelist[a],dloads[0,load]]-1)+i)

# %
# %    Compute the element load vector
# %
      for i in range(0,ndof):
          traction[i] = dloads[i+2,load]                    
      rel = eldload(ncoord,ndof,nfnodes,ident,lmncoord,traction)

# %
# %    Assemble the element load vector into global vector
# %
      for a in range(0,nfnodes):
          for i in range(0,ndof):
              rw = (connect[int(nodelist[a]),int(dloads[0,load])])*ndof + i
    
              r[rw] = r[rw] + rel[(a)*ndof+i]
  return r
      
###############

def globalresidual(ncoord,ndof,nnode,coords,nelem,maxnodes,elident,nelnodes,connect,materialprops,dofs):
   
  resid = np.zeros((ndof*nnode,1))
  lmncoord = np.zeros((ncoord,maxnodes))
  lmndof = np.zeros((ndof,maxnodes))
  rel = np.zeros((ndof*maxnodes,ndof*maxnodes))

  for lmn in range(0,nelem):

    for a in range(0,nelnodes):          # range(0,nelnodes[lmn]):
        
        for i in range(0,ncoord):
          lmncoord[i,a] = coords[i,connect[a,lmn]]
        for i in range (0,ndof):
          lmncoord[i,a] = dofs[ndof*(connect[a,lmn]-1)+i]

    n = nelnodes                         # n = nelnodes[lmn]
    ident = elident[lmn]
    rel = elresid(ncoord,ndof,n,ident,lmncoord,materialprops,lmndof)
   
    for a in range(0,nelnodes):
      for i in range(0,ndof):
          rw = ndof*(connect[a,lmn]-1)+i
          resid[rw] = resid[rw] + rel[ndof*(a-1)+i]

  return resid

########################

w = np.zeros((nnode*ndof,1))

nsteps = 5
tol = 0.0001
maxit = 1
relax = 1.0

forcevdisp = np.zeros((2,nsteps+1)) 
forcevdisp[0,0] = 0.0
forcevdisp[1,0] = 0.0


for step in range(0,nsteps):

  loadfactor = (step+1)/nsteps    # +1 added 

  err1 = 1.0
  nit = 0
  
  while (err1 > tol and nit < maxit):
    nit += 1

    K = globalstiffness(ncoord, ndof, nnode, coords, nelem, maxnodes, elident, nelnodes, connect, materialprops, w)
    F = globaltraction(ncoord, ndof, nnode, ndload, coords, nelnodes, elident, connect, dloads, w)
    R = globalresidual(ncoord, ndof, nnode, coords, nelem, maxnodes, elident, nelnodes, connect, materialprops, w)

    b = loadfactor * F - R

    # Fix constrained nodes
    for n in range(nfix):
        rw = ndof * (fixnodes[0][n] - 1) + fixnodes[1][n]
        for cl in range(ndof * nnode):
            K[rw][cl] = 0
        K[rw][rw] = 1.
        b[rw] = loadfactor * fixnodes[2][n] - w[rw]

    # Solve for the correction
    dw = np.linalg.solve(K, b)

    # Check convergence
    w += relax * dw
    wnorm = np.dot(w.flatten(), w.flatten())
    err1 = np.dot(dw.flatten(), dw.flatten())
    err2 = np.dot(b.flatten(), b.flatten())
    err1 = math.sqrt(err1 / wnorm)
    err2 = math.sqrt(err2) / (ndof * nnode)
    print("Iteration number", nit, "Correction", err1, "Residual", err2, "tolerance", tol)


  with open("FEM_results.txt", "w") as outfile:
    outfile.write('\n\n Step {} Load {}\n'.format(step, loadfactor))


  forcevdisp[0, step+1] = loadfactor*dloads[2, 0]
  forcevdisp[1, step+1] = w[2]


##########################

plt.plot(forcevdisp[0], forcevdisp[1], 'r', linewidth=3)
plt.xlabel('Displacement', fontsize=16)
plt.ylabel('Force', fontsize=16)
plt.show()

#########################


# def plotmesh(coords, ncoord, nnode, connect, nelem, elident, nelnodes, color):
#     f2D_3 = [1, 2, 3]
#     f2D_4 = [1, 2, 3, 4]
#     f2D_6 = [1, 4, 2, 5, 3, 6]
#     f2D_8 = [1, 5, 2, 6, 3, 7, 4, 8]

#     if ncoord == 2:
#         for lmn in range(nelem):
#             x = []
#             for i in range(nelnodes):
#                 x.append(coords[:, connect[i][lmn]])

#             plt.scatter(x[0], x[1], marker='o', color='r')

#             if nelnodes == 3:
#                 plt.triplot(x[0], x[1], f2D_3, color=color)
#             elif nelnodes == 4:
#                 plt.triplot(x[0], x[1], f2D_4, color=color)
#             elif nelnodes == 6:
#                 plt.triplot(x[0], x[1], f2D_6, color=color)
#             elif nelnodes == 8 or nelnodes == 9:
#                 plt.triplot(x[0], x[1], f2D_8, color=color)

#     plt.axis('equal')


###############################

defcoords = np.zeros((ndof,nnode))
scalefactor = 1.0
for i in range(0,nnode):
  for j in range(0,ndof):
      defcoords[j][i] = coords[j][i] + scalefactor*w[ndof*(i-1)+j-1]  # Check?

# plt.figure()
# plt.gca().set_aspect('equal')
# plotmesh(coords,ncoord,nnode,connect,nelem,elident,nelnodes,'g')
# plt.hold(True)
# plotmesh(defcoords,ncoord,nnode,connect,nelem,elident,nelnodes,'r')
# plt.show()
#########################

# outfile= open("FEM_results.txt","w+")

# def print_results(outfile, nprops,materialprops,ncoord,ndof,nnode,coords,nelem,maxnodes,connect,nelnodes,elident,nfix,fixnodes,ndload,dloads,dofs):



  # print(outfile,'Nodal Displacements: \n')
#    if (ndof == 2) 
#      fprintf(outfile,' Node      Coords         u1       u2 \n');
#      for i = 1:nnode
#       fprintf(outfile,'%3d %8.4f %8.4f %8.4f %8.4f\n', ...
#                                i,coords(1,i),coords(2,i),dofs(2*i-1),dofs(2*i));
#      end
#    elseif (ndof == 3) 
#      fprintf(outfile,' Node            Coords            u1       u2       u3 \n');
#      for i = 1:nnode
#       fprintf(outfile,'%3d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n', ...
#                     i,coords(1,i),coords(2,i),coords(3,i),dofs(3*i-2),dofs(3*i-1),dofs(3*i));
#      end
#    end

#    fprintf(outfile,'\n\n Strains and Stresses \n');


#    lmncoord = zeros(ncoord,maxnodes);
#    displacements = zeros(ndof,maxnodes);

# %
# %   Loop over all the elements
# %
#    for lmn = 1:nelem

#     fprintf(outfile,' \n Element; %d ',lmn);
#     if (ncoord == 2)   
#     fprintf(outfile,'  \n int pt    Coords          B_11      B_22     B_12      s_11       s_22      s_12 \n');

#     elseif (ncoord == 3) 
#     fprintf(outfile,'\n int pt         Coords            B_11      B_22     B_33      B_12       B_13      B_23      s_11      s_22      s_33      s_12      s_13      s_23 \n');
#     end
# %
# %   Extract coords of nodes, DOF for the current element
# %
#       for a = 1:nelnodes(lmn)
#         for i = 1:ncoord
#           lmncoord(i,a) = coords(i,connect(a,lmn));
#         end
#         for i = 1:ndof
#           displacements(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
#         end
#       end
#       n = nelnodes(lmn);
#       ident = elident(lmn);
 
#       npoints = numberofintegrationpoints(ncoord,n);
#       dNdx = zeros(n,ncoord);
#       dxdxi = zeros(ncoord,ncoord);
#       xi = zeros(ncoord,1);
#       x = zeros(ncoord,1);
# %
# %  Set up integration points 
# %
#       xilist = integrationpoints(ncoord,n,npoints);
# %
# %  Loop over the integration points
# %
#      for intpt = 1:npoints

# %     Compute shape functions && derivatives wrt local coords
# %
#        for i = 1:ncoord
#          xi(i) = xilist(i,intpt);
#        end
#        N = shapefunctions(n,ncoord,ident,xi);      
#        dNdxi = shapefunctionderivs(n,ncoord,ident,xi);
# %
# %     Compute the coords of the integration point
# %
#       for i = 1:ncoord
#         x(i) = 0.;
#         for a = 1:n
#           x(i) = x(i) + lmncoord(i,a)*N(a);
#         end
#       end
# %
# %     Compute the jacobian matrix && its determinant
# %
#       for i = 1:ncoord
#         for j = 1:ncoord
#           dxdxi(i,j) = 0.;
#           for a = 1:n
#             dxdxi(i,j) = dxdxi(i,j) + lmncoord(i,a)*dNdxi(a,j);
#           end
#         end
#       end

#       dxidx = inv(dxdxi);
# %
# %     Convert shape function derivatives:derivatives wrt global coords
# %
#       for a = 1:n
#         for i = 1:ncoord
#           dNdx(a,i) = 0.;
#           for j = 1:ncoord
#             dNdx(a,i) = dNdx(a,i) + dNdxi(a,j)*dxidx(j,i);
#           end
#         end
#       end
# %
# %     Compute the deformation gradients by differentiating displacements
# %
#       for i = 1 : ncoord
#          for j = 1 : ncoord
#             F(i,j) = 0.;
#             if (i==j) F(i,i) = 1.; end
#             for a = 1 : nelnodes
#               F(i,j) = F(i,j) ...
#                          +(displacements(i,a)*dNdx(a,j));
#             end
#          end
#       end
# %
# %     Compute Bbar and J
# %
#       J = det(F);
#       B = F*transpose(F);
# %
# %     Convert shape function derivatives to derivatives wrt spatial coords
# %
#       Finv = inv(F);
#       for a = 1 : nelnodes
#         for i = 1 : ncoord
#           dNdxs(a,i) = 0.;
#           for j = 1 : ncoord
#           dNdxs(a,i) = dNdxs(a,i) + dNdx(a,j)*Finv(j,i);  
#           end
#         end
#       end
# %
# %     Compute the stress
# %
#       stress = Kirchhoffstress(ndof,ncoord,B,J,materialprops);
#       stress = stress/J;   % Compute Cauchy stress
      
#       if (ncoord == 2) 

#       fprintf(outfile,'%5d %7.4f %7.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f \n', ...
#         intpt,x(1),x(2),B(1,1),B(2,2),B(1,2),stress(1,1),stress(2,2),stress(1,2));


#       elseif (ncoord == 3) 

#       fprintf(outfile,'%5d %7.4f %7.4f %7.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f \n',...
#               intpt,x(1),x(2),x(3), ...
#               B(1,1),B(2,2),B(3,3),B(1,2),B(1,3),B(2,3), ...
#               stress(1,1),stress(2,2),stress(3,3),stress(1,2),stress(1,3),stress(2,3));
