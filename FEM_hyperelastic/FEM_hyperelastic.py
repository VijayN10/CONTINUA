import numpy as np

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
elident = np.array([1])
connect = np.array([[1],[2],[3],[4]])

# No. nodes with prescribed displacements, with the prescribed displacements

nfix = 4
fixnodes = np.array([[1, 1, 2, 4],[1, 2, 2, 1],[0, 0, 0, 0]])

# No. loaded element faces, with the loads

ndload = 1
dloads = np.array([[1],[2],[3],[0]])


##############

def materialstiffness(ndof,ncoord,B,J,materialprops):
    
    mul = materialprops[0]
    K1 = materialprops[1]

    dl = np.array([[1,0,0],[0,1,0],[0,0,1]])

    C = np.zeros((ndof,ncoord,ndof,ncoord))
    
    if ncoord == 2:
        Bqq = B[0,1]+B[0,2] + 1
        for i in range(0,2):
            for j in range(0,2):
                for k in range(0,2):
                    for l in range(0,2):
                        C[i,j,k,l] =  mu1*( dl[i,k]*B[j,l]+B[i,l]*dl[j,k]\
                                    - (2/3)*(B[i,j]*dl[k,l]+dl[i,j]*B[k,l])\
                                    + (2/3)*Bqq*dl[i,j]*dl[k,l]/3 )/J^(2/3)\
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
            stress[i,j] = mu1*(B[i,j] - Bkk*dl[i,j]/3)/J^(2/3) + K1*J*(J-1)*dl[i,j]

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
      N[0] = 0.5*(1+xi[0])
      N[1] = 0.5*(1-xi[0])
    elif (nelnodes == 3):
      N[0] = -0.5*xi[0]*(1-xi[0])
      N[1] =  0.5*xi[0]*(1+xi[0])
      N[2] = (1-xi[0])*(1+xi[0])

  # %
  # %  2D elements
  # %
  elif (ncoord == 2):
  # 
  #     Triangular element
  # 
    if ( nelnodes == 3 ):
      N[0] = xi[0]
      N[1] = xi[0]
      N[2] = 1-xi[0]-xi[0]              
    elif ( nelnodes == 6 ):
      xi3 = 1-xi[0]-xi[0]
      N[0] = (2*xi[0]-1)*xi[0]
      N[1] = (2*xi[0]-1)*xi[0]
      N[2] = (2*xi3-1)*xi3
      N[3] = 4*xi[0]*xi[0]
      N[4] = 4*xi[0]*xi3
      N[5] = 4*xi3*xi[0]
  # 
  #     Rectangular element
  #                   
    elif ( nelnodes == 4 ):
      N[0] = 0.25*(1-xi[0])*(1-xi[0])
      N[1] = 0.25*(1+xi[0])*(1-xi[0])
      N[2] = 0.25*(1+xi[0])*(1+xi[0])
      N[4] = 0.25*(1-xi[0])*(1+xi[0])
    elif (nelnodes == 8):
      N[0] = -0.25*(1-xi[0])*(1-xi[0])*(1+xi[0]+xi[0])
      N[1] = 0.25*(1+xi[0])*(1-xi[0])*(xi[0]-xi[0]-1)
      N[2] = 0.25*(1+xi[0])*(1+xi[0])*(xi[0]+xi[0]-1)
      N[3] = 0.25*(1-xi[0])*(1+xi[0])*(xi[0]-xi[0]-1)
      N[4] = 0.5*(1-xi[0]*xi[0])*(1-xi[0])
      N[5] = 0.5*(1+xi[0])*(1-xi[0]*xi[0])
      N[6] = 0.5*(1-xi[0]*xi[0])*(1+xi[0])
      N[7] = 0.5*(1-xi[0])*(1-xi[0]*xi[0])

  return N

############################

def shapefunctionderivs(nelnodes,ncoord,elident,xi):
    
  dNdxi = np.zeros(nelnodes,ncoord)

# 
#  1D elements
# 
  if (ncoord == 1): 
      if (nelnodes==2): 
          dNdxi[0,0] = 0.5
          dNdxi[1,0] = -0.5
      elif (nelnodes == 3):
          dNdxi[0,0] = -0.5+xi[0]
          dNdxi[1,0] =  0.5+xi[0]
          dNdxi[2,0] = -2.*xi[0]

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
      xi3 = 1-xi[0]-xi[1]
      dNdxi[0,0] = 4*xi[0]-1
      dNdxi[1,1] = 4*xi[1]-1
      dNdxi[2,0] = -(4*xi3-1)
      dNdxi[2,1] = -(4*xi3-1)
      dNdxi[3,0] = 4*xi[1]
      dNdxi[3,1] = 4*xi[0]
      dNdxi[4,0] = -4*xi[1]
      dNdxi[4,1] = -4*xi[0]
      dNdxi[5,0] = 4*xi3 - 4*xi[0]
      dNdxi[5,1] = 4*xi3 - 4*xi[1]
# 
#     Rectangular element
#                   
    elif ( nelnodes == 4 ):
      dNdxi[0,0] = -0.25*(1-xi[1])
      dNdxi[1,2] = -0.25*(1-xi[0])
      dNdxi[1,0] = 0.25*(1-xi[1])
      dNdxi[1,1] = -0.25*(1+xi[0])
      dNdxi[2,0] = 0.25*(1+xi[1])
      dNdxi[2,1] = 0.25*(1+xi[0])
      dNdxi[3,0] = -0.25*(1+xi[1])
      dNdxi[3,1] = 0.25*(1-xi[0])
    elif (nelnodes == 8):
      dNdxi[0,0] = 0.25*(1-xi[1])*(2.*xi[0]+xi[1])
      dNdxi[1,2] = 0.25*(1-xi[0])*(xi[0]+2.*xi[1])
      dNdxi[1,0] = 0.25*(1-xi[1])*(2.*xi[0]-xi[1])
      dNdxi[1,1] = 0.25*(1+xi[0])*(2.*xi[1]-xi[0])
      dNdxi[2,0] = 0.25*(1+xi[1])*(2.*xi[0]+xi[1])
      dNdxi[2,1] = 0.25*(1+xi[0])*(2.*xi[1]+xi[0])
      dNdxi[3,0] = 0.25*(1+xi[1])*(2.*xi[0]-xi[1])
      dNdxi[3,1] = 0.25*(1-xi[0])*(2.*xi[1]-xi[0])
      dNdxi[4,0] = -xi[0]*(1-xi[1])
      dNdxi[4,1] = -0.5*(1-xi[0]*xi[0])
      dNdxi[5,0] = 0.5*(1-xi[1]*xi[1])
      dNdxi[5,1] = -(1+xi[0])*xi[1]
      dNdxi[6,0] = -xi[0]*(1+xi[1])
      dNdxi[6,1] = 0.5*(1-xi[0]*xi[0])
      dNdxi[7,0] = -0.5*(1-xi[1]*xi[1])
      dNdxi[7,1] = -(1-xi[0])*xi[1]
    
  return dNdxi

####################

def elresid(ncoord,ndof,nelnodes,elident,coord,materialprops,displacement):

  npoints = numberofintegrationpoints((ncoord,nelnodes,elident))
  dxdxi = np.zeros((ncoord,ncoord))
  dxidx = np.zeros((ncoord,ncoord))
  dNdxs = np.zeros((nelnodes,ncoord))
  rel = np.zeros((ndof*nelnodes,1))

  xilist = integrationpoints(ncoord,nelnodes,npoints,elident)
  w = integrationweights(ncoord,nelnodes,npoints,elident)

  for intpt in range(0,npoints):

    for i in range(0,ncoord):
      xi[i] = xilist[i,intpt]
   
    N = shapefunctions(nelnodes,ncoord,elident,xi)
    dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi)

    for i in range(0,ncoord):
      for j in range(0,ncoord):
        dxdxi[i,j] = 0.0
        for a in range(0,nelnodes):
          dxdxi[i,j] = dxdxi[i,j] + coord[i,a]*dNdxi[a,j]


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
            F[i,i] = F[i,j] + (displacement[i,a]*dNdx[a,j])

    J = np.linalg.det(F)
    B = F*(F.transpose())


    Finv = np.linalg.inv(F)
    for a in range(0,nelnodes):
      for i in range(0,ncoord):
        dNdxs[a.i] = 0.0
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

def elstif(ncoord,ndof,nelnodes,elident,coords,materialprops,displacement,xi,F, dNdxs):
    
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

        for i in range(0,ncoord):
            for j in range(0,ncoord):
                F[i,j] = 0.0
                if i==j: 
                    F[i,i] = 1
                    for a in range(0,nelnodes):
                        F[i,j] = F[i,j] + (displacement[i,a]*dNdx[a,j])

        # Compute Bbar and J

        J = np.linalg.det(F)
        B = F*(F.transpose())       


        # Convert shape function derivatives to derivatives wrt spatial coords

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

    list = np.zeros((nfacenodes(ncoord,nelnodes,face),1))
    
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

  npoints = numberofintegrationpoints(ncoord-1,nfacenodes)

  xi =  np.zeros((ncoord-1,1))
  dxdxi = np.zeros((ncoord,ncoord-1))
  r = np.zeros((ndof*nfacenodes,1))
  
  xilist = integrationpoints(ncoord-1,nfacenodes,npoints)
  w = integrationweights(ncoord-1,nfacenodes,npoints)

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
                lmncoord[i,a] = coords[i,connect[a,lmn]]
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

def globaltraction(ncoord,ndof,nnodes,ndloads,coords,
                   nelnodes,elident,connect,dloads,dofs):

  r = np.zeros((ndof*nnodes,1))
  traction = np.zeros((ndof,1))

  for load in range(0,ndloads):
# %
# %     Extract the coords of the nodes on the appropriate element face
# %
      lmn = int(dloads[0,load])
      face = int(dloads[1,load])
      n = nelnodes
      ident = elident[int(lmn)]
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
#                 print(rw)      
              r[rw] = r[rw] + rel[(a)*ndof+i]
  return F
      

###############

def globalresidual(ncoord,ndof,nnode,coords,nelem,maxnodes,elident,nelnodes,connect,materialprops,dofs):
   
  resid = np.zeros((ndof*nnode,1))
  lmncoord = np.zeros((ncoord,maxnodes))
  lmndof = np.zeros((ndof,maxnodes))
  rel = np.zeros((ndof*maxnodes,ndof*maxnodes))

  for lmn in range(0,nelem):

    for a in range(0,nelnodes[lmn]):
        
        for i in range(0,ncoord):
          lmncoord[i,a] = coords[i,connect[a,lmn]]
        for i in range (0,ndof):
          lmncoord[i,a] = dofs[ndof*(connect[a,lmn]-1)+i]

    n = nelnodes[lmn]
    ident = elident[lmn]
    rel = elresid(ncoord,ndof,n,ident,lmncoord,materialprops,lmndof)
   
    for a in range(0,nelnodes[lmn]):
      for i in range(0,ndof):
          rw = ndof*(connect[a,lmn]-1)+i
          resid[rw] = resid[rw] + rel[ndof*(a-1)+i]

  return resid

########################


w = np.zeros((nnode*ndof,1))


K = globalstiffness(ncoord,ndof,nnode,coords,
        nelem,maxnodes,elident,nelnodes,connect,materialprops,w)

F = globaltraction(ncoord,ndof,nnode,ndload,coords,
                nelnodes,elident,connect,dloads,w)