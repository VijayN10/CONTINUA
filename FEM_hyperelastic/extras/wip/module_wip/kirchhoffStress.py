import numpy as np

def kirchhoffStress(ndof,ncoord,B,J,materialprops):
  stress = np.zeros((ndof,ncoord))
  dl = [[1,0,0],[0,1,0],[0,0,1]]

  mu1 = materialprops[0][0]
  K1 = materialprops[1][0]
  Bkk = sum([B[i][i] for i in range(ndof)])
  if (ndof==2):
      Bkk = Bkk + 1
  for i in range(ndof):
      for j in range(ncoord):
          stress[i][j] = mu1*(B[i][j] - Bkk*dl[i][j]/3.)/J**(2/3) + K1*J*(J-1)*dl[i][j]

  return stress