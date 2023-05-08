import numpy as np

def Voigt2normal(Dnumerical):
    
    dsde = np.zeros((2,2,2,2))
    
    dsde[0,0,0,0] = Dnumerical[0,0]
    dsde[0,0,1,1] = Dnumerical[0,1]
    dsde[0,0,0,1] = Dnumerical[0,2]
    dsde[0,0,1,0] = Dnumerical[0,2]

    dsde[1,1,0,0] = Dnumerical[1,0]
    dsde[1,1,1,1] = Dnumerical[1,1]
    dsde[1,1,0,1] = Dnumerical[1,2]
    dsde[1,1,1,0] = Dnumerical[1,2]

    dsde[0,1,0,0] = Dnumerical[2,0]
    dsde[0,1,1,1] = Dnumerical[2,1]
    dsde[0,1,0,1] = Dnumerical[2,2]
    dsde[0,1,1,0] = Dnumerical[2,2]

    dsde[1,0,0,0] = Dnumerical[2,0]
    dsde[1,0,1,1] = Dnumerical[2,1]
    dsde[1,0,0,1] = Dnumerical[2,2]
    dsde[1,0,1,0] = Dnumerical[2,2]

    return dsde