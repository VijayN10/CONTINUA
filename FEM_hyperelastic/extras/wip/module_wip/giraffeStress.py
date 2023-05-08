import numpy as np

#@jit(nopython=False, parallel=True)
def giraffeStress(Lxx, Lyy, t, inp):
      
    # # Changing t temporarily    # Ask
    # t = 1
    
    # File for left side
    # monfilepath = "4x4Job/monitors/monitor_nodeset_2.txt"
    monfilepath = 'Giraffe/' + inp + '/monitors/monitor_nodeset_2.txt'
    data_left = np.genfromtxt(monfilepath, skip_header=1)
    posx = np.abs(data_left[:,1])
    posx = posx - posx[0]
    posy = np.abs(data_left[:,2])
    posy = posy - posy[0]
    forcex = np.abs(data_left[:,7])
    forcey = np.abs(data_left[:,8])
    data_left = np.column_stack((posx, posy, forcex, forcey))

    # File for top side
    # monfilepath = "4x4Job/monitors/monitor_nodeset_4.txt"
    monfilepath = 'Giraffe/' + inp + '/monitors/monitor_nodeset_4.txt'
    data_top = np.genfromtxt(monfilepath, skip_header=1)
    posx = np.abs(data_top[:,1])
    posx = posx - posx[0]
    posy = np.abs(data_top[:,2])
    posy = posy - posy[0]
    forcex = np.abs(data_top[:,7])
    forcey = np.abs(data_top[:,8])
    data_top = np.column_stack((posx, posy, forcex, forcey))

    # Get the forces
    Fxx = data_left[-1,2]
    Fyy = data_top[-1,3]
    Fxy = data_left[-1,3]
    Fyx = data_top[-1,2]

    # Stress
    stress = np.zeros((2,2))
    stress[0,0] = Fxx/(1000*Lyy*t)
    stress[1,1] = Fyy/(1000*Lxx*t)
    stress[0,1] = Fxy/(1000*Lyy*t)
    stress[1,0] = Fyx/(1000*Lxx*t)

    return stress