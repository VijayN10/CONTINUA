def numberOfIntegrationPoints(ncoord, nelnodes, elident):
    if (ncoord == 1):
        n = nelnodes
    elif (ncoord == 2):
        if (nelnodes == 3):
            n = 1
        elif (nelnodes == 6):
            n = 3
        elif (nelnodes == 4):
            n = 4
        elif (nelnodes == 8):
            n = 9
    elif (ncoord == 3):
        if (nelnodes == 4):
            n = 1
        elif (nelnodes == 10):
            n = 4
        elif (nelnodes == 8):
            n = 8
        elif (nelnodes == 20):
            n = 27
    return n