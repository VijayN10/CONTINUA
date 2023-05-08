def plotMesh(coords, ncoord, nnode, connect, nelem, elident, nelnodes, color):
    f2D_3 = np.array([1,2,3])
    f2D_4 = np.array([1,2,3,4])
    f2D_6 = np.array([1,4,2,5,3,6])
    f2D_8 = np.array([1,5,2,6,3,7,4,8])

    fig, ax = plt.subplots()
    for lmn in range(nelem):
        x = []
        for i in range(nelnodes):
            x.append(coords[:, connect[i, lmn]-1])
        x = np.array(x)
        ax.scatter(x[:, 0], x[:, 1], color='r', edgecolors='none')
        if nelnodes == 3:
            ax.add_patch(Polygon(x[f2D_3-1], facecolor='none', edgecolor=color))
        elif nelnodes == 4:
            ax.add_patch(Polygon(x[f2D_4-1], facecolor='none', edgecolor=color))
        elif nelnodes == 6:
            ax.add_patch(Polygon(x[f2D_6-1], facecolor='none', edgecolor=color))
        elif nelnodes == 8 or nelnodes == 9:
            ax.add_patch(Polygon(x[f2D_8-1], facecolor='none', edgecolor=color))
    ax.autoscale(True)
    ax.set_aspect('equal')
    plt.show()