import matplotlib.pyplot as plt
import matplotlib.animation as anim

from TBPSimulate import *
from TBPData import *
from TBPConsts import *
from TBPSetups import *

def display(timeseries, dimension=2, spac=1):
    fig = plt.figure('')
    if dimension==2:
        ax = fig.add_subplot()
        ax.set_aspect('equal', adjustable='box')
    elif dimension==3:
        ax = fig.add_subplot(projection='3d')
    for p in range(timeseries.data.n):
        ts = np.array(timeseries.ts)
        poss = timeseries.pos(p,spac=spac)
        xs = [i[0] for i in poss]
        ys = [i[1] for i in poss]
        zs = [i[2] for i in poss]
        if dimension==2:
            ax.scatter(xs,ys, sizes=10/np.arange(len(xs),0,-1)+0.5)
        elif dimension==3:
            ax.scatter(xs,ys,zs, sizes=10/np.arange(len(xs),0,-1)+0.5)
    plt.show()