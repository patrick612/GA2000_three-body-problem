import matplotlib.pyplot as plt
import matplotlib.animation as anim

from TBPSimulate import *
from TBPData import *
from TBPConsts import *

ps = [Particle(0,[-1,0,0],[0,0.5,0],1,0),Particle(0,[1,0,0],[0,-0.5,0],1,1),Particle(0,[1.1,0,0],[0.1,-1/0.1**0.5,0],0.00001,2)]
solv = EulerFirst(ps)

simulate(ps, solv, 0.001, tfinal=6*np.pi)

for p in ps:
    ts = p.timeseries.ts
    xs = [q[0] for q in p.timeseries.poss]
    ys = [q[1] for q in p.timeseries.poss]
    plt.scatter(xs,ys, sizes=ts)

plt.axis('equal')
plt.show()