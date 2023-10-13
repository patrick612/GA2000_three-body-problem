import matplotlib.pyplot as plt
import matplotlib.animation as anim

from TBPSimulate import *
from TBPData import *
from TBPConsts import *

ps = ParticleData()
ps+=[1,[-1,0,0],[0,0.5,0]]
ps+=[1,[1,0,0],[0,-0.5,0]]

solv = RK4(ps)

simulate(ps, solv, 0.001, tfinal=2*np.pi)

for p in range(ps.state.shape[0]//2):
    ts = ps.timeseries.ts
    xs = [i[2*p,0] for i in ps.timeseries.states]
    ys = [i[2*p,1] for i in ps.timeseries.states]
    plt.scatter(xs,ys)

plt.axis('equal')
plt.show()