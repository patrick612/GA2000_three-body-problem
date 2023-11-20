import matplotlib.pyplot as plt
import matplotlib.animation as anim

from TBPSimulate import *
from TBPData import *
from TBPConsts import *

ps = ParticleData()
ps+=[10,[-1/11,0,0],[0,1/11**0.5,0]]
ps+=[1,[10/11,0,0],[0,-10/11**0.5,0]]
ps+=[0,[1.25608290849,0,0],[0,-1.25608290849*11**0.5,0]]
ps+=[0,[0.408859,0.866013,0],[2.87224018457,-1.35603189516,0]]

# ps+=[1,[-1.0,0.0,0.0],[0.306893,0.125507,0.0]]
# ps+=[1,[1.0,0.0,0.0],[0.306893,0.125507,0.0]]
# ps+=[1,[0.0,0.0,0.0],[(-2)*0.306893,(-2)*0.125507,0.0]]

print(ps.state)

solv = SympI4(ps)

simulate(ps, solv, 0.0001, tfinal=np.pi*6)

for p in range(ps.state.shape[0]//2):
    ts = np.array(ps.timeseries.ts)
    poss = ps.timeseries.pos(p,spac=30)
    xs = [i[0] for i in poss]
    ys = [i[1] for i in poss]
    plt.scatter(xs,ys, sizes=ts/max(ts))


plt.axis('equal')
plt.show()

E = ps.timeseries.energy(spac=20)
E = ps.timeseries.energy(spac=20)
plt.plot(E)
plt.show()