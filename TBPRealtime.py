import matplotlib.pyplot as plt
import matplotlib.animation as anim

from TBPSimulate import *
from TBPData import *
from TBPConsts import *

ps = ParticleData()
ps+=[10,[-1/11,0,0],[0,1/11**0.5,0]]
ps+=[1,[10/11,0,0],[0,-10/11**0.5,0]]
ps+=[0.001,[10/11+0.2,0,0],[0,-10/11**0.5*1.9,0]]
#ps+=[0,[1.25608290849,0,0],[0,-1.25608290849*11**0.5,0]]
#ps+=[0,[0.408859,0.866013,0],[2.87224018457,-1.35603189516,0]]

# ps+=[1,[-1.0,0.0,0.0],[0.306893,0.125507,0.0]]
# ps+=[1,[1.0,0.0,0.0],[0.306893,0.125507,0.0]]
# ps+=[1,[0.0,0.0,0.0],[(-2)*0.306893,(-2)*0.125507,0.0]]

print(ps.state)

solv = RK45(ps,0.0001)
fig = plt.figure()
ax = fig.add_subplot()
ax.set_aspect('equal', adjustable='box')

def update(frame):
    realtimesim(ps, solv, 0.0005)
    ax.cla()
    for p in range(ps.state.shape[0]//2):
        #ts = np.array(ps.timeseries.ts)
        poss = ps.timeseries.pos(p,spac=10)
        poss = poss[max(0,len(poss)-1000):]
        xs = [i[0] for i in poss]
        ys = [i[1] for i in poss]
        zs = [i[2] for i in poss]
        ax.scatter(xs,ys, sizes=10/np.arange(len(xs),0,-1)+0.5)
    return ()

ani = anim.FuncAnimation(fig, update, interval=0.001,blit=True)
plt.show()
