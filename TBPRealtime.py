import matplotlib.pyplot as plt
import matplotlib.animation as anim

from TBPSimulate import *
from TBPData import *
from TBPConsts import *
from TBPSetups import setup

ps = ParticleData()

for i in setup["TwoEqualOneSmall"]()["init"]:
    ps+=i

print(ps.state)

solv = RK45(ps,eps=0.00001)
#solv = RK45(ps,0.00001,frameomega=np.sqrt(G*(1e2+1)))
fig = plt.figure()
ax = fig.add_subplot()
ax.set_aspect('equal', adjustable='box')

def update(frame):
    realtimesim(ps, solv, 0.00005)
    ax.cla()
    for p in range(ps.state.shape[0]//2):
        #ts = np.array(ps.timeseries.ts)
        poss = ps.timeseries.pos(p,spac=10)
        poss = poss[max(0,len(poss)-1000):]
        xs = [i[0] for i in poss]
        ys = [i[1] for i in poss]
        zs = [i[2] for i in poss]
        ax.scatter(xs,ys, sizes=10/np.arange(len(xs),0,-1)+0.5)
        ax.set_title(f"t={ps.time}")
    return ()

ani = anim.FuncAnimation(fig, update, interval=0.001,blit=True)
plt.show()
