import numpy as np
from TBPConsts import *
class Solver:
    def __init__(self,particlelist):
        self.particles=particlelist
    
    def accel(self,n, pos=None):
        f = np.zeros(3)
        for i in range(len(self.particles)):
            if i!=n:
                currp = self.particles[n] if not pos else pos
                pi = self.particles[i]
                r = pi.position-currp.position
                #print(np.linalg.norm(r), r, n, i, pi.position, currp.position)
                f+=G*pi.mass*r/(np.linalg.norm(r))**3
        return f
    
class EulerFirst(Solver):
    def __call__(self,i,dt):
        p = self.particles[i]
        return (p.position+p.velocity*dt, p.velocity+self.accel(i)*dt, p.time+dt)
    
#Add More methods here

def simulate(particles, solver, dt, tfinal=1.0):
    #i=-1
    while particles[0].time<tfinal:
        tmp = []
        for p in range(len(particles)):
            tmp+=[solver(p,dt)]
        for p in range(len(particles)):
            particles[p].update(tmp[p])