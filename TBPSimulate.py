import numpy as np
from TBPConsts import *
class Solver:
    def __init__(self,stateobj):
        self.data = stateobj

    def ds(self,state):
        f = np.zeros(state.shape)
        for p in range(state.shape[0]//2):
            f[2*p]=state[2*p+1]
            for q in range(state.shape[0]//2):
                if p!=q:
                    dqp = state[2*q]-state[2*p]
                    f[2*p+1]+=G*self.data.masses[q]*dqp/(np.linalg.norm(dqp))**3
        return f
    
class EulerFirst(Solver):
    def __call__(self,dt):
        k1 = self.ds(self.data.state)
        return dt*k1
    
class RK4(Solver):
    def __call__(self,dt):
        k1 = self.ds(self.data.state)
        k2 = self.ds(self.data.state+dt/2*k1)
        k3 = self.ds(self.data.state+dt/2*k2)
        k4 = self.ds(self.data.state+dt*k3)
        return dt/6*(k1+2*k2+2*k3+k4)
    
#Add More methods here

def simulate(particles, solver, dt, tfinal=1.0):
    while particles.time<tfinal:
        particles.state+=solver(dt)
        particles.time+=dt
        particles.update_timeseries()