import numpy as np
import time
from TBPConsts import *
class Solver:
    def __init__(self,stateobj):
        self.data = stateobj

    def x(self,state):
        return state[:self.data.n]

    def v(self,state):
        return state[self.data.n:]
    
    def a(self,state):
        xs=self.x(state)
        f = np.zeros((self.data.n,3),dtype=np.float64)
        for i in range(self.data.n):
            for j in range(self.data.n):
                if i!=j:
                    dij = xs[j]-xs[i]
                    f[i]+=G*self.data.masses[j]*dij/(np.linalg.norm(dij))**3
        return f

    def ds(self,state):
        return np.concatenate((self.v(state),self.a(state)), axis=0)
    
class EulerFirst(Solver):
    def __call__(self,dt):
        k1 = self.ds(self.data.state)
        return (dt, dt*k1)
    
class RK4(Solver):
    def __call__(self,dt):
        k1 = self.ds(self.data.state)
        k2 = self.ds(self.data.state+dt/2*k1)
        k3 = self.ds(self.data.state+dt/2*k2)
        k4 = self.ds(self.data.state+dt*k3)
        return (dt, dt/6*(k1+2*k2+2*k3+k4))
    
class SympEuler(Solver):
    def __call__(self,dt):
        s = self.data.state
        dv2 = self.a(s)*dt
        v2 = self.v(s)+dv2
        dx2 = v2*dt
        return (dt,np.concatenate((dx2,dv2), axis=0))
    
class Verlet(Solver):
    def __call__(self,dt):
        d1 = d2 = 1/2
        c1 = 1
        c2 = 0

        s = self.data.state
        v1 = self.v(s)
        x1 = self.x(s)
        v2 = v1+d1*self.a(x1)*dt
        x2 = x1+c1*v2*dt
        v3 = v2+d2*self.a(x2)*dt
        x3 = x2+c2*v3*dt

        return (dt,np.concatenate((x3-x1,v3-v1), axis=0))
    
class SympI4(Solver):
    def __init__(self,stateobj):
        super().__init__(stateobj)
        cbt2 = 2**(1/3)
        self.d = [0,1/(2-cbt2),-cbt2/(2-cbt2),1/(2-cbt2)]
        self.c = [self.d[3]/2,(1-cbt2)*self.d[3]/2,(1-cbt2)*self.d[3]/2,self.d[3]/2]

    def __call__(self,dt):

        s = self.data.state
        vi = v = self.v(s)
        xi = x = self.x(s)
        for i in range(4):
            v = v+self.d[i]*self.a(x)*dt
            x = x+self.c[i]*v*dt

        return (dt,np.concatenate((x-xi,v-vi), axis=0))
    
class RK45(Solver):
    def __init__(self,stateobj,eps):
        super().__init__(stateobj)
        self.eps = eps

    def __call__(self,dt):
        return
    
#Add More methods here

def simulate(particles, solver, dt, tfinal=1.0):
    while particles.time<tfinal:
        #print(particles.state)
        tmp=solver(dt)
        particles.time+=dt
        particles.state+=tmp[1]
        dt=tmp[0]
        particles.update_timeseries()

def realtimesim(particles, solver, dt):
    t  = time.time_ns()
    while time.time_ns()-t<1e9/40:
        tmp=solver(dt)
        particles.time+=dt
        particles.state+=tmp[1]
        dt=tmp[0]
        particles.update_timeseries()