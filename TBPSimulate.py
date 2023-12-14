import numpy as np
import time
from TBPConsts import *
from TBPData import *
class Solver:
    def __init__(self,stateobj,frameomega=0.0):
        self.data = stateobj
        self.frameomega=frameomega

    def x(self,state):
        return state[:self.data.n]

    def v(self,state):
        return state[self.data.n:]
    
    def a(self,state):
        xs=self.x(state)
        #Reduced force AKA acceleration
        f = np.zeros((self.data.n,3),dtype=np.float64)
        for i in range(self.data.n):
            for j in range(self.data.n):
                if i!=j:
                    dij = xs[j]-xs[i]
                    f[i]+=G*self.data.masses[j]*dij/(np.linalg.norm(dij))**3
            if self.frameomega:
                vs=self.v(state)
                #Coriolis
                f[i]-=2*self.frameomega*np.array([-vs[i][1],vs[i][0],0])
                #Centrifugal
                f[i]+=self.frameomega**2*np.array([xs[i][0],xs[i][1],0])
        return f

    def ds(self,state):
        return np.concatenate((self.v(state),self.a(state)), axis=0)
    
class EulerFirst(Solver):
    def __call__(self,dt):
        k1 = self.ds(self.data.state)
        return (dt, dt*k1)
    
class RK4(Solver):
    def __call__(self,dt):
        s=self.data.state
        k1 = self.ds(s)
        k2 = self.ds(s+dt/2*k1)
        k3 = self.ds(s+dt/2*k2)
        k4 = self.ds(s+dt*k3)
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
    def __init__(self,stateobj,frameomega=0.0):
        super().__init__(stateobj,frameomega=frameomega)
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
    def __init__(self,stateobj,eps,frameomega=0.0):
        super().__init__(stateobj,frameomega=frameomega)
        self.eps = eps
        self.acoeff = np.array([[1/4,0,0,0,0,0],
                           [3/32,9/32,0,0,0,0],
                           [1932/2197,-7200/2197,7296/2197,0,0,0],
                           [439/216,-8,3680/513,-845/4104,0,0],
                           [-8/27,2,3544/2565,1859/4104,-11/40,0]],
                           dtype=np.float64)

        self.b4th = np.array([25/216,0,1408/2565,2197/4104,-1/5,0],dtype=np.float64)
        self.b5th = np.array([16/135,0,6656/12825,28561/56430,-9/50,2/55],dtype=np.float64)

    def compnext(self,dt):
        s = self.data.state
        k = np.zeros((6,*s.shape))
        k[0] = self.ds(self.data.state)
        for i in range(1,6):
            k[i] = self.ds(s+np.sum(k*self.acoeff[i-1][:,None,None],axis=0)*dt)
        return (np.sum(k*self.b4th[:,None,None],axis=0)*dt,np.sum(k*self.b5th[:,None,None],axis=0)*dt)

    def __call__(self,dt):
        er = 1e70
        while er>self.eps:
            ds4th,ds5th = self.compnext(dt)
            er = np.sum((ds4th-ds5th)**2)**0.5
            #In case error is zero
            if er:
                dt *= 0.9*(self.eps/er)**0.2
        return (dt,ds4th)

    
#Add More methods here

def simulate(particles, solver, dt, tfinal=1.0, updatespac=1000, saveinter = 0, savefile="tseriesdata"):
    i=0
    while particles.time<tfinal:
        tmp=solver(dt)
        particles.time+=dt
        particles.state+=tmp[1]
        dt=tmp[0]
        particles.update_timeseries()
        i+=1
        if i%updatespac==0:
            print(f't = {round(particles.time,5)}, {particles.time/tfinal*100}% done')
            print ("\033[A\033[A")
        if saveinter and i%saveinter==0:
            write_timeseries(particles.timeseries, savefile)
    print(f't = {round(particles.time,5)}, {particles.time/tfinal*100}% done')
    print ("\033[A\033[A")
    write_timeseries(particles.timeseries, savefile)

def realtimesim(particles, solver, dt):
    t  = time.time_ns()
    while time.time_ns()-t<1e9/40:
        tmp=solver(dt)
        particles.time+=dt
        particles.state+=tmp[1]
        dt=tmp[0]
        particles.update_timeseries()