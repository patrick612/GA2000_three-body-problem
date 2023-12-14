import numpy as np
from TBPConsts import *
class Timeseries:
    #Inputs: 
    #   dataobj, type ParticleData, some particle data
    #Outputs: 
    #   None
    #Creates Timeseries object
    def __init__(self,dataobj):
        self.ts     = []
        self.states = []
        self.data = dataobj

    #Inputs: 
    #   tstate, type List, a list in the format [mass, position, velocity]
    #Outputs: 
    #   self
    #Add a particle to those being tracked in the timeseries with initial params
    def __add__(self, tstate):
        self.ts     += [tstate[0]]
        self.states += [tstate[1]]
        return self
    
    #Inputs: 
    #   n, Integer, which particle to get the position data for, None if all
    #   spac, Integer, spacing in which to slice the data, gets a chunk of it
    #Outputs: 
    #   type List, list of positions. If n is None, list of list, list[i] is the ith particle's data
    #Gives position data
    def pos(self, n=None, spac = 100):
        if n==None:
            return [[j[i] for j in self.states[::spac]] for i in range(len(self.states[0])//2)]
        else:
            return [j[n] for j in self.states[::spac]]
    
    #Inputs: 
    #   n, Integer, which particle to get the velocity data for, None if all
    #   spac, Integer, spacing in which to slice the data, gets a chunk of it
    #Outputs: 
    #   type List, list of velocities. If n is None, list of list, list[i] is the ith particle's data
    #Gives velocity data
    def vel(self, n=None, spac = 100):
        if n==None:
            return [[j[j.shape[0]//2+i] for j in self.states[::spac]] for i in range(len(self.states[0])//2)]
        else:
            return [j[j.shape[0]//2+n] for j in self.states[::spac]]

    #Inputs: 
    #   n, Integer, which particle to get the potential energy data for, None if all
    #   spac, Integer, spacing in which to slice the data, gets a chunk of it
    #Outputs: 
    #   type List, list of potentials. If n is None, list of list, list[i] is the ith particle's data
    #Gives potential energy data  
    def U(self, n=None, spac = 100):
        U=[]
        for s in self.states[::spac]:
            xs = s[:self.data.n]
            u = 0
            if n==None: 
                for i in range(self.data.n):
                    for j in range(i+1,self.data.n):
                        dij = xs[j]-xs[i]
                        u-=G*self.data.masses[j]*self.data.masses[i]/(np.linalg.norm(dij))
            else:
                for i in range(self.data.n):
                    din = xs[n]-xs[i]
                    u-=G*self.data.masses[n]*self.data.masses[i]/(np.linalg.norm(din))
            U+=[u]
        return U
        
    #Inputs: 
    #   n, Integer, which particle to get the kinetic energy data for, None if all
    #   spac, Integer, spacing in which to slice the data, gets a chunk of it
    #Outputs: 
    #   type List, list of kinetic energies. If n is None, list of list, list[i] is the ith particle's data
    #Gives kinetic energy data 
    def T(self, n=None, spac = 100):
        K=[]
        for s in self.states[::spac]:
            vs = s[self.data.n:]
            if n==None:
                k = sum([0.5*np.sum(vs[i]*vs[i])*self.data.masses[i] for i in range(self.data.n)])
            else:
                k = 0.5*np.sum(vs[n]*vs[n])*self.data.masses[n]
            K+=[k]
        return K
    
    #Inputs: 
    #   n, Integer, which particle to get the total energy data for, None if all
    #   spac, Integer, spacing in which to slice the data, gets a chunk of it
    #Outputs: 
    #   type List, list of total energies. If n is None, list of list, list[i] is the ith particle's data
    #Gives total energy data 
    def energy(self, n=None, spac = 100):
        U=self.U(n,spac=spac)
        T=self.T(n,spac=spac)
        return [U[i]+T[i] for i in range(0,len(U),spac)]




class ParticleData:
    def __init__(self):
        self.state  = np.zeros((0,3))
        self.masses = []
        self.time   = 0
        self.timeseries = Timeseries(self)
        self.n = 0

    def __add__(self,particleparams):
        self.masses+=[particleparams[0]]
        if self.n==0:
            self.state = np.array([particleparams[1],particleparams[2]])
        else:
            self.state = np.concatenate((self.state[:self.n], [particleparams[1]], self.state[self.n:],[particleparams[2]]), axis=0)
        self.n+=1
        return self

    def update_timeseries(self):
        self.timeseries+=[self.time,np.copy(self.state)]


def write_timeseries(tsobj, file):
    np.savez(file,ts = np.array(tsobj.ts),states = np.stack(tsobj.states),masses = np.array(tsobj.data.masses))

def load_timeseries(file):
    dat = np.load(file+".npz")
    tmpdat = ParticleData()
    tmpdat.masses = list(dat["masses"])
    tmpdat.n = len(dat["masses"])
    tseries = Timeseries(tmpdat)
    tseries.ts = list(dat["ts"])
    tseries.states = list(dat["states"])
    return tseries
