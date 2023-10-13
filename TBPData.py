import numpy as np
class Timeseries:
    def __init__(self):
        self.ts     = []
        self.states = []

    def __add__(self, tstate):
        self.ts     += [tstate[0]]
        self.states += [tstate[1]]
        return self


class ParticleData:
    def __init__(self):
        self.state  = np.zeros((0,3))
        self.masses = []
        self.time   = 0
        self.timeseries = Timeseries()

    def __add__(self,particleparams):
        self.masses+=[particleparams[0]]
        self.state = np.concatenate((self.state, [particleparams[1],particleparams[2]]), axis=0)
        return self

    def update_timeseries(self):
        self.timeseries+=[self.time,np.copy(self.state)]