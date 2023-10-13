import numpy as np
class Timeseries:
    def __init__(self, t, pos, vel):
        self.ts =   [t]
        self.poss = [pos]
        self.vels = [vel]

    def __add__(self, triple):
        self.ts +=   [triple[0]]
        self.poss += [triple[1]]
        self.vels += [triple[2]]
        return self


class Particle:
    def __init__(self,t,pos,vel,m,ind):
        self.position = np.array(pos)
        self.velocity = np.array(vel)
        self.time     = t
        self.mass     = m
        self.index    = ind

        self.timeseries = Timeseries(self.time, self.position, self.velocity)

    def update_timeseries(self):
        self.timeseries+=[self.time, self.position, self.velocity]

    def update(self, calced):
        self.position,self.velocity,self.time = calced
        self.update_timeseries()