from TBPSimulate import *
from TBPData import *
from TBPVisualize import *
import matplotlib.pyplot as plt


def runsetup(setupname,setupargs,solverclass,solverargs,initdt,tfinal,savefile):
    ps = ParticleData()
    s = setup[setupname](*setupargs)
    for i in s["init"]:
        ps+=i
    solv = solverclass(ps,*solverargs,frameomega=s["omega"])
    simulate(ps, solv, initdt, tfinal=tfinal, savefile=savefile)

runsetup("L4",[],RK45,[0.0001],0.0001,12*np.pi,"L4data")
ts = load_timeseries("L4data")
display(ts, dimension=2)