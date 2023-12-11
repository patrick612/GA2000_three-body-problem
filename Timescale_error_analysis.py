import numpy as np
import matplotlib.pyplot as plt
from TBPConsts import *
from TBPData import *
from TBPSimulate import *



def circular_orbit_err(m1, r1, m2, r2, time_scale, t_f, method = None, spac = None, perc = False):
    r = r1 + r2
    omega = np.sqrt(m1/(r**2 * r2))
    if spac ==None:
        spac = 1


    ps1 = ParticleData()
    ps1 += [m1,[r1,0,0],[0,-1*r1*omega,0]]
    ps1 += [m2, [-r2, 0, 0], [0, r2*omega, 0]]
    
    if method == None or method == 'RK4':
        solve = RK4(ps1)
    elif method == 'EulerFirst':
        solve = EulerFirst(ps1)
    elif method == 'SympEuler':
        solve = SympEuler(ps1)
    elif method == 'SympI4':
        solve = SympI4(ps1)
        
    
    simulate(ps1, solve, time_scale, tfinal=t_f)

    ts = np.array(ps1.timeseries.ts)
    poss1 = ps1.timeseries.pos(0,spac=spac)
    x1 = [i[0] for i in poss1]
    x1 = np.array(x1)
    y1 = [i[1] for i in poss1]
    y1 = np.array(y1)

    poss2 = ps1.timeseries.pos(1,spac=spac)
    x2 = [i[0] for i in poss2]
    x2 = np.array(x2)
    y2 = [i[1] for i in poss2]
    y2 = np.array(y2)


    times = ts
    x1_true = r1 * np.cos(omega * times)
    y1_true = r1 * -1*np.sin(omega * times)
    x1_true = x1_true[::spac]
    y1_true = y1_true[::spac]

    x2_true = -1*r2 * np.cos(omega * times)
    y2_true = r2 * np.sin(omega * times)
    x2_true = x2_true[::spac]
    y2_true = y2_true[::spac]


    x1_err = np.abs(x1_true - x1)
    y1_err = np.abs(y1_true - y1)

    x2_err = np.abs(x2_true - x2)
    y2_err = np.abs(y2_true - y2)

    if perc == False:
        p1_err = np.sqrt(x1_err**2 + y1_err**2)
        p2_err = np.sqrt(x2_err**2 + y2_err**2)
    else:
        p1_err = np.sqrt(x1_err**2 + y1_err**2)/np.sqrt(x1_true**2 + y1_true**2)
        p2_err = np.sqrt(x2_err**2 + y2_err**2)/np.sqrt(x2_true**2 + y1_true**2)

    return p1_err, p2_err, times[::spac]
    
def time_step_error(time_step_list, t_f, method):
    m1 = 10
    r1 = 10
    m2 = 10
    r2 = 10
    r = r1 + r2
    omega = np.sqrt(m1/(r**2 * r2))
    period = np.pi*2/omega
    t_f = t_f * period
    
    p1_err_list = []
    p2_err_list = []
    for step_size in time_step_list:
        p1err, p2err, times = circular_orbit_err(m1, r1, m2, r2, step_size, t_f, perc=True, method=method, spac = 1)
        p1_err_list.append(p1err[-1])
        p2_err_list.append(p2err[-1])
    
    return p1_err_list, p2_err_list
