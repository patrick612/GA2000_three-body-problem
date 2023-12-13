import numpy as np
import matplotlib.pyplot as plt
from fractions import Fraction
import time
import numpy as np
from TBPConsts import *
from TBPData import *
from TBPSimulate import *

def init_velo_circular(init_list, omega = False):
    m1 = init_list[0][0]
    r1x = init_list[0][1][0]
    r1y = init_list[0][1][1]
    r1z = init_list[0][1][2]

    m2 = init_list[1][0]
    r2x = init_list[1][1][0]
    r2y = init_list[1][1][1]
    r2z = init_list[1][1][2]

    m3 = init_list[2][0]
    r3x = init_list[2][1][0]
    r3y = init_list[2][1][1]
    r3z = init_list[2][1][2]


    r12 = np.sqrt((r1x - r2x)**2 + (r1y - r2y)**2 + (r1z - r2z)**2)
    r23 = np.sqrt((r3x - r2x)**2 + (r3y - r2y)**2 + (r3z - r2z)**2)
    r1_norm = np.sqrt(r1x**2 + r1y**2 + r1z**2)
    r2_norm = np.sqrt(r2x**2 + r2y**2 + r2z**2)
    r3_norm = np.sqrt(r3x**2 + r3y**2 + r3z**2)
    omega12 = np.sqrt(m2/(r12**2 * r1_norm))
    omega32 = np.sqrt(m2/(r23**2 * r3_norm))


    rot_mat = [[np.cos(np.pi/2), -np.sin(np.pi/2)], [np.sin(np.pi/2), np.cos(np.pi/2)]]
    rot_mat = np.array(rot_mat)
    v1 = omega12 * np.matmul(rot_mat, np.array(init_list[0][1][:-1]))
    v2 = omega12 * np.matmul(rot_mat, np.array(init_list[1][1][:-1]))
    v3 = omega32 * np.matmul(rot_mat, np.array(init_list[2][1][:-1]))
    v1 = list(v1)
    v1.append(0)
    v2 = list(v2)
    v2.append(0)
    v3 = list(v3)
    v3.append(0)

    init_l = init_list.copy()
    init_l[0].append(v1)
    init_l[1].append(v2)
    init_l[2].append(v3)
    if omega == False:
        return init_l
    if omega == True:
        return init_l, omega12, omega32

def TBP_circular(init_list, time_scale, t_f, method = None, spac = None):

    if spac ==None:
        spac = 1
    
    ps1 = ParticleData()
    ps1 += init_list[0]
    ps1 += init_list[1]
    ps1 += init_list[2]
    
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
    
    poss3 = ps1.timeseries.pos(2,spac=spac)
    x3 = [i[0] for i in poss3]
    x3 = np.array(x3)
    y3 = [i[1] for i in poss3]
    y3 = np.array(y3)
    
    p1_pos = [x1, y1]
    p1_pos = np.array(p1_pos)
    p2_pos = [x2, y2]
    p2_pos = np.array(p2_pos)
    p3_pos = [x3, y3]
    p3_pos = np.array(p3_pos)
    
    return p1_pos, p2_pos, p3_pos, ts

def ejected(x3_start, x3_end, steps, m3, period_factor, init_list, pos_data =False):
    step_list = np.linspace(x3_start, x3_end, steps)
    q1_list = []
    q3_list = []
    q1_radii = []
    q3_radii = []
    w_list = []
    time_list = []
    for i in step_list:
        init_l = init_list.copy()
        
        init_l.append([m3, [i, 0, 0]])

        
        init_l, w12, w23 = init_velo_circular(init_l, omega =True)
        if w12 >= w23:
            w = w12
        else:
            w = w23
        
        w_list.append(w)
        period = 2*np.pi/w
        q1, q2, q3, ts = TBP_circular(init_l, period/100, period_factor * period, spac=1)
        
        if len(ts) >= len(time_list):
            time_list = ts
            
        q1_list.append(q1)
        q3_list.append(q3)
        
        q1_r = np.sqrt( np.array(q1[0])**2 + np.array(q1[1])**2)
        q3_r = np.sqrt( np.array(q3[0])**2 + np.array(q3[1])**2)
        
        q1_radii.append(q1_r)
        q3_radii.append(q3_r)
    
    for i in range(len(q1_radii)):
        if len(q1_radii[i]) < len(time_list):
            len_diff = len(time_list) - len(q1_radii[i])
            a = np.concatenate((q1_radii[i], np.zeros(len_diff)))
            q1_radii[i] = a
            
            
    for i in range(len(q3_radii)):
        if len(q3_radii[i]) < len(time_list):
            len_diff = len(time_list) - len(q3_radii[i])
            a = np.concatenate((q3_radii[i], np.zeros(len_diff)))
            q3_radii[i] = a
    
    if pos_data == True:
        return q1_list, q3_list, step_list, time_list
    else:
        return q1_radii, q3_radii, step_list, time_list
    


def r_vs_w(q1_radii, q3_radii, step_list):
    q1_max = []
    q3_max = []
    
    for i in range(len(q1_radii)):
        if len(q1_radii[i]) == 0 or len(q3_radii[i]) == 0:
            st_list = list(step_list)
            st_list.pop(i)
            continue
        max1 = np.max(q1_radii[i])
        q1_max.append(max1)
        
        max3 = np.max(q3_radii[i])
        q3_max.append(max3)
    
    return q1_max, q3_max, st_list
        
