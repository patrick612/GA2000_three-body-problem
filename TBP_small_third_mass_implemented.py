import numpy as np
import matplotlib.pyplot as plt
from fractions import Fraction
import time
import numpy as np
from TBPConsts import *
from TBPData import *
from TBPSimulate import *
from TBP_small_third_mass import *

#try plotting orbits o 3 bodies

init_list = [[3, [200, 0, 0]] , [300, [-2, 0, 0]], [3,[600,0,0]]]
init_list = init_velo_circular(init_list)
period = 2*np.pi * np.sqrt(202**2 * 2/3)
q1, q2, q3, ts = TBP_circular(init_list, period/100, 3*period)
plt.plot(q1[0], q1[1])
plt.plot(q2[0], q2[1])
plt.plot(q3[0], q3[1])
plt.xlabel('x')
plt.ylabel('y')
plt.title('Three Body Orbit')

#try placing the third body far away
init_list = [[3, [200, 0, 0]] , [300, [-2, 0, 0]], [3,[600,0,0]]]
init_list = init_velo_circular(init_list)
period = 2*np.pi * np.sqrt(202**2 * 2/3)
q1, q2, q3, ts = TBP_circular(init_list, period/100, 8*period)
plt.plot(q1[0], q1[1])
plt.plot(q2[0], q2[1])
plt.plot(q3[0], q3[1])
plt.xlabel('x')
plt.ylabel('y')
plt.title('Three Body Orbit')

#try different initial positions and see how radii evolve

init_list = [[3, [200, 0, 0]] , [300, [-2, 0, 0]]]
q1_radii, q3_radii, step_list, time_list = ejected(0, 800, 150, 3, 15, init_list)

for i in range(len(q1_radii)):
    if len(q1_radii[i]) == 0:
        continue
    plt.plot(time_list[:600], q1_radii[i][:600])
    plt.xlabel("Time (G=1)")
    plt.ylabel("Radius (G=1)")
    plt.title("Radi of small mass 1")
    
for i in range(len(q3_radii)):
    if len(q3_radii[i]) == 0:
        continue
    plt.plot(time_list[:600], q3_radii[i][:600])
    plt.xlabel("Time (G=1)")
    plt.ylabel("Radius (G=1)")
    plt.title("Radi of small mass 2")

#plot max radius after 15 period of orbit for different initial conditions

q1_max, q3_max, st_list = r_vs_w((q1_radii), q3_radii, step_list)
plt.plot(st_list, q3_max)
plt.xlabel('Initial Position in x (G=1)')
plt.ylabel('Maximum Radius (G=1)')
plt.title('Maximum radius of small mass 2 within 15 periods')

plt.plot(st_list, q1_max)
plt.xlabel('Initial Position in x (G=1)')
plt.ylabel('Maximum Radius (G=1)')
plt.title('Maximum radius of small mass 1 within 15 periods')   
