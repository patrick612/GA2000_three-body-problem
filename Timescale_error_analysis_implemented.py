import numpy as np
import matplotlib.pyplot as plt
import Timescale_error_analysis
from fractions import Fraction


#########RK4
#RK4 5e-5 to 1e-1 of period
time_steplist = []
for i in np.arange(20*2*np.pi * 1/500000, 20*2*np.pi * 1/10, 20*2*np.pi * 1/200):
    time_steplist.append(i)
    
p1err_list, p2err_list = time_step_error(time_steplist, 1, 'RK4')
        

time_steplist_frac = []
for i in range(len(time_steplist)):
    frac = Fraction(np.arange(1/500000, 1/10, 1/200)[i]).limit_denominator()
    time_steplist_frac.append(frac)
    

plt.plot(time_steplist_frac, p1err_list)
plt.xlabel('Time Step (fraction of period)')
plt.ylabel('Error (%)')
plt.title('RK4 Error of different time steps')

#RK4 5e-6 to 1e-4
time_steplist = []
for i in np.arange(20*2*np.pi * 1/5000000, 20*2*np.pi * 1/10000, 20*2*np.pi * 1/100000):
    time_steplist.append(i)
    
p1err_list, p2err_list = time_step_error(time_steplist, 1, 'RK4')
        

time_steplist_frac = []
for i in range(len(time_steplist)):
    frac = Fraction(np.arange(1/5000000, 1/10000,  1/100000)[i]).limit_denominator()
    time_steplist_frac.append(frac)
    

plt.plot(time_steplist_frac, p1err_list)
plt.xlabel('Time Step (fraction of period)')
plt.ylabel('Error (%)')
plt.title('RK4 Error of different time steps')

#########EulerFirst
#EulerFirst 5e-5 to 1e-1 of period
time_steplist = []
for i in np.arange(20*2*np.pi * 1/500000, 20*2*np.pi * 1/10, 20*2*np.pi * 1/200):
    time_steplist.append(i)
    
p1err_list, p2err_list = time_step_error(time_steplist, 1, 'EulerFirst')
        

time_steplist_frac = []
for i in range(len(time_steplist)):
    frac = Fraction(np.arange(1/500000, 1/10, 1/200)[i]).limit_denominator()
    time_steplist_frac.append(frac)
    

plt.plot(time_steplist_frac, p1err_list)
plt.xlabel('Time Step (fraction of period)')
plt.ylabel('Error (%)')
plt.title('EulerFirst Error of different time steps')

#EulerFirst 5e-6 to 1e-4
time_steplist = []
for i in np.arange(20*2*np.pi * 1/5000000, 20*2*np.pi * 1/10000, 20*2*np.pi * 1/100000):
    time_steplist.append(i)
    
p1err_list, p2err_list = time_step_error(time_steplist, 1, 'EulerFirst')
        

time_steplist_frac = []
for i in range(len(time_steplist)):
    frac = Fraction(np.arange(1/5000000, 1/10000,  1/100000)[i]).limit_denominator()
    time_steplist_frac.append(frac)
    

plt.plot(time_steplist_frac, p1err_list)
plt.xlabel('Time Step (fraction of period)')
plt.ylabel('Error (%)')
plt.title('EulerFirst Error of different time steps')

#########SympEuler
#EulerFirst 5e-5 to 1e-1 of period
time_steplist = []
for i in np.arange(20*2*np.pi * 1/500000, 20*2*np.pi * 1/10, 20*2*np.pi * 1/200):
    time_steplist.append(i)
    
p1err_list, p2err_list = time_step_error(time_steplist, 1, 'SympEuler')
        

time_steplist_frac = []
for i in range(len(time_steplist)):
    frac = Fraction(np.arange(1/500000, 1/10, 1/200)[i]).limit_denominator()
    time_steplist_frac.append(frac)
    

plt.plot(time_steplist_frac, p1err_list)
plt.xlabel('Time Step (fraction of period)')
plt.ylabel('Error (%)')
plt.title('SympEuler Error of different time steps')

#SympEuler 5e-6 to 1e-4
time_steplist = []
for i in np.arange(20*2*np.pi * 1/5000000, 20*2*np.pi * 1/10000, 20*2*np.pi * 1/100000):
    time_steplist.append(i)
    
p1err_list, p2err_list = time_step_error(time_steplist, 1, 'SympEuler')
        

time_steplist_frac = []
for i in range(len(time_steplist)):
    frac = Fraction(np.arange(1/5000000, 1/10000,  1/100000)[i]).limit_denominator()
    time_steplist_frac.append(frac)
    

plt.plot(time_steplist_frac, p1err_list)
plt.xlabel('Time Step (fraction of period)')
plt.ylabel('Error (%)')
plt.title('SympEuler Error of different time steps')

#########SympI4
#EulerFirst 5e-5 to 1e-1 of period
time_steplist = []
for i in np.arange(20*2*np.pi * 1/500000, 20*2*np.pi * 1/10, 20*2*np.pi * 1/200):
    time_steplist.append(i)
    
p1err_list, p2err_list = time_step_error(time_steplist, 1, 'SympI4')
        

time_steplist_frac = []
for i in range(len(time_steplist)):
    frac = Fraction(np.arange(1/500000, 1/10, 1/200)[i]).limit_denominator()
    time_steplist_frac.append(frac)
    

plt.plot(time_steplist_frac, p1err_list)
plt.xlabel('Time Step (fraction of period)')
plt.ylabel('Error (%)')
plt.title('SympI4 Error of different time steps')

#SympI4 5e-6 to 1e-4
time_steplist = []
for i in np.arange(20*2*np.pi * 1/5000000, 20*2*np.pi * 1/10000, 20*2*np.pi * 1/100000):
    time_steplist.append(i)
    
p1err_list, p2err_list = time_step_error(time_steplist, 1, 'SympI4')
        

time_steplist_frac = []
for i in range(len(time_steplist)):
    frac = Fraction(np.arange(1/5000000, 1/10000,  1/100000)[i]).limit_denominator()
    time_steplist_frac.append(frac)
    

plt.plot(time_steplist_frac, p1err_list)
plt.xlabel('Time Step (fraction of period)')
plt.ylabel('Error (%)')
plt.title('SympI4 Error of different time steps')
