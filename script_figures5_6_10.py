#
#Figures 5, 6 and 10 of the paper can be reproduced using this script
#

from scipy import interpolate
import numpy as N 
import matplotlib.pyplot as plt
import os 
import netCDF4

import numpy
from scipy.stats import pearsonr
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt


L_vap                = 2.501e6     #J/kg
C_p                  = 1005.0      #J/kgK
R                    = 287.058 
kappa                = R/C_p    
gra                  = 9.8
L_sub=2.85e6

time_10min_CTRL_f_ens_mean_R=np.load('time_10min_CTRL_f_ens_mean.npy')
time_10min_MHALF_f_ens_mean_R=np.load('time_10min_S65L200_ens_mean.npy')
time_10min_PHALF_f_ens_mean_R=np.load('time_10min_S195L600_ens_mean.npy')
time_10min_RTHQV_ens_mean_R=np.load('time_10min_RTHQV_ens_mean.npy')

nr=36+93
time_10min_CTRL_shift0=N.zeros((nr))
time_10min_PHALF_shift0=N.zeros((nr))
time_10min_MHALF_shift0=N.zeros((nr))
time_10min_RTHQV_shift0=N.zeros((nr))

# t_0=time after triggering
time_10min_PHALF_shift0[45:100]=time_10min_PHALF_f_ens_mean_R[45:100]- 2.5
time_10min_MHALF_shift0[48:100]=time_10min_MHALF_f_ens_mean_R[48:100] -3.25
time_10min_CTRL_shift0[46:100]=time_10min_CTRL_f_ens_mean_R[46:100]- 2.75
time_10min_RTHQV_shift0[46:100]=time_10min_RTHQV_ens_mean_R[46:100]-2.75


nlag=55
lag=N.zeros((nlag))
zero_lag=N.zeros((nlag))
lag[0]=0
for j in N.arange(nlag-1):
     lag[j+1]=lag[j]-0.25
zero_lag[:]=0.0

#
#values of P[R(A,to)] as a function of time
#
proR_lag_ens_mean_PHALF = np.load('proR_lag_ens_mean_PHALF.npy')
proR_lag_ens_mean_MHALF = np.load('proR_lag_ens_mean_MHALF.npy')
proR_lag_ens_mean_RTHQV = np.load('proR_lag_ens_mean_RTHQV.npy')
proR_lag_ens_mean_CTRL = np.load('proR_lag_ens_mean_CTRL.npy')

#
#values of M[R(A,t0,Dt)] as a function of time lag Dt and to: (Dt, tO)
#
proRR_ProR2t_lag_dt_RTHQV=np.load('proRR_ProR2t_lag_dt_RTHQV.npy')
proRR_ProR2t_lag_dt_MHALF=np.load('proRR_ProR2t_lag_dt_MHALF.npy')
proRR_ProR2t_lag_dt_PHALF=np.load('proRR_ProR2t_lag_dt_PHALF.npy')

proRR_ProR2t_lag_dt_CTRL=np.load('proRR_ProR2t_lag_dt_CTRL.npy')
proRR_ProR2t_lag_dt_CTRL_xy=np.load('proRR_ProR2t_lag_dt_CTRL_xy.npy')
proRR_ProR2t_lag_dt_CTRL_22=np.load('proRR_ProR2t_lag_dt_CTRL_22.npy')
proRR_ProR2t_lag_dt_CTRL_1010=np.load('proRR_ProR2t_lag_dt_CTRL_1010.npy')
proRR_ProR2t_lag_dt_CTRL_1515=np.load('proRR_ProR2t_lag_dt_CTRL_1515.npy')
proRR_ProR2t_lag_dt_CTRL_2525=np.load('proRR_ProR2t_lag_dt_CTRL_2525.npy')
proRR_ProR2t_lag_dt_CTRL_5050=np.load('proRR_ProR2t_lag_dt_CTRL_5050.npy')








fig = plt.figure(figsize=(20,13))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
plt.subplot(231)
plt.plot(time_10min_CTRL_shift0[38:107], proR_lag_ens_mean_CTRL[1:70], marker='D', color="black", label="C")
plt.plot(time_10min_PHALF_shift0[38:107], proR_lag_ens_mean_PHALF[1:70], marker='D', color="red", label="S")
plt.plot(time_10min_MHALF_shift0[38:107], proR_lag_ens_mean_MHALF[1:70], marker='D', color="blue", label="W")
plt.legend()

plt.ylabel(' P[R(A,t$_0$)]', fontsize='x-large')
plt.xlabel('Time after triggering (hour)', fontsize='x-large')
plt.xlim(0,12)
plt.ylim(0,0.8)
axes = plt.gca()
axes.set_xlim([0, 12])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10,11,12])

axes = plt.gca()
axes.set_ylim([0, 0.8])
axes.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
axes.set_yticklabels([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.7, r'a)', fontdict=font, color="black", fontsize='x-large')


plt.subplot(232)
plt.plot(-lag[1:7], proRR_ProR2t_lag_dt_CTRL[1:7,15], marker='D', color="black", label='Ctrl')
plt.plot(-lag[1:7], proRR_ProR2t_lag_dt_PHALF[1:7,14], marker='D', color="red", label='Strong')
plt.plot(-lag[1:7], proRR_ProR2t_lag_dt_MHALF[1:7,17], marker='D', color="blue", label='Weak')


plt.plot(-lag[0:12], zero_lag[0:12],linestyle=":",lw=2, color="black")

plt.ylabel('M(A,t$_0$,$\Delta$t)', fontsize='x-large')
plt.xlabel('$\Delta$t (hour)', fontsize='x-large')
plt.xlim(0,3)
plt.ylim(-0.1,0.15)


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.3, 0.13, r'b) t$_0$=1.5 h', fontdict=font, color="black", fontsize='x-large')


plt.subplot(233)
plt.plot(-lag[1:9], proRR_ProR2t_lag_dt_CTRL[1:9,17], marker='D', color="black", label='Ctrl')
plt.plot(-lag[1:9], proRR_ProR2t_lag_dt_PHALF[1:9,16], marker='D', color="red", label='Strong')
plt.plot(-lag[1:9], proRR_ProR2t_lag_dt_MHALF[1:9,19], marker='D', color="blue", label='Weak')
plt.plot(-lag[0:33], zero_lag[0:33],linestyle=":",lw=2, color="black")

plt.ylabel('M(A,t$_0$,$\Delta$t)', fontsize='x-large')
plt.xlabel('$\Delta$t (hour)', fontsize='x-large')
plt.xlim(0,4)
plt.ylim(-0.1,0.15)


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.3, 0.13, r'c) t$_0$=2 h', fontdict=font, color="black", fontsize='x-large')


plt.subplot(234)
plt.plot(-lag[1:13], proRR_ProR2t_lag_dt_CTRL[1:13,21], marker='D', color="black", label='C')
plt.plot(-lag[1:13], proRR_ProR2t_lag_dt_PHALF[1:13,20], marker='D', color="red", label='S')
plt.plot(-lag[1:10], proRR_ProR2t_lag_dt_MHALF[1:10,23], marker='D', color="blue", label='W')
plt.plot(-lag[0:33], zero_lag[0:33],linestyle=":",lw=2, color="black")

plt.ylabel('M(A,t$_0$,$\Delta$t)', fontsize='x-large')
plt.xlabel('$\Delta$t (hour)', fontsize='x-large')
plt.xlim(0,4)
plt.ylim(-0.1,0.15)


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.3, 0.13, r'd) t$_0$=3 h', fontdict=font, color="black", fontsize='x-large')


plt.subplot(2,3,5)
plt.plot(-lag[1:19], proRR_ProR2t_lag_dt_CTRL[1:19,27], marker='D', color="black", label='Ctrl')
plt.plot(-lag[1:19], proRR_ProR2t_lag_dt_PHALF[1:19,26], marker='D', color="red", label='Strong')
plt.plot(-lag[1:19], proRR_ProR2t_lag_dt_MHALF[1:19,29], marker='D', color="blue", label='Weak')
plt.plot(-lag[0:33], zero_lag[0:33],linestyle=":",lw=2, color="black")

plt.ylabel('M(A,t$_0$,$\Delta$t)', fontsize='x-large')
plt.xlabel('$\Delta$t (hour)')
plt.xlim(0,6)
plt.ylim(-0.1,0.15)


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.13, r'e) t$_0$=4.5 h', fontdict=font, color="black", fontsize='x-large')



plt.subplot(2,3,6)
plt.plot(-lag[1:33], proRR_ProR2t_lag_dt_CTRL[1:33,41], marker='D', color="black", label='Ctrl')
plt.plot(-lag[1:33], proRR_ProR2t_lag_dt_PHALF[1:33,40], marker='D', color="red", label='Strong')
plt.plot(-lag[1:33], proRR_ProR2t_lag_dt_MHALF[1:33,43], marker='D', color="blue", label='Weak')
plt.plot(-lag[0:33], zero_lag[0:33],linestyle=":",lw=2, color="black")

plt.ylabel('M(A,t$_0$,$\Delta$t)', fontsize='x-large')
plt.xlabel('$\Delta$t (hour)', fontsize='x-large')
plt.xlim(0,8)
plt.ylim(-0.1,0.15)


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.13, r'f) t$_0$=8 h', fontdict=font, color="black", fontsize='x-large')
plt.savefig('Figure6.png')



fig = plt.figure(figsize=(20,13))
ax = fig.add_axes([0.1,0.1,0.8,0.8])

plt.subplot(221)
plt.plot(time_10min_CTRL_shift0[38:107], proR_lag_ens_mean_CTRL[1:70], marker='D', color="black", label="C")
plt.ylabel(' P[R(A,t$_0$)]', fontsize='x-large')
plt.xlabel('Time after triggering (hour)', fontsize='x-large')
plt.xlim(0,12)
plt.ylim(0,0.8)
axes = plt.gca()
axes.set_xlim([0, 12])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10,11,12])

axes = plt.gca()
axes.set_ylim([0, 0.8])
axes.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
axes.set_yticklabels([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.7, r'a)', fontdict=font, color="black", fontsize='x-large')


plt.subplot(2,2,2)
plt.plot(-lag[1:6], proRR_ProR2t_lag_dt_CTRL[1:6,15], marker='D', color="black", label='t$_0$=1.5 h')
plt.plot(-lag[1:6], proRR_ProR2t_lag_dt_CTRL[1:6,17], marker='D', color="blue", label='t$_0$=2 h')
plt.plot(-lag[1:9],proRR_ProR2t_lag_dt_CTRL[1:9,18] , marker='D', color="cyan", label='t$_0$=2.25 h')
plt.plot(-lag[1:14], proRR_ProR2t_lag_dt_CTRL[1:14,22], marker='D', color="green", label='t$_0$=3.25 h')
plt.plot(-lag[1:21], proRR_ProR2t_lag_dt_CTRL[1:21,29], marker='D', color="purple", label='t$_0$=5 h')
plt.plot(-lag[1:24], proRR_ProR2t_lag_dt_CTRL[1:24,32], marker='D', color="grey", label='t$_0$=5.75 h')
plt.plot(-lag[1:25], proRR_ProR2t_lag_dt_CTRL[1:25,33], marker='D', color="red", label='t$_0$=6 h')
plt.plot(-lag[1:33], proRR_ProR2t_lag_dt_CTRL[1:33,45], marker='D', color="brown", label='t$_0$=9 h')
plt.legend( fontsize='14')
plt.plot(-lag[0:33], zero_lag[0:33],linestyle=":",lw=2, color="black")

plt.ylabel('M(A,t$_0$,$\Delta$t)', fontsize='x-large')
plt.xlabel('$\Delta$t (hour)', fontsize='x-large')
plt.xlim(0,8)
plt.ylim(-0.1,0.15)


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.13, r'b)', fontdict=font, color="black", fontsize='x-large')

plt.subplot(223)
plt.plot(-lag[1:13], proRR_ProR2t_lag_dt_CTRL_xy[1:13,21], marker='D', color="blue", label='A=0.2X0.2 km$^2$')
plt.plot(-lag[1:13], proRR_ProR2t_lag_dt_CTRL_22[1:13,21], marker='D', color="cyan", label='A=2X2 km$^2$')
plt.plot(-lag[1:13], proRR_ProR2t_lag_dt_CTRL[1:13,21], marker='D', color="black", label='A=4X4 km$^2$')
plt.plot(-lag[1:13], proRR_ProR2t_lag_dt_CTRL_1010[1:13,21], marker='D', color="green", label='A=10X10 km$^2$')
plt.plot(-lag[1:13], proRR_ProR2t_lag_dt_CTRL_1515[1:13,21], marker='D', color="purple", label='A=15X15 km$^2$')
plt.plot(-lag[1:13], proRR_ProR2t_lag_dt_CTRL_2525[1:13,21], marker='D', color="red", label='A=25X25 km$^2$')
plt.plot(-lag[1:13], proRR_ProR2t_lag_dt_CTRL_5050[1:13,21], marker='D', color="brown", label='A=50X50 km$^2$')
plt.legend( fontsize='14')
plt.plot(-lag[0:33], zero_lag[0:33],linestyle=":",lw=2, color="black")

plt.ylabel('M(A,t$_0$,$\Delta$t)', fontsize='x-large')
plt.xlabel('$\Delta$t (hour)')
plt.xlim(0,8)
plt.ylim(-0.1,0.15)


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.13, r'c) t$_0$=3 h', fontdict=font, color="black", fontsize='x-large')


plt.subplot(2,2,4)
plt.plot(-lag[1:25], proRR_ProR2t_lag_dt_CTRL_xy[1:25,33], marker='D', color="blue", label='A=0.2X0.2 km$^2$')
plt.plot(-lag[1:25], proRR_ProR2t_lag_dt_CTRL_22[1:25,33], marker='D', color="cyan", label='A=2X2 km$^2$')
plt.plot(-lag[1:25], proRR_ProR2t_lag_dt_CTRL[1:25,33], marker='D', color="black", label='A=4X4 km$^2$')
plt.plot(-lag[1:25], proRR_ProR2t_lag_dt_CTRL_1010[1:25,33], marker='D', color="green", label='A=10X10 km$^2$')
plt.plot(-lag[1:25], proRR_ProR2t_lag_dt_CTRL_1515[1:25,33], marker='D', color="purple", label='A=15X15 km$^2$')
plt.plot(-lag[1:25], proRR_ProR2t_lag_dt_CTRL_2525[1:25,33], marker='D', color="red", label='A=25X25 km$^2$')
plt.plot(-lag[1:25], proRR_ProR2t_lag_dt_CTRL_5050[1:25,33], marker='D', color="brown", label='A=50X50 km$^2$')
plt.legend( fontsize='14')
plt.plot(-lag[0:33], zero_lag[0:33],linestyle=":",lw=2, color="black")

plt.ylabel('M(A,t$_0$,$\Delta$t)', fontsize='x-large')
plt.xlabel('$\Delta$t (hour)')
plt.xlim(0,8)
plt.ylim(-0.1,0.15)


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.13, r'd) t$_0$=6 h', fontdict=font, color="black", fontsize='x-large')

plt.savefig('Figure5.png')



fig = plt.figure(figsize=(20,13))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
plt.subplot(231)
plt.plot(time_10min_CTRL_shift0[38:107], proR_lag_ens_mean_CTRL[1:70], marker='D', color="black", label="C")
plt.plot(time_10min_RTHQV_shift0[38:107], proR_lag_ens_mean_RTHQV[1:70], marker='D', color="red", label="H")

plt.legend()

plt.ylabel(' P[R(A,t$_0$)]', fontsize='x-large')
plt.xlabel('Time (hour)', fontsize='x-large')
plt.xlim(0,12)
plt.ylim(0,0.9)
axes = plt.gca()
axes.set_xlim([0, 12])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10,11,12])

axes = plt.gca()
axes.set_ylim([0, 0.9])
axes.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9])
axes.set_yticklabels([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9])


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.8, r'a)', fontdict=font, color="black", fontsize='x-large')

plt.subplot(232)
plt.plot(-lag[1:8], proRR_ProR2t_lag_dt_CTRL[1:8,16], marker='D', color="black", label='Ctrl')
plt.plot(-lag[1:8], proRR_ProR2t_lag_dt_RTHQV[1:8,16], marker='D', color="red", label='HP')

plt.plot(-lag[0:33], zero_lag[0:33],linestyle=":",lw=2, color="black")

plt.ylabel('M(A,t$_0$,$\Delta$t)', fontsize='x-large')
plt.xlabel('$\Delta$t (hour)', fontsize='x-large')
plt.xlim(0,3)
plt.ylim(-0.1,0.15)


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.2, 0.13, r'b) t$_0$=1.75 h', fontdict=font, color="black", fontsize='x-large')

plt.subplot(233)
plt.plot(-lag[1:15], proRR_ProR2t_lag_dt_CTRL[1:15,23], marker='D', color="black", label='Ctrl')
plt.plot(-lag[1:15], proRR_ProR2t_lag_dt_RTHQV[1:15,23], marker='D', color="red", label='HP')

plt.plot(-lag[0:33], zero_lag[0:33],linestyle=":",lw=2, color="black")

plt.ylabel('M(A,t$_0$,$\Delta$t)', fontsize='x-large')
plt.xlabel('$\Delta$t (hour)', fontsize='x-large')
plt.xlim(0,4)
plt.ylim(-0.1,0.15)


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.32, 0.13, r'c) t$_0$=3.5 h', fontdict=font, color="black", fontsize='x-large')

plt.subplot(234)
plt.plot(-lag[1:21], proRR_ProR2t_lag_dt_CTRL[1:21,29], marker='D', color="black", label='Ctrl')
plt.plot(-lag[1:21], proRR_ProR2t_lag_dt_RTHQV[1:21,29], marker='D', color="red", label='HP')

plt.plot(-lag[0:33], zero_lag[0:33],linestyle=":",lw=2, color="black")

plt.ylabel('M(A,t$_0$,$\Delta$t)', fontsize='x-large')
plt.xlabel('$\Delta$t (hour)', fontsize='x-large')
plt.xlim(0,8)
plt.ylim(-0.1,0.15)


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.13, r'd) t$_0$=5 h', fontdict=font, color="black", fontsize='x-large')


plt.subplot(235)
plt.plot(-lag[1:24], proRR_ProR2t_lag_dt_CTRL[1:24,33], marker='D', color="black", label='Ctrl')
plt.plot(-lag[1:24], proRR_ProR2t_lag_dt_RTHQV[1:24,33], marker='D', color="red", label='HP')

plt.plot(-lag[0:33], zero_lag[0:33],linestyle=":",lw=2, color="black")

plt.ylabel('M(A,t$_0$,$\Delta$t)', fontsize='x-large')
plt.xlabel('$\Delta$t (hour)', fontsize='x-large')
plt.xlim(0,8)
plt.ylim(-0.1,0.15)


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.13, r'e) t$_0$=6 h', fontdict=font, color="black", fontsize='x-large')


plt.subplot(236)
plt.plot(-lag[1:33], proRR_ProR2t_lag_dt_CTRL[1:33,41], marker='D', color="black", label='Ctrl')
plt.plot(-lag[1:33], proRR_ProR2t_lag_dt_RTHQV[1:33,41], marker='D', color="red", label='HP')

plt.plot(-lag[0:33], zero_lag[0:33],linestyle=":",lw=2, color="black")

plt.ylabel('M(A,t$_0$,$\Delta$t)', fontsize='x-large')
plt.xlabel('$\Delta$t (hour)', fontsize='x-large')
plt.xlim(0,8)
plt.ylim(-0.1,0.15)


font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.13, r'f) t$_0$=8 h', fontdict=font, color="black", fontsize='x-large')

plt.savefig('Figure10.png')


