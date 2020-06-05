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



L_vap                = 2.501e6     #J/kg
C_p                  = 1005.0      #J/kgK
R                    = 287.058 
kappa                = R/C_p    
gra                  = 9.8
L_sub=2.85e6 #2.834e6
nz=99
ny=512
nx=512
dx=200
dy=200
length_x=N.zeros(nx)
length_y=N.zeros(ny)

length_x[0]=00.
length_y[0]=00.
for j in N.arange(ny-1):
     length_x[j+1]=length_x[j]+dx
     length_y[j+1]=length_y[j]+dy

length_x=length_x/1000.
length_y=length_y/1000.

hein = np.load('hein.npy')

#2D surface precipitation field
precip_xy = np.load('surf_precip_xy_time_CTRL.npy')
precip_xy = np.load('surf_precip_xy_time_S65L200.npy')
precip_xy = np.load('surf_precip_xy_time_S195L600.npy')

precip_xy_rthqvd2 = np.load('surf_precip_xy_time_RTHQVD2.npy')
precip_xy_rthqvd3 = np.load('surf_precip_xy_time_RTHQVD3.npy')
precip_xy_rthqvd4 = np.load('surf_precip_xy_time_RTHQVD4.npy')
precip_xy_rthqvd5 = np.load('surf_precip_xy_time_RTHQVD5.npy')
precip_xy_rthqvd6 = np.load('surf_precip_xy_time_RTHQVD6.npy')
precip_xy_rthqvd7 = np.load('surf_precip_xy_time_RTHQVD7.npy')
precip_xy_rthqvd8 = np.load('surf_precip_xy_time_RTHQVD8.npy')
precip_xy_rthqvd9 = np.load('surf_precip_xy_time_RTHQVD9.npy')
precip_xy_rthqvd10 = np.load('surf_precip_xy_time_RTHQVD10.npy')


n_30min=93
ntr=1
MRR_thr=N.zeros((ntr))
Surf_precip_threshold_Day_mean=N.zeros((ntr))

precip_xy_d2=N.zeros((n_30min,nx,ny))
precip_xy_d3=N.zeros((n_30min,nx,ny))
precip_xy_d4=N.zeros((n_30min,nx,ny))
precip_xy_d5=N.zeros((n_30min,nx,ny))
precip_xy_d6=N.zeros((n_30min,nx,ny))
precip_xy_d7=N.zeros((n_30min,nx,ny))
precip_xy_d8=N.zeros((n_30min,nx,ny))
precip_xy_d9=N.zeros((n_30min,nx,ny))
precip_xy_d10=N.zeros((n_30min,nx,ny))

hour=N.zeros((n_30min))

#control simulation
precip_xy_d2[0:92,:,:]=3600.*precip_xy[90:182,:,:]
precip_xy_d3[0:92,:,:]=3600.*precip_xy[182:274,:,:]
precip_xy_d4[0:92,:,:]=3600.*precip_xy[273:365,:,:]
precip_xy_d5[0:92,:,:]=3600.*precip_xy[363:455,:,:]
precip_xy_d6[0:92,:,:]=3600.*precip_xy[455:547,:,:]
precip_xy_d7[0:92,:,:]=3600.*precip_xy[547:639,:,:]
precip_xy_d8[0:92,:,:]=3600.*precip_xy[637:729,:,:]
precip_xy_d9[0:92,:,:]=3600.*precip_xy[728:820,:,:]
precip_xy_d10[0:92,:,:]=3600.*precip_xy[819:911,:,:]

#weakly forcing simulation
#precip_xy_d2[0:92,:,:]=3600.*precip_xy[90:182,:,:]
#precip_xy_d3[0:92,:,:]=3600.*precip_xy[179:271,:,:]
#precip_xy_d4[0:92,:,:]=3600.*precip_xy[269:361,:,:]
#precip_xy_d5[0:92,:,:]=3600.*precip_xy[359:451,:,:]
#precip_xy_d6[0:92,:,:]=3600.*precip_xy[449:541,:,:]
#precip_xy_d7[0:92,:,:]=3600.*precip_xy[539:631,:,:]
#precip_xy_d8[0:92,:,:]=3600.*precip_xy[629:721,:,:]
#precip_xy_d9[0:92,:,:]=3600.*precip_xy[720:812,:,:]
#precip_xy_d10[0:92,:,:]=3600.*precip_xy[811:903,:,:]

#strongly forcing simulation
#precip_xy_d2[0:92,:,:]=3600.*precip_xy[91:183,:,:]
#precip_xy_d3[0:92,:,:]=3600.*precip_xy[183:275,:,:]
#precip_xy_d4[0:92,:,:]=3600.*precip_xy[275:367,:,:]
#precip_xy_d5[0:92,:,:]=3600.*precip_xy[367:459,:,:]
#precip_xy_d6[0:92,:,:]=3600.*precip_xy[459:551,:,:]
#precip_xy_d7[0:92,:,:]=3600.*precip_xy[550:642,:,:]
#precip_xy_d8[0:92,:,:]=3600.*precip_xy[642:734,:,:]
#precip_xy_d9[0:92,:,:]=3600.*precip_xy[733:825,:,:]
#precip_xy_d10[0:92,:,:]=3600.*precip_xy[825:917,:,:]

#rthqv: simulations following homogenization
#precip_xy_d2[0:92,:,:]=3600.*precip_xy_rthqvd2[88:180,:,:]
#precip_xy_d3[0:92,:,:]=3600.*precip_xy_rthqvd3[181:273,:,:]
#precip_xy_d4[0:92,:,:]=3600.*precip_xy_rthqvd4[273:365,:,:]
#precip_xy_d5[0:92,:,:]=3600.*precip_xy_rthqvd5[362:454,:,:]
#precip_xy_d6[0:92,:,:]=3600.*precip_xy_rthqvd6[453:545,:,:]
#precip_xy_d7[0:92,:,:]=3600.*precip_xy_rthqvd7[544:636,:,:]
#precip_xy_d8[0:92,:,:]=3600.*precip_xy_rthqvd8[638:730,:,:]
#precip_xy_d9[0:92,:,:]=3600.*precip_xy_rthqvd9[727:819,:,:]
#precip_xy_d10[0:92,:,:]=3600.*precip_xy_rthqvd10[819:911,:,:]


hour[0]=0.25
for ihr in N.arange(n_30min-1):
     hour[ihr+1]=hour[ihr]+0.25


#
#4*4KM
#
nx4=26
ny4=26
length_x4=N.zeros(nx4)
length_y4=N.zeros(ny4)
dx4=4000
dy4=4000


precip_44_ctrl_d2=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d3=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d4=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d5=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d6=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d7=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d8=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d9=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d10=N.zeros((n_30min,nx4,ny4))


precip_44_ctrl_d2_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d3_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d4_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d5_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d6_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d7_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d8_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d9_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d10_Day_mean=N.zeros((nx4,ny4))

surf_precip_44_contr_hist_Day_mean=N.zeros((81))
precip_44_ctrl_d6_hist_Day_mean=N.zeros((81))
precip_44_ctrl_d4_hist_Day_mean=N.zeros((81))
precip_44_ctrl_d5_hist_Day_mean=N.zeros((81))
precip_44_ctrl_d3_hist_Day_mean=N.zeros((81))


precip_44_ctrl_d2_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR_lag=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRNR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRNR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRNR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d2_inst_PRNR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRNR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRNR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRNR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRNR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRNR_lag1=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag1=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag2=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag3=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag4=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag6=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag5=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag7=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag8=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag9=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag10=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag11=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag12=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag13=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag14=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag15=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag16=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag17=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag18=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag19=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag20=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag21=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag22=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag23=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag24=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag25=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag26=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag27=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag28=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag29=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag30=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag31=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag32=N.zeros((ntr, n_30min))


precip_44_ctrl_d2_inst_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag33=N.zeros((ntr, n_30min))



precip_44_ctrl_d2_inst_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag34=N.zeros((ntr, n_30min))




precip_44_ctrl_d2_inst_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag35=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag36=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag37=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag38=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag39=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag40=N.zeros((ntr, n_30min))




precip_44_ctrl_d2_inst_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag41=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag42=N.zeros((ntr, n_30min))


precip_44_ctrl_d2_inst_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag43=N.zeros((ntr, n_30min))



precip_44_ctrl_d2_inst_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag44=N.zeros((ntr, n_30min))




precip_44_ctrl_d2_inst_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag45=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag46=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag47=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag48=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag49=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag50=N.zeros((ntr, n_30min))


precip_44_ctrl_d2_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))




precip_44_ctrl_d2_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))



precip_44_ctrl_d2_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))




precip_44_ctrl_d2_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))


precip_44_ctrl_d2_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag1=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag2=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag3=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag4=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag5=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag6=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag7=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag8=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag9=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag10=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag11=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag12=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag13=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag14=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag15=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag16=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag17=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag18=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag19=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag20=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag21=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag22=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag23=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag24=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag25=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag26=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag27=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag28=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag29=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag30=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag31=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag32=N.zeros((ntr, n_30min))


precip_44_ctrl_d2_inst_PRR_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag33=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag34=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag35=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag36=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag37=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag38=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag39=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag40=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag41=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag42=N.zeros((ntr, n_30min))




precip_44_ctrl_d2_inst_PRR_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag43=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag44=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag45=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag46=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag47=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag48=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag49=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag50=N.zeros((ntr, n_30min))

proRR_lag1_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag2_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag3_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag4_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag5_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag6_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag7_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag8_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag9_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag10_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag11_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))

length_x4[0]=0.#1000.
length_y4[0]=0.#1000.

for j in N.arange(ny4-1):
     length_x4[j+1]=length_x4[j]+dx4
     length_y4[j+1]=length_y4[j]+dy4

length_x4=length_x4/1000.
length_y4=length_y4/1000.

for j in N.arange(60):
     for ix in N.arange(nx4):
          for iy in N.arange(ny4):
#this is an example for an average over of 4X4 km square
               precip_44_ctrl_d2[j,ix,iy]=precip_xy_d2[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d3[j,ix,iy]=precip_xy_d3[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d4[j,ix,iy]=precip_xy_d4[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d5[j,ix,iy]=precip_xy_d5[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d6[j,ix,iy]=precip_xy_d6[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d7[j,ix,iy]=precip_xy_d7[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d8[j,ix,iy]=precip_xy_d8[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d9[j,ix,iy]=precip_xy_d9[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d10[j,ix,iy]=precip_xy_d10[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               if (ix == 25) :
                    precip_44_ctrl_d2[j,ix,iy]=precip_xy_d2[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d3[j,ix,iy]=precip_xy_d3[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d4[j,ix,iy]=precip_xy_d4[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d5[j,ix,iy]=precip_xy_d5[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d6[j,ix,iy]=precip_xy_d6[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d7[j,ix,iy]=precip_xy_d7[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d8[j,ix,iy]=precip_xy_d8[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d9[j,ix,iy]=precip_xy_d9[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d10[j,ix,iy]=precip_xy_d10[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
               if (iy == 25) :
                    precip_44_ctrl_d2[j,ix,iy]=precip_xy_d2[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d3[j,ix,iy]=precip_xy_d3[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d4[j,ix,iy]=precip_xy_d4[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d5[j,ix,iy]=precip_xy_d5[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d6[j,ix,iy]=precip_xy_d6[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d7[j,ix,iy]=precip_xy_d7[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d8[j,ix,iy]=precip_xy_d8[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d9[j,ix,iy]=precip_xy_d9[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d10[j,ix,iy]=precip_xy_d10[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()

#this is an example grid scale analysis
#               precip_44_ctrl_d2[j,ix,iy]=precip_xy_d2[j,ix,iy]
 #              precip_44_ctrl_d3[j,ix,iy]=precip_xy_d3[j,ix,iy]
  #             precip_44_ctrl_d4[j,ix,iy]=precip_xy_d4[j,ix,iy]
   #            precip_44_ctrl_d5[j,ix,iy]=precip_xy_d5[j,ix,iy]
    #           precip_44_ctrl_d6[j,ix,iy]=precip_xy_d6[j,ix,iy]
     #          precip_44_ctrl_d7[j,ix,iy]=precip_xy_d7[j,ix,iy]
      #         precip_44_ctrl_d8[j,ix,iy]=precip_xy_d8[j,ix,iy]
       #        precip_44_ctrl_d9[j,ix,iy]=precip_xy_d9[j,ix,iy]
        #       precip_44_ctrl_d10[j,ix,iy]=precip_xy_d10[j,ix,iy]



precip_44_mean=( precip_44_ctrl_d2[0:80,:,:].mean()+precip_44_ctrl_d3[0:80,:,:].mean()+precip_44_ctrl_d4[0:80,:,:].mean()+precip_44_ctrl_d5[0:80,:,:].mean()+precip_44_ctrl_d6[0:80,:,:].mean()+precip_44_ctrl_d7[0:80,:,:].mean()+precip_44_ctrl_d8[0:80,:,:].mean()+precip_44_ctrl_d9[0:80,:,:].mean()+precip_44_ctrl_d10[0:80,:,:].mean())/9.


Surf_precip_threshold_Day_mean[i]=precip_44_mean*0.5 #50% of the domain-mean daily-mean rain rate
     

for ix in N.arange(nx4):
     for iy in N.arange(ny4):
          for  j in N.arange(n_30min-1):
               for itr in N.arange(ntr):
                    if ( precip_44_ctrl_d2[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d2_inst_PR_lag[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d3[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d3_inst_PR_lag[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d4[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d4_inst_PR_lag[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d5[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d5_inst_PR_lag[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d6[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d6_inst_PR_lag[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d7[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d7_inst_PR_lag[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d8[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d8_inst_PR_lag[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d9[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d9_inst_PR_lag[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d10[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d10_inst_PR_lag[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j]+1.

                    if ( j > 0 ):
                         if ( precip_44_ctrl_d2[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d2[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag1[itr,j]=precip_44_ctrl_d2_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d2[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag2[itr,j]=precip_44_ctrl_d2_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d2[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag3[itr,j]=precip_44_ctrl_d2_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d2[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag4[itr,j]=precip_44_ctrl_d2_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d2[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag5[itr,j]=precip_44_ctrl_d2_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d2[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag6[itr,j]=precip_44_ctrl_d2_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d2[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag7[itr,j]=precip_44_ctrl_d2_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d2[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag8[itr,j]=precip_44_ctrl_d2_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d2[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag9[itr,j]=precip_44_ctrl_d2_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d2[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag10[itr,j]=precip_44_ctrl_d2_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d2[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag11[itr,j]=precip_44_ctrl_d2_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d2[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag12[itr,j]=precip_44_ctrl_d2_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d2[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag13[itr,j]=precip_44_ctrl_d2_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d2[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag14[itr,j]=precip_44_ctrl_d2_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d2[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag15[itr,j]=precip_44_ctrl_d2_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d2[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag16[itr,j]=precip_44_ctrl_d2_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d2[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag17[itr,j]=precip_44_ctrl_d2_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d2[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag18[itr,j]=precip_44_ctrl_d2_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d2[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag19[itr,j]=precip_44_ctrl_d2_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d2[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag20[itr,j]=precip_44_ctrl_d2_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d2[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag21[itr,j]=precip_44_ctrl_d2_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d2[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag22[itr,j]=precip_44_ctrl_d2_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d2[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag23[itr,j]=precip_44_ctrl_d2_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d2[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag24[itr,j]=precip_44_ctrl_d2_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d2[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag25[itr,j]=precip_44_ctrl_d2_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d2[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag26[itr,j]=precip_44_ctrl_d2_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d2[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag27[itr,j]=precip_44_ctrl_d2_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d2[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag28[itr,j]=precip_44_ctrl_d2_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d2[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag29[itr,j]=precip_44_ctrl_d2_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d2[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag30[itr,j]=precip_44_ctrl_d2_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d2[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag31[itr,j]=precip_44_ctrl_d2_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d2[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag32[itr,j]=precip_44_ctrl_d2_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d3[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d3[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag1[itr,j]=precip_44_ctrl_d3_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d3[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag2[itr,j]=precip_44_ctrl_d3_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d3[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag3[itr,j]=precip_44_ctrl_d3_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d3[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag4[itr,j]=precip_44_ctrl_d3_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d3[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag5[itr,j]=precip_44_ctrl_d3_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d3[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag6[itr,j]=precip_44_ctrl_d3_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d3[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag7[itr,j]=precip_44_ctrl_d3_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d3[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag8[itr,j]=precip_44_ctrl_d3_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d3[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag9[itr,j]=precip_44_ctrl_d3_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d3[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag10[itr,j]=precip_44_ctrl_d3_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d3[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag11[itr,j]=precip_44_ctrl_d3_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d3[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag12[itr,j]=precip_44_ctrl_d3_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d3[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag13[itr,j]=precip_44_ctrl_d3_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d3[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag14[itr,j]=precip_44_ctrl_d3_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d3[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag15[itr,j]=precip_44_ctrl_d3_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d3[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag16[itr,j]=precip_44_ctrl_d3_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d3[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag17[itr,j]=precip_44_ctrl_d3_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d3[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag18[itr,j]=precip_44_ctrl_d3_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d3[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag19[itr,j]=precip_44_ctrl_d3_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d3[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag20[itr,j]=precip_44_ctrl_d3_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d3[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag21[itr,j]=precip_44_ctrl_d3_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d3[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag22[itr,j]=precip_44_ctrl_d3_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d3[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag23[itr,j]=precip_44_ctrl_d3_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d3[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag24[itr,j]=precip_44_ctrl_d3_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d3[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag25[itr,j]=precip_44_ctrl_d3_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d3[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag26[itr,j]=precip_44_ctrl_d3_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d3[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag27[itr,j]=precip_44_ctrl_d3_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d3[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag28[itr,j]=precip_44_ctrl_d3_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d3[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag29[itr,j]=precip_44_ctrl_d3_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d3[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag30[itr,j]=precip_44_ctrl_d3_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d3[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag31[itr,j]=precip_44_ctrl_d3_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d3[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag32[itr,j]=precip_44_ctrl_d3_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d4[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d4[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag1[itr,j]=precip_44_ctrl_d4_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d4[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag2[itr,j]=precip_44_ctrl_d4_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d4[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag3[itr,j]=precip_44_ctrl_d4_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d4[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag4[itr,j]=precip_44_ctrl_d4_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d4[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag5[itr,j]=precip_44_ctrl_d4_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d4[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag6[itr,j]=precip_44_ctrl_d4_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d4[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag7[itr,j]=precip_44_ctrl_d4_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d4[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag8[itr,j]=precip_44_ctrl_d4_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d4[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag9[itr,j]=precip_44_ctrl_d4_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d4[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag10[itr,j]=precip_44_ctrl_d4_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d4[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag11[itr,j]=precip_44_ctrl_d4_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d4[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag12[itr,j]=precip_44_ctrl_d4_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d4[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag13[itr,j]=precip_44_ctrl_d4_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d4[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag14[itr,j]=precip_44_ctrl_d4_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d4[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag15[itr,j]=precip_44_ctrl_d4_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d4[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag16[itr,j]=precip_44_ctrl_d4_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d4[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag17[itr,j]=precip_44_ctrl_d4_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d4[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag18[itr,j]=precip_44_ctrl_d4_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d4[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag19[itr,j]=precip_44_ctrl_d4_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d4[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag20[itr,j]=precip_44_ctrl_d4_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d4[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag21[itr,j]=precip_44_ctrl_d4_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d4[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag22[itr,j]=precip_44_ctrl_d4_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d4[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag23[itr,j]=precip_44_ctrl_d4_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d4[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag24[itr,j]=precip_44_ctrl_d4_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d4[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag25[itr,j]=precip_44_ctrl_d4_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d4[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag26[itr,j]=precip_44_ctrl_d4_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d4[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag27[itr,j]=precip_44_ctrl_d4_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d4[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag28[itr,j]=precip_44_ctrl_d4_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d4[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag29[itr,j]=precip_44_ctrl_d4_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d4[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag30[itr,j]=precip_44_ctrl_d4_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d4[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag31[itr,j]=precip_44_ctrl_d4_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d4[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag32[itr,j]=precip_44_ctrl_d4_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d5[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d5[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag1[itr,j]=precip_44_ctrl_d5_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d5[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag2[itr,j]=precip_44_ctrl_d5_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d5[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag3[itr,j]=precip_44_ctrl_d5_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d5[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag4[itr,j]=precip_44_ctrl_d5_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d5[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag5[itr,j]=precip_44_ctrl_d5_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d5[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag6[itr,j]=precip_44_ctrl_d5_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d5[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag7[itr,j]=precip_44_ctrl_d5_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d5[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag8[itr,j]=precip_44_ctrl_d5_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d5[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag9[itr,j]=precip_44_ctrl_d5_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d5[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag10[itr,j]=precip_44_ctrl_d5_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d5[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag11[itr,j]=precip_44_ctrl_d5_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d5[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag12[itr,j]=precip_44_ctrl_d5_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d5[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag13[itr,j]=precip_44_ctrl_d5_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d5[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag14[itr,j]=precip_44_ctrl_d5_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d5[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag15[itr,j]=precip_44_ctrl_d5_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d5[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag16[itr,j]=precip_44_ctrl_d5_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d5[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag17[itr,j]=precip_44_ctrl_d5_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d5[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag18[itr,j]=precip_44_ctrl_d5_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d5[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag19[itr,j]=precip_44_ctrl_d5_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d5[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag20[itr,j]=precip_44_ctrl_d5_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d5[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag21[itr,j]=precip_44_ctrl_d5_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d5[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag22[itr,j]=precip_44_ctrl_d5_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d5[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag23[itr,j]=precip_44_ctrl_d5_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d5[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag24[itr,j]=precip_44_ctrl_d5_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d5[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag25[itr,j]=precip_44_ctrl_d5_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d5[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag26[itr,j]=precip_44_ctrl_d5_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d5[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag27[itr,j]=precip_44_ctrl_d5_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d5[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag28[itr,j]=precip_44_ctrl_d5_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d5[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag29[itr,j]=precip_44_ctrl_d5_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d5[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag30[itr,j]=precip_44_ctrl_d5_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d5[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag31[itr,j]=precip_44_ctrl_d5_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d5[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag32[itr,j]=precip_44_ctrl_d5_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d6[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d6[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag1[itr,j]=precip_44_ctrl_d6_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d6[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag2[itr,j]=precip_44_ctrl_d6_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d6[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag3[itr,j]=precip_44_ctrl_d6_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d6[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag4[itr,j]=precip_44_ctrl_d6_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d6[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag5[itr,j]=precip_44_ctrl_d6_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d6[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag6[itr,j]=precip_44_ctrl_d6_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d6[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag7[itr,j]=precip_44_ctrl_d6_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d6[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag8[itr,j]=precip_44_ctrl_d6_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d6[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag9[itr,j]=precip_44_ctrl_d6_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d6[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag10[itr,j]=precip_44_ctrl_d6_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d6[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag11[itr,j]=precip_44_ctrl_d6_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d6[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag12[itr,j]=precip_44_ctrl_d6_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d6[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag13[itr,j]=precip_44_ctrl_d6_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d6[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag14[itr,j]=precip_44_ctrl_d6_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d6[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag15[itr,j]=precip_44_ctrl_d6_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d6[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag16[itr,j]=precip_44_ctrl_d6_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d6[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag17[itr,j]=precip_44_ctrl_d6_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d6[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag18[itr,j]=precip_44_ctrl_d6_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d6[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag19[itr,j]=precip_44_ctrl_d6_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d6[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag20[itr,j]=precip_44_ctrl_d6_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d6[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag21[itr,j]=precip_44_ctrl_d6_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d6[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag22[itr,j]=precip_44_ctrl_d6_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d6[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag23[itr,j]=precip_44_ctrl_d6_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d6[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag24[itr,j]=precip_44_ctrl_d6_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d6[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag25[itr,j]=precip_44_ctrl_d6_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d6[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag26[itr,j]=precip_44_ctrl_d6_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d6[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag27[itr,j]=precip_44_ctrl_d6_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d6[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag28[itr,j]=precip_44_ctrl_d6_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d6[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag29[itr,j]=precip_44_ctrl_d6_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d6[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag30[itr,j]=precip_44_ctrl_d6_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d6[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag31[itr,j]=precip_44_ctrl_d6_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d6[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag32[itr,j]=precip_44_ctrl_d6_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d7[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d7[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag1[itr,j]=precip_44_ctrl_d7_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d7[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag2[itr,j]=precip_44_ctrl_d7_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d7[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag3[itr,j]=precip_44_ctrl_d7_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d7[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag4[itr,j]=precip_44_ctrl_d7_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d7[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag5[itr,j]=precip_44_ctrl_d7_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d7[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag6[itr,j]=precip_44_ctrl_d7_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d7[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag7[itr,j]=precip_44_ctrl_d7_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d7[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag8[itr,j]=precip_44_ctrl_d7_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d7[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag9[itr,j]=precip_44_ctrl_d7_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d7[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag10[itr,j]=precip_44_ctrl_d7_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d7[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag11[itr,j]=precip_44_ctrl_d7_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d7[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag12[itr,j]=precip_44_ctrl_d7_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d7[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag13[itr,j]=precip_44_ctrl_d7_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d7[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag14[itr,j]=precip_44_ctrl_d7_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d7[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag15[itr,j]=precip_44_ctrl_d7_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d7[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag16[itr,j]=precip_44_ctrl_d7_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d7[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag17[itr,j]=precip_44_ctrl_d7_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d7[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag18[itr,j]=precip_44_ctrl_d7_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d7[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag19[itr,j]=precip_44_ctrl_d7_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d7[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag20[itr,j]=precip_44_ctrl_d7_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d7[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag21[itr,j]=precip_44_ctrl_d7_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d7[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag22[itr,j]=precip_44_ctrl_d7_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d7[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag23[itr,j]=precip_44_ctrl_d7_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d7[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag24[itr,j]=precip_44_ctrl_d7_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d7[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag25[itr,j]=precip_44_ctrl_d7_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d7[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag26[itr,j]=precip_44_ctrl_d7_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d7[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag27[itr,j]=precip_44_ctrl_d7_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d7[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag28[itr,j]=precip_44_ctrl_d7_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d7[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag29[itr,j]=precip_44_ctrl_d7_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d7[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag30[itr,j]=precip_44_ctrl_d7_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d7[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag31[itr,j]=precip_44_ctrl_d7_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d7[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag32[itr,j]=precip_44_ctrl_d7_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d8[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d8[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag1[itr,j]=precip_44_ctrl_d8_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d8[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag2[itr,j]=precip_44_ctrl_d8_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d8[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag3[itr,j]=precip_44_ctrl_d8_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d8[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag4[itr,j]=precip_44_ctrl_d8_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d8[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag5[itr,j]=precip_44_ctrl_d8_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d8[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag6[itr,j]=precip_44_ctrl_d8_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d8[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag7[itr,j]=precip_44_ctrl_d8_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d8[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag8[itr,j]=precip_44_ctrl_d8_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d8[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag9[itr,j]=precip_44_ctrl_d8_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d8[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag10[itr,j]=precip_44_ctrl_d8_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d8[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag11[itr,j]=precip_44_ctrl_d8_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d8[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag12[itr,j]=precip_44_ctrl_d8_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d8[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag13[itr,j]=precip_44_ctrl_d8_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d8[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag14[itr,j]=precip_44_ctrl_d8_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d8[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag15[itr,j]=precip_44_ctrl_d8_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d8[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag16[itr,j]=precip_44_ctrl_d8_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d8[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag17[itr,j]=precip_44_ctrl_d8_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d8[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag18[itr,j]=precip_44_ctrl_d8_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d8[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag19[itr,j]=precip_44_ctrl_d8_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d8[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag20[itr,j]=precip_44_ctrl_d8_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d8[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag21[itr,j]=precip_44_ctrl_d8_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d8[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag22[itr,j]=precip_44_ctrl_d8_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d8[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag23[itr,j]=precip_44_ctrl_d8_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d8[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag24[itr,j]=precip_44_ctrl_d8_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d8[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag25[itr,j]=precip_44_ctrl_d8_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d8[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag26[itr,j]=precip_44_ctrl_d8_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d8[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag27[itr,j]=precip_44_ctrl_d8_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d8[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag28[itr,j]=precip_44_ctrl_d8_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d8[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag29[itr,j]=precip_44_ctrl_d8_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d8[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag30[itr,j]=precip_44_ctrl_d8_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d8[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag31[itr,j]=precip_44_ctrl_d8_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d8[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag32[itr,j]=precip_44_ctrl_d8_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d9[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d9[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag1[itr,j]=precip_44_ctrl_d9_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d9[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag2[itr,j]=precip_44_ctrl_d9_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d9[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag3[itr,j]=precip_44_ctrl_d9_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d9[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag4[itr,j]=precip_44_ctrl_d9_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d9[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag5[itr,j]=precip_44_ctrl_d9_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d9[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag6[itr,j]=precip_44_ctrl_d9_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d9[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag7[itr,j]=precip_44_ctrl_d9_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d9[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag8[itr,j]=precip_44_ctrl_d9_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d9[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag9[itr,j]=precip_44_ctrl_d9_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d9[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag10[itr,j]=precip_44_ctrl_d9_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d9[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag11[itr,j]=precip_44_ctrl_d9_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d9[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag12[itr,j]=precip_44_ctrl_d9_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d9[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag13[itr,j]=precip_44_ctrl_d9_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d9[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag14[itr,j]=precip_44_ctrl_d9_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d9[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag15[itr,j]=precip_44_ctrl_d9_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d9[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag16[itr,j]=precip_44_ctrl_d9_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d9[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag17[itr,j]=precip_44_ctrl_d9_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d9[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag18[itr,j]=precip_44_ctrl_d9_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d9[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag19[itr,j]=precip_44_ctrl_d9_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d9[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag20[itr,j]=precip_44_ctrl_d9_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d9[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag21[itr,j]=precip_44_ctrl_d9_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d9[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag22[itr,j]=precip_44_ctrl_d9_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d9[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag23[itr,j]=precip_44_ctrl_d9_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d9[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag24[itr,j]=precip_44_ctrl_d9_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d9[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag25[itr,j]=precip_44_ctrl_d9_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d9[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag26[itr,j]=precip_44_ctrl_d9_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d9[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag27[itr,j]=precip_44_ctrl_d9_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d9[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag28[itr,j]=precip_44_ctrl_d9_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d9[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag29[itr,j]=precip_44_ctrl_d9_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d9[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag30[itr,j]=precip_44_ctrl_d9_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d9[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag31[itr,j]=precip_44_ctrl_d9_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d9[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag32[itr,j]=precip_44_ctrl_d9_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d10[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d10[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag1[itr,j]=precip_44_ctrl_d10_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d10[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag2[itr,j]=precip_44_ctrl_d10_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d10[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag3[itr,j]=precip_44_ctrl_d10_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d10[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag4[itr,j]=precip_44_ctrl_d10_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d10[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag5[itr,j]=precip_44_ctrl_d10_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d10[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag6[itr,j]=precip_44_ctrl_d10_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d10[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag7[itr,j]=precip_44_ctrl_d10_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d10[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag8[itr,j]=precip_44_ctrl_d10_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d10[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag9[itr,j]=precip_44_ctrl_d10_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d10[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag10[itr,j]=precip_44_ctrl_d10_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d10[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag11[itr,j]=precip_44_ctrl_d10_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d10[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag12[itr,j]=precip_44_ctrl_d10_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d10[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag13[itr,j]=precip_44_ctrl_d10_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d10[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag14[itr,j]=precip_44_ctrl_d10_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d10[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag15[itr,j]=precip_44_ctrl_d10_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d10[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag16[itr,j]=precip_44_ctrl_d10_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d10[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag17[itr,j]=precip_44_ctrl_d10_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d10[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag18[itr,j]=precip_44_ctrl_d10_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d10[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag19[itr,j]=precip_44_ctrl_d10_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d10[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag20[itr,j]=precip_44_ctrl_d10_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d10[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag21[itr,j]=precip_44_ctrl_d10_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d10[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag22[itr,j]=precip_44_ctrl_d10_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d10[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag23[itr,j]=precip_44_ctrl_d10_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d10[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag24[itr,j]=precip_44_ctrl_d10_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d10[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag25[itr,j]=precip_44_ctrl_d10_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d10[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag26[itr,j]=precip_44_ctrl_d10_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d10[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag27[itr,j]=precip_44_ctrl_d10_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d10[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag28[itr,j]=precip_44_ctrl_d10_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d10[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag29[itr,j]=precip_44_ctrl_d10_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d10[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag30[itr,j]=precip_44_ctrl_d10_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d10[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag31[itr,j]=precip_44_ctrl_d10_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d10[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag32[itr,j]=precip_44_ctrl_d10_inst_PRR_lag32[itr,j]+1.




precip_44_ctrl_d2_inst_PR_lag[:,:]=precip_44_ctrl_d2_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag1[:,:]=precip_44_ctrl_d2_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag2[:,:]=precip_44_ctrl_d2_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag3[:,:]=precip_44_ctrl_d2_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag4[:,:]=precip_44_ctrl_d2_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag5[:,:]=precip_44_ctrl_d2_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag6[:,:]=precip_44_ctrl_d2_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag7[:,:]=precip_44_ctrl_d2_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag8[:,:]=precip_44_ctrl_d2_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag9[:,:]=precip_44_ctrl_d2_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag10[:,:]=precip_44_ctrl_d2_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag11[:,:]=precip_44_ctrl_d2_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag12[:,:]=precip_44_ctrl_d2_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag13[:,:]=precip_44_ctrl_d2_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag14[:,:]=precip_44_ctrl_d2_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag15[:,:]=precip_44_ctrl_d2_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag16[:,:]=precip_44_ctrl_d2_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag17[:,:]=precip_44_ctrl_d2_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag18[:,:]=precip_44_ctrl_d2_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag19[:,:]=precip_44_ctrl_d2_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag20[:,:]=precip_44_ctrl_d2_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag21[:,:]=precip_44_ctrl_d2_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag22[:,:]=precip_44_ctrl_d2_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag23[:,:]=precip_44_ctrl_d2_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag24[:,:]=precip_44_ctrl_d2_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag25[:,:]=precip_44_ctrl_d2_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag26[:,:]=precip_44_ctrl_d2_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag27[:,:]=precip_44_ctrl_d2_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag28[:,:]=precip_44_ctrl_d2_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag29[:,:]=precip_44_ctrl_d2_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag30[:,:]=precip_44_ctrl_d2_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag31[:,:]=precip_44_ctrl_d2_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag32[:,:]=precip_44_ctrl_d2_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)

precip_44_ctrl_d3_inst_PR_lag[:,:]=precip_44_ctrl_d3_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag1[:,:]=precip_44_ctrl_d3_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag2[:,:]=precip_44_ctrl_d3_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag3[:,:]=precip_44_ctrl_d3_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag4[:,:]=precip_44_ctrl_d3_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag5[:,:]=precip_44_ctrl_d3_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag6[:,:]=precip_44_ctrl_d3_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag7[:,:]=precip_44_ctrl_d3_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag8[:,:]=precip_44_ctrl_d3_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag9[:,:]=precip_44_ctrl_d3_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag10[:,:]=precip_44_ctrl_d3_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag11[:,:]=precip_44_ctrl_d3_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag12[:,:]=precip_44_ctrl_d3_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag13[:,:]=precip_44_ctrl_d3_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag14[:,:]=precip_44_ctrl_d3_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag15[:,:]=precip_44_ctrl_d3_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag16[:,:]=precip_44_ctrl_d3_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag17[:,:]=precip_44_ctrl_d3_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag18[:,:]=precip_44_ctrl_d3_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag19[:,:]=precip_44_ctrl_d3_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag20[:,:]=precip_44_ctrl_d3_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag21[:,:]=precip_44_ctrl_d3_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag22[:,:]=precip_44_ctrl_d3_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag23[:,:]=precip_44_ctrl_d3_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag24[:,:]=precip_44_ctrl_d3_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag25[:,:]=precip_44_ctrl_d3_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag26[:,:]=precip_44_ctrl_d3_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag27[:,:]=precip_44_ctrl_d3_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag28[:,:]=precip_44_ctrl_d3_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag29[:,:]=precip_44_ctrl_d3_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag30[:,:]=precip_44_ctrl_d3_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag31[:,:]=precip_44_ctrl_d3_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag32[:,:]=precip_44_ctrl_d3_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)


precip_44_ctrl_d4_inst_PR_lag[:,:]=precip_44_ctrl_d4_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag1[:,:]=precip_44_ctrl_d4_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRNR_lag1[:,:]=precip_44_ctrl_d4_inst_PRNR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag2[:,:]=precip_44_ctrl_d4_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag3[:,:]=precip_44_ctrl_d4_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag4[:,:]=precip_44_ctrl_d4_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag5[:,:]=precip_44_ctrl_d4_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag6[:,:]=precip_44_ctrl_d4_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag7[:,:]=precip_44_ctrl_d4_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag8[:,:]=precip_44_ctrl_d4_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag9[:,:]=precip_44_ctrl_d4_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag10[:,:]=precip_44_ctrl_d4_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag11[:,:]=precip_44_ctrl_d4_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag12[:,:]=precip_44_ctrl_d4_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag13[:,:]=precip_44_ctrl_d4_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag14[:,:]=precip_44_ctrl_d4_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag15[:,:]=precip_44_ctrl_d4_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag16[:,:]=precip_44_ctrl_d4_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag17[:,:]=precip_44_ctrl_d4_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag18[:,:]=precip_44_ctrl_d4_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag19[:,:]=precip_44_ctrl_d4_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag20[:,:]=precip_44_ctrl_d4_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag21[:,:]=precip_44_ctrl_d4_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag22[:,:]=precip_44_ctrl_d4_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag23[:,:]=precip_44_ctrl_d4_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag24[:,:]=precip_44_ctrl_d4_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag25[:,:]=precip_44_ctrl_d4_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag26[:,:]=precip_44_ctrl_d4_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag27[:,:]=precip_44_ctrl_d4_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag28[:,:]=precip_44_ctrl_d4_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag29[:,:]=precip_44_ctrl_d4_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag30[:,:]=precip_44_ctrl_d4_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag31[:,:]=precip_44_ctrl_d4_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag32[:,:]=precip_44_ctrl_d4_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)

precip_44_ctrl_d5_inst_PR_lag[:,:]=precip_44_ctrl_d5_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag1[:,:]=precip_44_ctrl_d5_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRNR_lag1[:,:]=precip_44_ctrl_d5_inst_PRNR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag2[:,:]=precip_44_ctrl_d5_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag3[:,:]=precip_44_ctrl_d5_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag4[:,:]=precip_44_ctrl_d5_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag5[:,:]=precip_44_ctrl_d5_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag6[:,:]=precip_44_ctrl_d5_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag7[:,:]=precip_44_ctrl_d5_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag8[:,:]=precip_44_ctrl_d5_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag9[:,:]=precip_44_ctrl_d5_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag10[:,:]=precip_44_ctrl_d5_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag11[:,:]=precip_44_ctrl_d5_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag12[:,:]=precip_44_ctrl_d5_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag13[:,:]=precip_44_ctrl_d5_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag14[:,:]=precip_44_ctrl_d5_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag15[:,:]=precip_44_ctrl_d5_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag16[:,:]=precip_44_ctrl_d5_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag17[:,:]=precip_44_ctrl_d5_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag18[:,:]=precip_44_ctrl_d5_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag19[:,:]=precip_44_ctrl_d5_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag20[:,:]=precip_44_ctrl_d5_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag21[:,:]=precip_44_ctrl_d5_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag22[:,:]=precip_44_ctrl_d5_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag23[:,:]=precip_44_ctrl_d5_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag24[:,:]=precip_44_ctrl_d5_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag25[:,:]=precip_44_ctrl_d5_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag26[:,:]=precip_44_ctrl_d5_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag27[:,:]=precip_44_ctrl_d5_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag28[:,:]=precip_44_ctrl_d5_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag29[:,:]=precip_44_ctrl_d5_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag30[:,:]=precip_44_ctrl_d5_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag31[:,:]=precip_44_ctrl_d5_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag32[:,:]=precip_44_ctrl_d5_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)

precip_44_ctrl_d6_inst_PR_lag[:,:]=precip_44_ctrl_d6_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag1[:,:]=precip_44_ctrl_d6_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRNR_lag1[:,:]=precip_44_ctrl_d6_inst_PRNR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag2[:,:]=precip_44_ctrl_d6_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag3[:,:]=precip_44_ctrl_d6_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag4[:,:]=precip_44_ctrl_d6_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag5[:,:]=precip_44_ctrl_d6_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag6[:,:]=precip_44_ctrl_d6_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag7[:,:]=precip_44_ctrl_d6_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag8[:,:]=precip_44_ctrl_d6_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag9[:,:]=precip_44_ctrl_d6_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag10[:,:]=precip_44_ctrl_d6_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag11[:,:]=precip_44_ctrl_d6_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag12[:,:]=precip_44_ctrl_d6_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag13[:,:]=precip_44_ctrl_d6_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag14[:,:]=precip_44_ctrl_d6_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag15[:,:]=precip_44_ctrl_d6_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag16[:,:]=precip_44_ctrl_d6_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag17[:,:]=precip_44_ctrl_d6_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag18[:,:]=precip_44_ctrl_d6_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag19[:,:]=precip_44_ctrl_d6_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag20[:,:]=precip_44_ctrl_d6_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag21[:,:]=precip_44_ctrl_d6_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag22[:,:]=precip_44_ctrl_d6_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag23[:,:]=precip_44_ctrl_d6_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag24[:,:]=precip_44_ctrl_d6_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag25[:,:]=precip_44_ctrl_d6_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag26[:,:]=precip_44_ctrl_d6_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag27[:,:]=precip_44_ctrl_d6_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag28[:,:]=precip_44_ctrl_d6_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag29[:,:]=precip_44_ctrl_d6_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag30[:,:]=precip_44_ctrl_d6_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag31[:,:]=precip_44_ctrl_d6_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag32[:,:]=precip_44_ctrl_d6_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)

precip_44_ctrl_d7_inst_PR_lag[:,:]=precip_44_ctrl_d7_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag1[:,:]=precip_44_ctrl_d7_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRNR_lag1[:,:]=precip_44_ctrl_d7_inst_PRNR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag2[:,:]=precip_44_ctrl_d7_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag3[:,:]=precip_44_ctrl_d7_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag4[:,:]=precip_44_ctrl_d7_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag5[:,:]=precip_44_ctrl_d7_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag6[:,:]=precip_44_ctrl_d7_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag7[:,:]=precip_44_ctrl_d7_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag8[:,:]=precip_44_ctrl_d7_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag9[:,:]=precip_44_ctrl_d7_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag10[:,:]=precip_44_ctrl_d7_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag11[:,:]=precip_44_ctrl_d7_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag12[:,:]=precip_44_ctrl_d7_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag13[:,:]=precip_44_ctrl_d7_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag14[:,:]=precip_44_ctrl_d7_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag15[:,:]=precip_44_ctrl_d7_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag16[:,:]=precip_44_ctrl_d7_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag17[:,:]=precip_44_ctrl_d7_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag18[:,:]=precip_44_ctrl_d7_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag19[:,:]=precip_44_ctrl_d7_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag20[:,:]=precip_44_ctrl_d7_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag21[:,:]=precip_44_ctrl_d7_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag22[:,:]=precip_44_ctrl_d7_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag23[:,:]=precip_44_ctrl_d7_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag24[:,:]=precip_44_ctrl_d7_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag25[:,:]=precip_44_ctrl_d7_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag26[:,:]=precip_44_ctrl_d7_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag27[:,:]=precip_44_ctrl_d7_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag28[:,:]=precip_44_ctrl_d7_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag29[:,:]=precip_44_ctrl_d7_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag30[:,:]=precip_44_ctrl_d7_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag31[:,:]=precip_44_ctrl_d7_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag32[:,:]=precip_44_ctrl_d7_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)

precip_44_ctrl_d8_inst_PR_lag[:,:]=precip_44_ctrl_d8_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag1[:,:]=precip_44_ctrl_d8_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRNR_lag1[:,:]=precip_44_ctrl_d8_inst_PRNR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag2[:,:]=precip_44_ctrl_d8_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag3[:,:]=precip_44_ctrl_d8_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag4[:,:]=precip_44_ctrl_d8_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag5[:,:]=precip_44_ctrl_d8_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag6[:,:]=precip_44_ctrl_d8_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag7[:,:]=precip_44_ctrl_d8_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag8[:,:]=precip_44_ctrl_d8_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag9[:,:]=precip_44_ctrl_d8_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag10[:,:]=precip_44_ctrl_d8_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag11[:,:]=precip_44_ctrl_d8_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag12[:,:]=precip_44_ctrl_d8_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag13[:,:]=precip_44_ctrl_d8_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag14[:,:]=precip_44_ctrl_d8_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag15[:,:]=precip_44_ctrl_d8_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag16[:,:]=precip_44_ctrl_d8_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag17[:,:]=precip_44_ctrl_d8_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag18[:,:]=precip_44_ctrl_d8_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag19[:,:]=precip_44_ctrl_d8_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag20[:,:]=precip_44_ctrl_d8_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag21[:,:]=precip_44_ctrl_d8_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag22[:,:]=precip_44_ctrl_d8_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag23[:,:]=precip_44_ctrl_d8_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag24[:,:]=precip_44_ctrl_d8_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag25[:,:]=precip_44_ctrl_d8_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag26[:,:]=precip_44_ctrl_d8_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag27[:,:]=precip_44_ctrl_d8_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag28[:,:]=precip_44_ctrl_d8_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag29[:,:]=precip_44_ctrl_d8_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag30[:,:]=precip_44_ctrl_d8_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag31[:,:]=precip_44_ctrl_d8_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag32[:,:]=precip_44_ctrl_d8_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)

precip_44_ctrl_d9_inst_PR_lag[:,:]=precip_44_ctrl_d9_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag1[:,:]=precip_44_ctrl_d9_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRNR_lag1[:,:]=precip_44_ctrl_d9_inst_PRNR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag2[:,:]=precip_44_ctrl_d9_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag3[:,:]=precip_44_ctrl_d9_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag4[:,:]=precip_44_ctrl_d9_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag5[:,:]=precip_44_ctrl_d9_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag6[:,:]=precip_44_ctrl_d9_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag7[:,:]=precip_44_ctrl_d9_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag8[:,:]=precip_44_ctrl_d9_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag9[:,:]=precip_44_ctrl_d9_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag10[:,:]=precip_44_ctrl_d9_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag11[:,:]=precip_44_ctrl_d9_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag12[:,:]=precip_44_ctrl_d9_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag13[:,:]=precip_44_ctrl_d9_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag14[:,:]=precip_44_ctrl_d9_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag15[:,:]=precip_44_ctrl_d9_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag16[:,:]=precip_44_ctrl_d9_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag17[:,:]=precip_44_ctrl_d9_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag18[:,:]=precip_44_ctrl_d9_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag19[:,:]=precip_44_ctrl_d9_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag20[:,:]=precip_44_ctrl_d9_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag21[:,:]=precip_44_ctrl_d9_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag22[:,:]=precip_44_ctrl_d9_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag23[:,:]=precip_44_ctrl_d9_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag24[:,:]=precip_44_ctrl_d9_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag25[:,:]=precip_44_ctrl_d9_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag26[:,:]=precip_44_ctrl_d9_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag27[:,:]=precip_44_ctrl_d9_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag28[:,:]=precip_44_ctrl_d9_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag29[:,:]=precip_44_ctrl_d9_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag30[:,:]=precip_44_ctrl_d9_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag31[:,:]=precip_44_ctrl_d9_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag32[:,:]=precip_44_ctrl_d9_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)

precip_44_ctrl_d10_inst_PR_lag[:,:]=precip_44_ctrl_d10_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag1[:,:]=precip_44_ctrl_d10_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRNR_lag1[:,:]=precip_44_ctrl_d10_inst_PRNR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag2[:,:]=precip_44_ctrl_d10_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag3[:,:]=precip_44_ctrl_d10_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag4[:,:]=precip_44_ctrl_d10_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag5[:,:]=precip_44_ctrl_d10_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag6[:,:]=precip_44_ctrl_d10_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag7[:,:]=precip_44_ctrl_d10_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag8[:,:]=precip_44_ctrl_d10_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag9[:,:]=precip_44_ctrl_d10_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag10[:,:]=precip_44_ctrl_d10_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag11[:,:]=precip_44_ctrl_d10_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag12[:,:]=precip_44_ctrl_d10_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag13[:,:]=precip_44_ctrl_d10_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag14[:,:]=precip_44_ctrl_d10_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag15[:,:]=precip_44_ctrl_d10_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag16[:,:]=precip_44_ctrl_d10_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag17[:,:]=precip_44_ctrl_d10_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag18[:,:]=precip_44_ctrl_d10_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag19[:,:]=precip_44_ctrl_d10_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag20[:,:]=precip_44_ctrl_d10_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag21[:,:]=precip_44_ctrl_d10_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag22[:,:]=precip_44_ctrl_d10_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag23[:,:]=precip_44_ctrl_d10_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag24[:,:]=precip_44_ctrl_d10_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag25[:,:]=precip_44_ctrl_d10_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag26[:,:]=precip_44_ctrl_d10_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag27[:,:]=precip_44_ctrl_d10_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag28[:,:]=precip_44_ctrl_d10_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag29[:,:]=precip_44_ctrl_d10_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag30[:,:]=precip_44_ctrl_d10_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag31[:,:]=precip_44_ctrl_d10_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag32[:,:]=precip_44_ctrl_d10_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)


proR_lag_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag1_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag2_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag3_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag4_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag5_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag6_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag7_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag8_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag9_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag10_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag11_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag12_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag13_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag14_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag15_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag16_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag17_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag18_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag19_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag20_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag21_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag22_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag23_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag24_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag25_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag26_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag27_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag28_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag29_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag30_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag31_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proR2_lag32_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))

proRR_lag1_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag2_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag3_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag4_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag5_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag6_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag7_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag8_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag9_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag10_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag11_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag12_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag13_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag14_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag15_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag16_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag17_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag18_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag19_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag20_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag21_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag22_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag23_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag24_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag25_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag26_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag27_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag28_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag29_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag30_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag31_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_lag32_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))

proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag33_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag34_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag35_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag36_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag37_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag38_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag39_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag40_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag41_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag42_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag43_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag44_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag45_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag46_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag47_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag48_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag49_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))
proRR_proR2_lag50_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))

proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))


nlag=55
lag=N.zeros((nlag))
zero_lag=N.zeros((nlag))
lag[0]=0
for j in N.arange(nlag-1):
     lag[j+1]=lag[j]-0.25
zero_lag[:]=0.0

proR_lag_ens_mean_inst_44_vs_AveRR_thr=N.zeros((ntr,n_30min))

for itr in N.arange(ntr):
     for j in N.arange(n_30min-1):


          if ( j >=1 ) :
               precip_44_ctrl_d2_inst_PR2_lag1[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-1]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag1[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-1]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag1[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-1]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag1[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-1]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag1[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-1]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag1[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-1]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag1[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-1]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag1[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-1]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag1[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-1]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d2_inst_PRR_lag1[itr,j]-precip_44_ctrl_d2_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d3_inst_PRR_lag1[itr,j]-precip_44_ctrl_d3_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d4_inst_PRR_lag1[itr,j]-precip_44_ctrl_d4_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d5_inst_PRR_lag1[itr,j]-precip_44_ctrl_d5_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d6_inst_PRR_lag1[itr,j]-precip_44_ctrl_d6_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d7_inst_PRR_lag1[itr,j]-precip_44_ctrl_d7_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d8_inst_PRR_lag1[itr,j]-precip_44_ctrl_d8_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d9_inst_PRR_lag1[itr,j]-precip_44_ctrl_d9_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d10_inst_PRR_lag1[itr,j]-precip_44_ctrl_d10_inst_PR2_lag1[itr,j]

          if ( j >=2 ) :
               precip_44_ctrl_d2_inst_PR2_lag2[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-2]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag2[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-2]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag2[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-2]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag2[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-2]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag2[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-2]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag2[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-2]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag2[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-2]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag2[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-2]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag2[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-2]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d2_inst_PRR_lag2[itr,j]-precip_44_ctrl_d2_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d3_inst_PRR_lag2[itr,j]-precip_44_ctrl_d3_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d4_inst_PRR_lag2[itr,j]-precip_44_ctrl_d4_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d5_inst_PRR_lag2[itr,j]-precip_44_ctrl_d5_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d6_inst_PRR_lag2[itr,j]-precip_44_ctrl_d6_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d7_inst_PRR_lag2[itr,j]-precip_44_ctrl_d7_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d8_inst_PRR_lag2[itr,j]-precip_44_ctrl_d8_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d9_inst_PRR_lag2[itr,j]-precip_44_ctrl_d9_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d10_inst_PRR_lag2[itr,j]-precip_44_ctrl_d10_inst_PR2_lag2[itr,j]

          if ( j >=3 ) :
               precip_44_ctrl_d2_inst_PR2_lag3[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-3]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag3[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-3]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag3[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-3]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag3[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-3]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag3[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-3]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag3[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-3]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag3[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-3]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag3[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-3]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag3[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-3]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d2_inst_PRR_lag3[itr,j]-precip_44_ctrl_d2_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d3_inst_PRR_lag3[itr,j]-precip_44_ctrl_d3_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d4_inst_PRR_lag3[itr,j]-precip_44_ctrl_d4_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d5_inst_PRR_lag3[itr,j]-precip_44_ctrl_d5_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d6_inst_PRR_lag3[itr,j]-precip_44_ctrl_d6_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d7_inst_PRR_lag3[itr,j]-precip_44_ctrl_d7_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d8_inst_PRR_lag3[itr,j]-precip_44_ctrl_d8_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d9_inst_PRR_lag3[itr,j]-precip_44_ctrl_d9_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d10_inst_PRR_lag3[itr,j]-precip_44_ctrl_d10_inst_PR2_lag3[itr,j]

          if ( j >=4 ) :
               precip_44_ctrl_d2_inst_PR2_lag4[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-4]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag4[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-4]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag4[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-4]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag4[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-4]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag4[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-4]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag4[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-4]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag4[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-4]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag4[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-4]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag4[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-4]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d2_inst_PRR_lag4[itr,j]-precip_44_ctrl_d2_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d3_inst_PRR_lag4[itr,j]-precip_44_ctrl_d3_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d4_inst_PRR_lag4[itr,j]-precip_44_ctrl_d4_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d5_inst_PRR_lag4[itr,j]-precip_44_ctrl_d5_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d6_inst_PRR_lag4[itr,j]-precip_44_ctrl_d6_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d7_inst_PRR_lag4[itr,j]-precip_44_ctrl_d7_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d8_inst_PRR_lag4[itr,j]-precip_44_ctrl_d8_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d9_inst_PRR_lag4[itr,j]-precip_44_ctrl_d9_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d10_inst_PRR_lag4[itr,j]-precip_44_ctrl_d10_inst_PR2_lag4[itr,j]

 
          if ( j >=5 ) :
               precip_44_ctrl_d2_inst_PR2_lag5[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-5]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag5[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-5]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag5[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-5]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag5[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-5]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag5[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-5]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag5[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-5]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag5[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-5]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag5[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-5]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag5[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-5]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d2_inst_PRR_lag5[itr,j]-precip_44_ctrl_d2_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d3_inst_PRR_lag5[itr,j]-precip_44_ctrl_d3_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d4_inst_PRR_lag5[itr,j]-precip_44_ctrl_d4_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d5_inst_PRR_lag5[itr,j]-precip_44_ctrl_d5_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d6_inst_PRR_lag5[itr,j]-precip_44_ctrl_d6_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d7_inst_PRR_lag5[itr,j]-precip_44_ctrl_d7_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d8_inst_PRR_lag5[itr,j]-precip_44_ctrl_d8_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d9_inst_PRR_lag5[itr,j]-precip_44_ctrl_d9_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d10_inst_PRR_lag5[itr,j]-precip_44_ctrl_d10_inst_PR2_lag5[itr,j]

          if ( j >=6 ) :
               precip_44_ctrl_d2_inst_PR2_lag6[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-6]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag6[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-6]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag6[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-6]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag6[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-6]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag6[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-6]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag6[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-6]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag6[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-6]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag6[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-6]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag6[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-6]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d2_inst_PRR_lag6[itr,j]-precip_44_ctrl_d2_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d3_inst_PRR_lag6[itr,j]-precip_44_ctrl_d3_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d4_inst_PRR_lag6[itr,j]-precip_44_ctrl_d4_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d5_inst_PRR_lag6[itr,j]-precip_44_ctrl_d5_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d6_inst_PRR_lag6[itr,j]-precip_44_ctrl_d6_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d7_inst_PRR_lag6[itr,j]-precip_44_ctrl_d7_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d8_inst_PRR_lag6[itr,j]-precip_44_ctrl_d8_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d9_inst_PRR_lag6[itr,j]-precip_44_ctrl_d9_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d10_inst_PRR_lag6[itr,j]-precip_44_ctrl_d10_inst_PR2_lag6[itr,j]

          if ( j >=7 ) :
               precip_44_ctrl_d2_inst_PR2_lag7[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-7]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag7[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-7]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag7[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-7]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag7[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-7]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag7[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-7]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag7[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-7]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag7[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-7]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag7[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-7]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag7[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-7]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d2_inst_PRR_lag7[itr,j]-precip_44_ctrl_d2_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d3_inst_PRR_lag7[itr,j]-precip_44_ctrl_d3_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d4_inst_PRR_lag7[itr,j]-precip_44_ctrl_d4_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d5_inst_PRR_lag7[itr,j]-precip_44_ctrl_d5_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d6_inst_PRR_lag7[itr,j]-precip_44_ctrl_d6_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d7_inst_PRR_lag7[itr,j]-precip_44_ctrl_d7_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d8_inst_PRR_lag7[itr,j]-precip_44_ctrl_d8_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d9_inst_PRR_lag7[itr,j]-precip_44_ctrl_d9_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d10_inst_PRR_lag7[itr,j]-precip_44_ctrl_d10_inst_PR2_lag7[itr,j]

          if ( j >=8 ) :
               precip_44_ctrl_d2_inst_PR2_lag8[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-8]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag8[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-8]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag8[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-8]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag8[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-8]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag8[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-8]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag8[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-8]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag8[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-8]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag8[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-8]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag8[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-8]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d2_inst_PRR_lag8[itr,j]-precip_44_ctrl_d2_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d3_inst_PRR_lag8[itr,j]-precip_44_ctrl_d3_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d4_inst_PRR_lag8[itr,j]-precip_44_ctrl_d4_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d5_inst_PRR_lag8[itr,j]-precip_44_ctrl_d5_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d6_inst_PRR_lag8[itr,j]-precip_44_ctrl_d6_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d7_inst_PRR_lag8[itr,j]-precip_44_ctrl_d7_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d8_inst_PRR_lag8[itr,j]-precip_44_ctrl_d8_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d9_inst_PRR_lag8[itr,j]-precip_44_ctrl_d9_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d10_inst_PRR_lag8[itr,j]-precip_44_ctrl_d10_inst_PR2_lag8[itr,j]

          if ( j >=9 ) :
               precip_44_ctrl_d2_inst_PR2_lag9[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-9]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag9[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-9]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag9[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-9]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag9[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-9]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag9[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-9]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag9[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-9]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag9[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-9]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag9[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-9]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag9[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-9]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d2_inst_PRR_lag9[itr,j]-precip_44_ctrl_d2_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d3_inst_PRR_lag9[itr,j]-precip_44_ctrl_d3_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d4_inst_PRR_lag9[itr,j]-precip_44_ctrl_d4_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d5_inst_PRR_lag9[itr,j]-precip_44_ctrl_d5_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d6_inst_PRR_lag9[itr,j]-precip_44_ctrl_d6_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d7_inst_PRR_lag9[itr,j]-precip_44_ctrl_d7_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d8_inst_PRR_lag9[itr,j]-precip_44_ctrl_d8_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d9_inst_PRR_lag9[itr,j]-precip_44_ctrl_d9_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d10_inst_PRR_lag9[itr,j]-precip_44_ctrl_d10_inst_PR2_lag9[itr,j]

          if ( j >=10 ) :
               precip_44_ctrl_d2_inst_PR2_lag10[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-10]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag10[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-10]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag10[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-10]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag10[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-10]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag10[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-10]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag10[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-10]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag10[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-10]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag10[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-10]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag10[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-10]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d2_inst_PRR_lag10[itr,j]-precip_44_ctrl_d2_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d3_inst_PRR_lag10[itr,j]-precip_44_ctrl_d3_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d4_inst_PRR_lag10[itr,j]-precip_44_ctrl_d4_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d5_inst_PRR_lag10[itr,j]-precip_44_ctrl_d5_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d6_inst_PRR_lag10[itr,j]-precip_44_ctrl_d6_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d7_inst_PRR_lag10[itr,j]-precip_44_ctrl_d7_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d8_inst_PRR_lag10[itr,j]-precip_44_ctrl_d8_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d9_inst_PRR_lag10[itr,j]-precip_44_ctrl_d9_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d10_inst_PRR_lag10[itr,j]-precip_44_ctrl_d10_inst_PR2_lag10[itr,j]

          if ( j >=11 ) :
               precip_44_ctrl_d2_inst_PR2_lag11[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-11]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag11[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-11]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag11[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-11]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag11[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-11]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag11[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-11]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag11[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-11]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag11[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-11]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag11[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-11]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag11[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-11]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d2_inst_PRR_lag11[itr,j]-precip_44_ctrl_d2_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d3_inst_PRR_lag11[itr,j]-precip_44_ctrl_d3_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d4_inst_PRR_lag11[itr,j]-precip_44_ctrl_d4_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d5_inst_PRR_lag11[itr,j]-precip_44_ctrl_d5_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d6_inst_PRR_lag11[itr,j]-precip_44_ctrl_d6_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d7_inst_PRR_lag11[itr,j]-precip_44_ctrl_d7_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d8_inst_PRR_lag11[itr,j]-precip_44_ctrl_d8_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d9_inst_PRR_lag11[itr,j]-precip_44_ctrl_d9_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d10_inst_PRR_lag11[itr,j]-precip_44_ctrl_d10_inst_PR2_lag11[itr,j]

          if ( j >=12 ) :
               precip_44_ctrl_d2_inst_PR2_lag12[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-12]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag12[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-12]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag12[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-12]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag12[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-12]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag12[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-12]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag12[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-12]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag12[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-12]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag12[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-12]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag12[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-12]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d2_inst_PRR_lag12[itr,j]-precip_44_ctrl_d2_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d3_inst_PRR_lag12[itr,j]-precip_44_ctrl_d3_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d4_inst_PRR_lag12[itr,j]-precip_44_ctrl_d4_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d5_inst_PRR_lag12[itr,j]-precip_44_ctrl_d5_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d6_inst_PRR_lag12[itr,j]-precip_44_ctrl_d6_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d7_inst_PRR_lag12[itr,j]-precip_44_ctrl_d7_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d8_inst_PRR_lag12[itr,j]-precip_44_ctrl_d8_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d9_inst_PRR_lag12[itr,j]-precip_44_ctrl_d9_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d10_inst_PRR_lag12[itr,j]-precip_44_ctrl_d10_inst_PR2_lag12[itr,j]

          if ( j >=13 ) :
               precip_44_ctrl_d2_inst_PR2_lag13[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-13]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag13[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-13]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag13[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-13]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag13[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-13]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag13[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-13]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag13[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-13]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag13[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-13]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag13[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-13]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag13[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-13]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d2_inst_PRR_lag13[itr,j]-precip_44_ctrl_d2_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d3_inst_PRR_lag13[itr,j]-precip_44_ctrl_d3_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d4_inst_PRR_lag13[itr,j]-precip_44_ctrl_d4_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d5_inst_PRR_lag13[itr,j]-precip_44_ctrl_d5_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d6_inst_PRR_lag13[itr,j]-precip_44_ctrl_d6_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d7_inst_PRR_lag13[itr,j]-precip_44_ctrl_d7_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d8_inst_PRR_lag13[itr,j]-precip_44_ctrl_d8_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d9_inst_PRR_lag13[itr,j]-precip_44_ctrl_d9_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d10_inst_PRR_lag13[itr,j]-precip_44_ctrl_d10_inst_PR2_lag13[itr,j]

          if ( j >=14 ) :
               precip_44_ctrl_d2_inst_PR2_lag14[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-14]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag14[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-14]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag14[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-14]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag14[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-14]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag14[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-14]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag14[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-14]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag14[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-14]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag14[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-14]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag14[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-14]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d2_inst_PRR_lag14[itr,j]-precip_44_ctrl_d2_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d3_inst_PRR_lag14[itr,j]-precip_44_ctrl_d3_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d4_inst_PRR_lag14[itr,j]-precip_44_ctrl_d4_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d5_inst_PRR_lag14[itr,j]-precip_44_ctrl_d5_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d6_inst_PRR_lag14[itr,j]-precip_44_ctrl_d6_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d7_inst_PRR_lag14[itr,j]-precip_44_ctrl_d7_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d8_inst_PRR_lag14[itr,j]-precip_44_ctrl_d8_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d9_inst_PRR_lag14[itr,j]-precip_44_ctrl_d9_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d10_inst_PRR_lag14[itr,j]-precip_44_ctrl_d10_inst_PR2_lag14[itr,j]

          if ( j >=15 ) :
               precip_44_ctrl_d2_inst_PR2_lag15[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-15]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag15[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-15]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag15[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-15]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag15[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-15]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag15[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-15]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag15[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-15]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag15[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-15]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag15[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-15]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag15[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-15]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d2_inst_PRR_lag15[itr,j]-precip_44_ctrl_d2_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d3_inst_PRR_lag15[itr,j]-precip_44_ctrl_d3_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d4_inst_PRR_lag15[itr,j]-precip_44_ctrl_d4_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d5_inst_PRR_lag15[itr,j]-precip_44_ctrl_d5_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d6_inst_PRR_lag15[itr,j]-precip_44_ctrl_d6_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d7_inst_PRR_lag15[itr,j]-precip_44_ctrl_d7_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d8_inst_PRR_lag15[itr,j]-precip_44_ctrl_d8_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d9_inst_PRR_lag15[itr,j]-precip_44_ctrl_d9_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d10_inst_PRR_lag15[itr,j]-precip_44_ctrl_d10_inst_PR2_lag15[itr,j]

          if ( j >=16 ) :
               precip_44_ctrl_d2_inst_PR2_lag16[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-16]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag16[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-16]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag16[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-16]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag16[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-16]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag16[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-16]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag16[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-16]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag16[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-16]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag16[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-16]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag16[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-16]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d2_inst_PRR_lag16[itr,j]-precip_44_ctrl_d2_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d3_inst_PRR_lag16[itr,j]-precip_44_ctrl_d3_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d4_inst_PRR_lag16[itr,j]-precip_44_ctrl_d4_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d5_inst_PRR_lag16[itr,j]-precip_44_ctrl_d5_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d6_inst_PRR_lag16[itr,j]-precip_44_ctrl_d6_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d7_inst_PRR_lag16[itr,j]-precip_44_ctrl_d7_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d8_inst_PRR_lag16[itr,j]-precip_44_ctrl_d8_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d9_inst_PRR_lag16[itr,j]-precip_44_ctrl_d9_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d10_inst_PRR_lag16[itr,j]-precip_44_ctrl_d10_inst_PR2_lag16[itr,j]

          if ( j >=17 ) :
               precip_44_ctrl_d2_inst_PR2_lag17[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-17]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag17[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-17]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag17[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-17]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag17[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-17]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag17[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-17]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag17[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-17]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag17[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-17]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag17[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-17]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag17[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-17]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d2_inst_PRR_lag17[itr,j]-precip_44_ctrl_d2_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d3_inst_PRR_lag17[itr,j]-precip_44_ctrl_d3_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d4_inst_PRR_lag17[itr,j]-precip_44_ctrl_d4_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d5_inst_PRR_lag17[itr,j]-precip_44_ctrl_d5_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d6_inst_PRR_lag17[itr,j]-precip_44_ctrl_d6_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d7_inst_PRR_lag17[itr,j]-precip_44_ctrl_d7_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d8_inst_PRR_lag17[itr,j]-precip_44_ctrl_d8_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d9_inst_PRR_lag17[itr,j]-precip_44_ctrl_d9_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d10_inst_PRR_lag17[itr,j]-precip_44_ctrl_d10_inst_PR2_lag17[itr,j]

          if ( j >=18 ) :
               precip_44_ctrl_d2_inst_PR2_lag18[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-18]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag18[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-18]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag18[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-18]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag18[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-18]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag18[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-18]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag18[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-18]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag18[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-18]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag18[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-18]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag18[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-18]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d2_inst_PRR_lag18[itr,j]-precip_44_ctrl_d2_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d3_inst_PRR_lag18[itr,j]-precip_44_ctrl_d3_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d4_inst_PRR_lag18[itr,j]-precip_44_ctrl_d4_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d5_inst_PRR_lag18[itr,j]-precip_44_ctrl_d5_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d6_inst_PRR_lag18[itr,j]-precip_44_ctrl_d6_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d7_inst_PRR_lag18[itr,j]-precip_44_ctrl_d7_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d8_inst_PRR_lag18[itr,j]-precip_44_ctrl_d8_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d9_inst_PRR_lag18[itr,j]-precip_44_ctrl_d9_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d10_inst_PRR_lag18[itr,j]-precip_44_ctrl_d10_inst_PR2_lag18[itr,j]

          if ( j >=19 ) :
               precip_44_ctrl_d2_inst_PR2_lag19[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-19]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag19[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-19]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag19[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-19]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag19[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-19]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag19[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-19]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag19[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-19]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag19[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-19]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag19[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-19]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag19[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-19]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d2_inst_PRR_lag19[itr,j]-precip_44_ctrl_d2_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d3_inst_PRR_lag19[itr,j]-precip_44_ctrl_d3_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d4_inst_PRR_lag19[itr,j]-precip_44_ctrl_d4_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d5_inst_PRR_lag19[itr,j]-precip_44_ctrl_d5_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d6_inst_PRR_lag19[itr,j]-precip_44_ctrl_d6_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d7_inst_PRR_lag19[itr,j]-precip_44_ctrl_d7_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d8_inst_PRR_lag19[itr,j]-precip_44_ctrl_d8_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d9_inst_PRR_lag19[itr,j]-precip_44_ctrl_d9_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d10_inst_PRR_lag19[itr,j]-precip_44_ctrl_d10_inst_PR2_lag19[itr,j]

          if ( j >=20 ) :
               precip_44_ctrl_d2_inst_PR2_lag20[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-20]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag20[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-20]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag20[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-20]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag20[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-20]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag20[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-20]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag20[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-20]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag20[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-20]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag20[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-20]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag20[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-20]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d2_inst_PRR_lag20[itr,j]-precip_44_ctrl_d2_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d3_inst_PRR_lag20[itr,j]-precip_44_ctrl_d3_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d4_inst_PRR_lag20[itr,j]-precip_44_ctrl_d4_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d5_inst_PRR_lag20[itr,j]-precip_44_ctrl_d5_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d6_inst_PRR_lag20[itr,j]-precip_44_ctrl_d6_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d7_inst_PRR_lag20[itr,j]-precip_44_ctrl_d7_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d8_inst_PRR_lag20[itr,j]-precip_44_ctrl_d8_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d9_inst_PRR_lag20[itr,j]-precip_44_ctrl_d9_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d10_inst_PRR_lag20[itr,j]-precip_44_ctrl_d10_inst_PR2_lag20[itr,j]

          if ( j >=21 ) :
               precip_44_ctrl_d2_inst_PR2_lag21[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-21]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag21[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-21]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag21[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-21]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag21[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-21]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag21[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-21]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag21[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-21]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag21[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-21]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag21[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-21]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag21[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-21]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d2_inst_PRR_lag21[itr,j]-precip_44_ctrl_d2_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d3_inst_PRR_lag21[itr,j]-precip_44_ctrl_d3_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d4_inst_PRR_lag21[itr,j]-precip_44_ctrl_d4_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d5_inst_PRR_lag21[itr,j]-precip_44_ctrl_d5_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d6_inst_PRR_lag21[itr,j]-precip_44_ctrl_d6_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d7_inst_PRR_lag21[itr,j]-precip_44_ctrl_d7_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d8_inst_PRR_lag21[itr,j]-precip_44_ctrl_d8_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d9_inst_PRR_lag21[itr,j]-precip_44_ctrl_d9_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d10_inst_PRR_lag21[itr,j]-precip_44_ctrl_d10_inst_PR2_lag21[itr,j]

          if ( j >=22 ) :
               precip_44_ctrl_d2_inst_PR2_lag22[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-22]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag22[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-22]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag22[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-22]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag22[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-22]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag22[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-22]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag22[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-22]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag22[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-22]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag22[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-22]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag22[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-22]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d2_inst_PRR_lag22[itr,j]-precip_44_ctrl_d2_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d3_inst_PRR_lag22[itr,j]-precip_44_ctrl_d3_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d4_inst_PRR_lag22[itr,j]-precip_44_ctrl_d4_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d5_inst_PRR_lag22[itr,j]-precip_44_ctrl_d5_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d6_inst_PRR_lag22[itr,j]-precip_44_ctrl_d6_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d7_inst_PRR_lag22[itr,j]-precip_44_ctrl_d7_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d8_inst_PRR_lag22[itr,j]-precip_44_ctrl_d8_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d9_inst_PRR_lag22[itr,j]-precip_44_ctrl_d9_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d10_inst_PRR_lag22[itr,j]-precip_44_ctrl_d10_inst_PR2_lag22[itr,j]

          if ( j >=23 ) :
               precip_44_ctrl_d2_inst_PR2_lag23[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-23]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag23[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-23]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag23[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-23]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag23[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-23]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag23[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-23]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag23[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-23]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag23[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-23]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag23[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-23]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag23[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-23]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d2_inst_PRR_lag23[itr,j]-precip_44_ctrl_d2_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d3_inst_PRR_lag23[itr,j]-precip_44_ctrl_d3_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d4_inst_PRR_lag23[itr,j]-precip_44_ctrl_d4_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d5_inst_PRR_lag23[itr,j]-precip_44_ctrl_d5_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d6_inst_PRR_lag23[itr,j]-precip_44_ctrl_d6_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d7_inst_PRR_lag23[itr,j]-precip_44_ctrl_d7_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d8_inst_PRR_lag23[itr,j]-precip_44_ctrl_d8_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d9_inst_PRR_lag23[itr,j]-precip_44_ctrl_d9_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d10_inst_PRR_lag23[itr,j]-precip_44_ctrl_d10_inst_PR2_lag23[itr,j]

          if ( j >=24 ) :
               precip_44_ctrl_d2_inst_PR2_lag24[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-24]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag24[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-24]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag24[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-24]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag24[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-24]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag24[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-24]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag24[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-24]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag24[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-24]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag24[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-24]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag24[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-24]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d2_inst_PRR_lag24[itr,j]-precip_44_ctrl_d2_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d3_inst_PRR_lag24[itr,j]-precip_44_ctrl_d3_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d4_inst_PRR_lag24[itr,j]-precip_44_ctrl_d4_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d5_inst_PRR_lag24[itr,j]-precip_44_ctrl_d5_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d6_inst_PRR_lag24[itr,j]-precip_44_ctrl_d6_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d7_inst_PRR_lag24[itr,j]-precip_44_ctrl_d7_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d8_inst_PRR_lag24[itr,j]-precip_44_ctrl_d8_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d9_inst_PRR_lag24[itr,j]-precip_44_ctrl_d9_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d10_inst_PRR_lag24[itr,j]-precip_44_ctrl_d10_inst_PR2_lag24[itr,j]

          if ( j >=25 ) :
               precip_44_ctrl_d2_inst_PR2_lag25[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-25]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag25[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-25]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag25[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-25]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag25[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-25]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag25[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-25]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag25[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-25]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag25[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-25]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag25[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-25]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag25[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-25]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d2_inst_PRR_lag25[itr,j]-precip_44_ctrl_d2_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d3_inst_PRR_lag25[itr,j]-precip_44_ctrl_d3_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d4_inst_PRR_lag25[itr,j]-precip_44_ctrl_d4_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d5_inst_PRR_lag25[itr,j]-precip_44_ctrl_d5_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d6_inst_PRR_lag25[itr,j]-precip_44_ctrl_d6_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d7_inst_PRR_lag25[itr,j]-precip_44_ctrl_d7_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d8_inst_PRR_lag25[itr,j]-precip_44_ctrl_d8_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d9_inst_PRR_lag25[itr,j]-precip_44_ctrl_d9_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d10_inst_PRR_lag25[itr,j]-precip_44_ctrl_d10_inst_PR2_lag25[itr,j]

          if ( j >=26 ) :
               precip_44_ctrl_d2_inst_PR2_lag26[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-26]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag26[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-26]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag26[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-26]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag26[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-26]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag26[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-26]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag26[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-26]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag26[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-26]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag26[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-26]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag26[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-26]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d2_inst_PRR_lag26[itr,j]-precip_44_ctrl_d2_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d3_inst_PRR_lag26[itr,j]-precip_44_ctrl_d3_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d4_inst_PRR_lag26[itr,j]-precip_44_ctrl_d4_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d5_inst_PRR_lag26[itr,j]-precip_44_ctrl_d5_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d6_inst_PRR_lag26[itr,j]-precip_44_ctrl_d6_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d7_inst_PRR_lag26[itr,j]-precip_44_ctrl_d7_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d8_inst_PRR_lag26[itr,j]-precip_44_ctrl_d8_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d9_inst_PRR_lag26[itr,j]-precip_44_ctrl_d9_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d10_inst_PRR_lag26[itr,j]-precip_44_ctrl_d10_inst_PR2_lag26[itr,j]

          if ( j >=27 ) :
               precip_44_ctrl_d2_inst_PR2_lag27[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-27]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag27[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-27]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag27[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-27]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag27[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-27]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag27[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-27]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag27[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-27]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag27[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-27]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag27[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-27]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag27[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-27]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d2_inst_PRR_lag27[itr,j]-precip_44_ctrl_d2_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d3_inst_PRR_lag27[itr,j]-precip_44_ctrl_d3_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d4_inst_PRR_lag27[itr,j]-precip_44_ctrl_d4_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d5_inst_PRR_lag27[itr,j]-precip_44_ctrl_d5_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d6_inst_PRR_lag27[itr,j]-precip_44_ctrl_d6_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d7_inst_PRR_lag27[itr,j]-precip_44_ctrl_d7_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d8_inst_PRR_lag27[itr,j]-precip_44_ctrl_d8_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d9_inst_PRR_lag27[itr,j]-precip_44_ctrl_d9_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d10_inst_PRR_lag27[itr,j]-precip_44_ctrl_d10_inst_PR2_lag27[itr,j]

          if ( j >=28 ) :
               precip_44_ctrl_d2_inst_PR2_lag28[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-28]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag28[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-28]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag28[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-28]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag28[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-28]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag28[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-28]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag28[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-28]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag28[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-28]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag28[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-28]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag28[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-28]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d2_inst_PRR_lag28[itr,j]-precip_44_ctrl_d2_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d3_inst_PRR_lag28[itr,j]-precip_44_ctrl_d3_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d4_inst_PRR_lag28[itr,j]-precip_44_ctrl_d4_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d5_inst_PRR_lag28[itr,j]-precip_44_ctrl_d5_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d6_inst_PRR_lag28[itr,j]-precip_44_ctrl_d6_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d7_inst_PRR_lag28[itr,j]-precip_44_ctrl_d7_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d8_inst_PRR_lag28[itr,j]-precip_44_ctrl_d8_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d9_inst_PRR_lag28[itr,j]-precip_44_ctrl_d9_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d10_inst_PRR_lag28[itr,j]-precip_44_ctrl_d10_inst_PR2_lag28[itr,j]

          if ( j >=29 ) :
               precip_44_ctrl_d2_inst_PR2_lag29[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-29]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag29[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-29]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag29[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-29]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag29[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-29]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag29[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-29]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag29[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-29]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag29[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-29]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag29[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-29]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag29[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-29]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d2_inst_PRR_lag29[itr,j]-precip_44_ctrl_d2_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d3_inst_PRR_lag29[itr,j]-precip_44_ctrl_d3_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d4_inst_PRR_lag29[itr,j]-precip_44_ctrl_d4_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d5_inst_PRR_lag29[itr,j]-precip_44_ctrl_d5_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d6_inst_PRR_lag29[itr,j]-precip_44_ctrl_d6_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d7_inst_PRR_lag29[itr,j]-precip_44_ctrl_d7_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d8_inst_PRR_lag29[itr,j]-precip_44_ctrl_d8_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d9_inst_PRR_lag29[itr,j]-precip_44_ctrl_d9_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d10_inst_PRR_lag29[itr,j]-precip_44_ctrl_d10_inst_PR2_lag29[itr,j]

          if ( j >=30 ) :
               precip_44_ctrl_d2_inst_PR2_lag30[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-30]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag30[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-30]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag30[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-30]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag30[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-30]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag30[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-30]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag30[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-30]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag30[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-30]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag30[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-30]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag30[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-30]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d2_inst_PRR_lag30[itr,j]-precip_44_ctrl_d2_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d3_inst_PRR_lag30[itr,j]-precip_44_ctrl_d3_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d4_inst_PRR_lag30[itr,j]-precip_44_ctrl_d4_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d5_inst_PRR_lag30[itr,j]-precip_44_ctrl_d5_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d6_inst_PRR_lag30[itr,j]-precip_44_ctrl_d6_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d7_inst_PRR_lag30[itr,j]-precip_44_ctrl_d7_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d8_inst_PRR_lag30[itr,j]-precip_44_ctrl_d8_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d9_inst_PRR_lag30[itr,j]-precip_44_ctrl_d9_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d10_inst_PRR_lag30[itr,j]-precip_44_ctrl_d10_inst_PR2_lag30[itr,j]

          if ( j >=31 ) :
               precip_44_ctrl_d2_inst_PR2_lag31[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-31]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag31[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-31]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag31[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-31]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag31[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-31]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag31[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-31]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag31[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-31]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag31[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-31]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag31[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-31]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag31[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-31]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d2_inst_PRR_lag31[itr,j]-precip_44_ctrl_d2_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d3_inst_PRR_lag31[itr,j]-precip_44_ctrl_d3_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d4_inst_PRR_lag31[itr,j]-precip_44_ctrl_d4_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d5_inst_PRR_lag31[itr,j]-precip_44_ctrl_d5_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d6_inst_PRR_lag31[itr,j]-precip_44_ctrl_d6_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d7_inst_PRR_lag31[itr,j]-precip_44_ctrl_d7_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d8_inst_PRR_lag31[itr,j]-precip_44_ctrl_d8_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d9_inst_PRR_lag31[itr,j]-precip_44_ctrl_d9_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d10_inst_PRR_lag31[itr,j]-precip_44_ctrl_d10_inst_PR2_lag31[itr,j]

          if ( j >=32 ) :
               precip_44_ctrl_d2_inst_PR2_lag32[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-32]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag32[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-32]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag32[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-32]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag32[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-32]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag32[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-32]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag32[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-32]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag32[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-32]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag32[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-32]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag32[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-32]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d2_inst_PRR_lag32[itr,j]-precip_44_ctrl_d2_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d3_inst_PRR_lag32[itr,j]-precip_44_ctrl_d3_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d4_inst_PRR_lag32[itr,j]-precip_44_ctrl_d4_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d5_inst_PRR_lag32[itr,j]-precip_44_ctrl_d5_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d6_inst_PRR_lag32[itr,j]-precip_44_ctrl_d6_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d7_inst_PRR_lag32[itr,j]-precip_44_ctrl_d7_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d8_inst_PRR_lag32[itr,j]-precip_44_ctrl_d8_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d9_inst_PRR_lag32[itr,j]-precip_44_ctrl_d9_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d10_inst_PRR_lag32[itr,j]-precip_44_ctrl_d10_inst_PR2_lag32[itr,j]


          if ( j >=33 ) :
               precip_44_ctrl_d2_inst_PR2_lag33[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-33]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag33[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-33]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag33[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-33]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag33[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-33]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag33[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-33]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag33[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-33]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag33[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-33]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag33[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-33]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag33[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-33]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag33[itr,j]=precip_44_ctrl_d2_inst_PRR_lag33[itr,j]-precip_44_ctrl_d2_inst_PR2_lag33[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag33[itr,j]=precip_44_ctrl_d3_inst_PRR_lag33[itr,j]-precip_44_ctrl_d3_inst_PR2_lag33[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag33[itr,j]=precip_44_ctrl_d4_inst_PRR_lag33[itr,j]-precip_44_ctrl_d4_inst_PR2_lag33[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag33[itr,j]=precip_44_ctrl_d5_inst_PRR_lag33[itr,j]-precip_44_ctrl_d5_inst_PR2_lag33[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag33[itr,j]=precip_44_ctrl_d6_inst_PRR_lag33[itr,j]-precip_44_ctrl_d6_inst_PR2_lag33[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag33[itr,j]=precip_44_ctrl_d7_inst_PRR_lag33[itr,j]-precip_44_ctrl_d7_inst_PR2_lag33[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag33[itr,j]=precip_44_ctrl_d8_inst_PRR_lag33[itr,j]-precip_44_ctrl_d8_inst_PR2_lag33[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag33[itr,j]=precip_44_ctrl_d9_inst_PRR_lag33[itr,j]-precip_44_ctrl_d9_inst_PR2_lag33[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag33[itr,j]=precip_44_ctrl_d10_inst_PRR_lag33[itr,j]-precip_44_ctrl_d10_inst_PR2_lag33[itr,j]


          if ( j >=34 ) :
               precip_44_ctrl_d2_inst_PR2_lag34[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-34]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag34[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-34]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag34[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-34]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag34[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-34]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag34[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-34]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag34[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-34]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag34[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-34]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag34[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-34]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag34[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-34]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag34[itr,j]=precip_44_ctrl_d2_inst_PRR_lag34[itr,j]-precip_44_ctrl_d2_inst_PR2_lag34[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag34[itr,j]=precip_44_ctrl_d3_inst_PRR_lag34[itr,j]-precip_44_ctrl_d3_inst_PR2_lag34[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag34[itr,j]=precip_44_ctrl_d4_inst_PRR_lag34[itr,j]-precip_44_ctrl_d4_inst_PR2_lag34[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag34[itr,j]=precip_44_ctrl_d5_inst_PRR_lag34[itr,j]-precip_44_ctrl_d5_inst_PR2_lag34[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag34[itr,j]=precip_44_ctrl_d6_inst_PRR_lag34[itr,j]-precip_44_ctrl_d6_inst_PR2_lag34[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag34[itr,j]=precip_44_ctrl_d7_inst_PRR_lag34[itr,j]-precip_44_ctrl_d7_inst_PR2_lag34[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag34[itr,j]=precip_44_ctrl_d8_inst_PRR_lag34[itr,j]-precip_44_ctrl_d8_inst_PR2_lag34[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag34[itr,j]=precip_44_ctrl_d9_inst_PRR_lag34[itr,j]-precip_44_ctrl_d9_inst_PR2_lag34[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag34[itr,j]=precip_44_ctrl_d10_inst_PRR_lag34[itr,j]-precip_44_ctrl_d10_inst_PR2_lag34[itr,j]

          if ( j >=35 ) :
               precip_44_ctrl_d2_inst_PR2_lag35[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-35]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag35[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-35]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag35[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-35]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag35[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-35]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag35[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-35]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag35[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-35]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag35[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-35]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag35[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-35]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag35[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-35]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag35[itr,j]=precip_44_ctrl_d2_inst_PRR_lag35[itr,j]-precip_44_ctrl_d2_inst_PR2_lag35[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag35[itr,j]=precip_44_ctrl_d3_inst_PRR_lag35[itr,j]-precip_44_ctrl_d3_inst_PR2_lag35[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag35[itr,j]=precip_44_ctrl_d4_inst_PRR_lag35[itr,j]-precip_44_ctrl_d4_inst_PR2_lag35[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag35[itr,j]=precip_44_ctrl_d5_inst_PRR_lag35[itr,j]-precip_44_ctrl_d5_inst_PR2_lag35[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag35[itr,j]=precip_44_ctrl_d6_inst_PRR_lag35[itr,j]-precip_44_ctrl_d6_inst_PR2_lag35[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag35[itr,j]=precip_44_ctrl_d7_inst_PRR_lag35[itr,j]-precip_44_ctrl_d7_inst_PR2_lag35[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag35[itr,j]=precip_44_ctrl_d8_inst_PRR_lag35[itr,j]-precip_44_ctrl_d8_inst_PR2_lag35[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag35[itr,j]=precip_44_ctrl_d9_inst_PRR_lag35[itr,j]-precip_44_ctrl_d9_inst_PR2_lag35[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag35[itr,j]=precip_44_ctrl_d10_inst_PRR_lag35[itr,j]-precip_44_ctrl_d10_inst_PR2_lag35[itr,j]

          if ( j >=36 ) :
               precip_44_ctrl_d2_inst_PR2_lag36[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-36]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag36[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-36]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag36[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-36]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag36[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-36]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag36[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-36]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag36[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-36]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag36[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-36]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag36[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-36]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag36[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-36]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag36[itr,j]=precip_44_ctrl_d2_inst_PRR_lag36[itr,j]-precip_44_ctrl_d2_inst_PR2_lag36[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag36[itr,j]=precip_44_ctrl_d3_inst_PRR_lag36[itr,j]-precip_44_ctrl_d3_inst_PR2_lag36[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag36[itr,j]=precip_44_ctrl_d4_inst_PRR_lag36[itr,j]-precip_44_ctrl_d4_inst_PR2_lag36[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag36[itr,j]=precip_44_ctrl_d5_inst_PRR_lag36[itr,j]-precip_44_ctrl_d5_inst_PR2_lag36[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag36[itr,j]=precip_44_ctrl_d6_inst_PRR_lag36[itr,j]-precip_44_ctrl_d6_inst_PR2_lag36[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag36[itr,j]=precip_44_ctrl_d7_inst_PRR_lag36[itr,j]-precip_44_ctrl_d7_inst_PR2_lag36[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag36[itr,j]=precip_44_ctrl_d8_inst_PRR_lag36[itr,j]-precip_44_ctrl_d8_inst_PR2_lag36[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag36[itr,j]=precip_44_ctrl_d9_inst_PRR_lag36[itr,j]-precip_44_ctrl_d9_inst_PR2_lag36[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag36[itr,j]=precip_44_ctrl_d10_inst_PRR_lag36[itr,j]-precip_44_ctrl_d10_inst_PR2_lag36[itr,j]

          if ( j >=37 ) :
               precip_44_ctrl_d2_inst_PR2_lag37[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-37]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag37[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-37]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag37[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-37]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag37[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-37]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag37[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-37]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag37[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-37]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag37[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-37]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag37[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-37]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag37[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-37]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag37[itr,j]=precip_44_ctrl_d2_inst_PRR_lag37[itr,j]-precip_44_ctrl_d2_inst_PR2_lag37[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag37[itr,j]=precip_44_ctrl_d3_inst_PRR_lag37[itr,j]-precip_44_ctrl_d3_inst_PR2_lag37[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag37[itr,j]=precip_44_ctrl_d4_inst_PRR_lag37[itr,j]-precip_44_ctrl_d4_inst_PR2_lag37[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag37[itr,j]=precip_44_ctrl_d5_inst_PRR_lag37[itr,j]-precip_44_ctrl_d5_inst_PR2_lag37[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag37[itr,j]=precip_44_ctrl_d6_inst_PRR_lag37[itr,j]-precip_44_ctrl_d6_inst_PR2_lag37[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag37[itr,j]=precip_44_ctrl_d7_inst_PRR_lag37[itr,j]-precip_44_ctrl_d7_inst_PR2_lag37[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag37[itr,j]=precip_44_ctrl_d8_inst_PRR_lag37[itr,j]-precip_44_ctrl_d8_inst_PR2_lag37[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag37[itr,j]=precip_44_ctrl_d9_inst_PRR_lag37[itr,j]-precip_44_ctrl_d9_inst_PR2_lag37[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag37[itr,j]=precip_44_ctrl_d10_inst_PRR_lag37[itr,j]-precip_44_ctrl_d10_inst_PR2_lag37[itr,j]

          if ( j >=38 ) :
               precip_44_ctrl_d2_inst_PR2_lag38[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-38]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag38[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-38]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag38[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-38]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag38[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-38]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag38[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-38]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag38[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-38]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag38[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-38]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag38[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-38]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag38[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-38]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag38[itr,j]=precip_44_ctrl_d2_inst_PRR_lag38[itr,j]-precip_44_ctrl_d2_inst_PR2_lag38[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag38[itr,j]=precip_44_ctrl_d3_inst_PRR_lag38[itr,j]-precip_44_ctrl_d3_inst_PR2_lag38[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag38[itr,j]=precip_44_ctrl_d4_inst_PRR_lag38[itr,j]-precip_44_ctrl_d4_inst_PR2_lag38[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag38[itr,j]=precip_44_ctrl_d5_inst_PRR_lag38[itr,j]-precip_44_ctrl_d5_inst_PR2_lag38[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag38[itr,j]=precip_44_ctrl_d6_inst_PRR_lag38[itr,j]-precip_44_ctrl_d6_inst_PR2_lag38[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag38[itr,j]=precip_44_ctrl_d7_inst_PRR_lag38[itr,j]-precip_44_ctrl_d7_inst_PR2_lag38[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag38[itr,j]=precip_44_ctrl_d8_inst_PRR_lag38[itr,j]-precip_44_ctrl_d8_inst_PR2_lag38[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag38[itr,j]=precip_44_ctrl_d9_inst_PRR_lag38[itr,j]-precip_44_ctrl_d9_inst_PR2_lag38[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag38[itr,j]=precip_44_ctrl_d10_inst_PRR_lag38[itr,j]-precip_44_ctrl_d10_inst_PR2_lag38[itr,j]

          if ( j >=39 ) :
               precip_44_ctrl_d2_inst_PR2_lag39[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-39]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag39[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-39]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag39[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-39]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag39[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-39]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag39[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-39]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag39[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-39]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag39[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-39]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag39[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-39]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag39[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-39]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag39[itr,j]=precip_44_ctrl_d2_inst_PRR_lag39[itr,j]-precip_44_ctrl_d2_inst_PR2_lag39[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag39[itr,j]=precip_44_ctrl_d3_inst_PRR_lag39[itr,j]-precip_44_ctrl_d3_inst_PR2_lag39[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag39[itr,j]=precip_44_ctrl_d4_inst_PRR_lag39[itr,j]-precip_44_ctrl_d4_inst_PR2_lag39[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag39[itr,j]=precip_44_ctrl_d5_inst_PRR_lag39[itr,j]-precip_44_ctrl_d5_inst_PR2_lag39[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag39[itr,j]=precip_44_ctrl_d6_inst_PRR_lag39[itr,j]-precip_44_ctrl_d6_inst_PR2_lag39[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag39[itr,j]=precip_44_ctrl_d7_inst_PRR_lag39[itr,j]-precip_44_ctrl_d7_inst_PR2_lag39[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag39[itr,j]=precip_44_ctrl_d8_inst_PRR_lag39[itr,j]-precip_44_ctrl_d8_inst_PR2_lag39[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag39[itr,j]=precip_44_ctrl_d9_inst_PRR_lag39[itr,j]-precip_44_ctrl_d9_inst_PR2_lag39[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag39[itr,j]=precip_44_ctrl_d10_inst_PRR_lag39[itr,j]-precip_44_ctrl_d10_inst_PR2_lag39[itr,j]

          if ( j >=40 ) :
               precip_44_ctrl_d2_inst_PR2_lag40[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-40]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag40[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-40]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag40[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-40]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag40[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-40]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag40[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-40]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag40[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-40]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag40[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-40]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag40[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-40]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag40[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-40]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag40[itr,j]=precip_44_ctrl_d2_inst_PRR_lag40[itr,j]-precip_44_ctrl_d2_inst_PR2_lag40[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag40[itr,j]=precip_44_ctrl_d3_inst_PRR_lag40[itr,j]-precip_44_ctrl_d3_inst_PR2_lag40[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag40[itr,j]=precip_44_ctrl_d4_inst_PRR_lag40[itr,j]-precip_44_ctrl_d4_inst_PR2_lag40[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag40[itr,j]=precip_44_ctrl_d5_inst_PRR_lag40[itr,j]-precip_44_ctrl_d5_inst_PR2_lag40[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag40[itr,j]=precip_44_ctrl_d6_inst_PRR_lag40[itr,j]-precip_44_ctrl_d6_inst_PR2_lag40[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag40[itr,j]=precip_44_ctrl_d7_inst_PRR_lag40[itr,j]-precip_44_ctrl_d7_inst_PR2_lag40[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag40[itr,j]=precip_44_ctrl_d8_inst_PRR_lag40[itr,j]-precip_44_ctrl_d8_inst_PR2_lag40[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag40[itr,j]=precip_44_ctrl_d9_inst_PRR_lag40[itr,j]-precip_44_ctrl_d9_inst_PR2_lag40[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag40[itr,j]=precip_44_ctrl_d10_inst_PRR_lag40[itr,j]-precip_44_ctrl_d10_inst_PR2_lag40[itr,j]

          if ( j >=41 ) :
               precip_44_ctrl_d2_inst_PR2_lag41[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-41]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag41[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-41]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag41[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-41]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag41[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-41]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag41[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-41]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag41[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-41]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag41[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-41]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag41[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-41]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag41[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-41]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag41[itr,j]=precip_44_ctrl_d2_inst_PRR_lag41[itr,j]-precip_44_ctrl_d2_inst_PR2_lag41[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag41[itr,j]=precip_44_ctrl_d3_inst_PRR_lag41[itr,j]-precip_44_ctrl_d3_inst_PR2_lag41[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag41[itr,j]=precip_44_ctrl_d4_inst_PRR_lag41[itr,j]-precip_44_ctrl_d4_inst_PR2_lag41[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag41[itr,j]=precip_44_ctrl_d5_inst_PRR_lag41[itr,j]-precip_44_ctrl_d5_inst_PR2_lag41[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag41[itr,j]=precip_44_ctrl_d6_inst_PRR_lag41[itr,j]-precip_44_ctrl_d6_inst_PR2_lag41[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag41[itr,j]=precip_44_ctrl_d7_inst_PRR_lag41[itr,j]-precip_44_ctrl_d7_inst_PR2_lag41[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag41[itr,j]=precip_44_ctrl_d8_inst_PRR_lag41[itr,j]-precip_44_ctrl_d8_inst_PR2_lag41[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag41[itr,j]=precip_44_ctrl_d9_inst_PRR_lag41[itr,j]-precip_44_ctrl_d9_inst_PR2_lag41[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag41[itr,j]=precip_44_ctrl_d10_inst_PRR_lag41[itr,j]-precip_44_ctrl_d10_inst_PR2_lag41[itr,j]

          if ( j >=42 ) :
               precip_44_ctrl_d2_inst_PR2_lag42[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-42]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag42[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-42]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag42[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-42]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag42[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-42]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag42[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-42]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag42[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-42]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag42[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-42]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag42[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-42]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag42[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-42]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag42[itr,j]=precip_44_ctrl_d2_inst_PRR_lag42[itr,j]-precip_44_ctrl_d2_inst_PR2_lag42[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag42[itr,j]=precip_44_ctrl_d3_inst_PRR_lag42[itr,j]-precip_44_ctrl_d3_inst_PR2_lag42[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag42[itr,j]=precip_44_ctrl_d4_inst_PRR_lag42[itr,j]-precip_44_ctrl_d4_inst_PR2_lag42[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag42[itr,j]=precip_44_ctrl_d5_inst_PRR_lag42[itr,j]-precip_44_ctrl_d5_inst_PR2_lag42[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag42[itr,j]=precip_44_ctrl_d6_inst_PRR_lag42[itr,j]-precip_44_ctrl_d6_inst_PR2_lag42[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag42[itr,j]=precip_44_ctrl_d7_inst_PRR_lag42[itr,j]-precip_44_ctrl_d7_inst_PR2_lag42[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag42[itr,j]=precip_44_ctrl_d8_inst_PRR_lag42[itr,j]-precip_44_ctrl_d8_inst_PR2_lag42[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag42[itr,j]=precip_44_ctrl_d9_inst_PRR_lag42[itr,j]-precip_44_ctrl_d9_inst_PR2_lag42[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag42[itr,j]=precip_44_ctrl_d10_inst_PRR_lag42[itr,j]-precip_44_ctrl_d10_inst_PR2_lag42[itr,j]

          if ( j >=43 ) :
               precip_44_ctrl_d2_inst_PR2_lag43[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-43]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag43[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-43]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag43[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-43]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag43[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-43]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag43[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-43]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag43[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-43]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag43[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-43]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag43[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-43]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag43[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-43]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag43[itr,j]=precip_44_ctrl_d2_inst_PRR_lag43[itr,j]-precip_44_ctrl_d2_inst_PR2_lag43[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag43[itr,j]=precip_44_ctrl_d3_inst_PRR_lag43[itr,j]-precip_44_ctrl_d3_inst_PR2_lag43[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag43[itr,j]=precip_44_ctrl_d4_inst_PRR_lag43[itr,j]-precip_44_ctrl_d4_inst_PR2_lag43[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag43[itr,j]=precip_44_ctrl_d5_inst_PRR_lag43[itr,j]-precip_44_ctrl_d5_inst_PR2_lag43[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag43[itr,j]=precip_44_ctrl_d6_inst_PRR_lag43[itr,j]-precip_44_ctrl_d6_inst_PR2_lag43[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag43[itr,j]=precip_44_ctrl_d7_inst_PRR_lag43[itr,j]-precip_44_ctrl_d7_inst_PR2_lag43[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag43[itr,j]=precip_44_ctrl_d8_inst_PRR_lag43[itr,j]-precip_44_ctrl_d8_inst_PR2_lag43[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag43[itr,j]=precip_44_ctrl_d9_inst_PRR_lag43[itr,j]-precip_44_ctrl_d9_inst_PR2_lag43[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag43[itr,j]=precip_44_ctrl_d10_inst_PRR_lag43[itr,j]-precip_44_ctrl_d10_inst_PR2_lag43[itr,j]

          if ( j >=44 ) :
               precip_44_ctrl_d2_inst_PR2_lag44[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-44]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag44[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-44]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag44[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-44]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag44[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-44]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag44[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-44]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag44[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-44]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag44[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-44]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag44[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-44]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag44[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-44]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag44[itr,j]=precip_44_ctrl_d2_inst_PRR_lag44[itr,j]-precip_44_ctrl_d2_inst_PR2_lag44[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag44[itr,j]=precip_44_ctrl_d3_inst_PRR_lag44[itr,j]-precip_44_ctrl_d3_inst_PR2_lag44[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag44[itr,j]=precip_44_ctrl_d4_inst_PRR_lag44[itr,j]-precip_44_ctrl_d4_inst_PR2_lag44[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag44[itr,j]=precip_44_ctrl_d5_inst_PRR_lag44[itr,j]-precip_44_ctrl_d5_inst_PR2_lag44[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag44[itr,j]=precip_44_ctrl_d6_inst_PRR_lag44[itr,j]-precip_44_ctrl_d6_inst_PR2_lag44[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag44[itr,j]=precip_44_ctrl_d7_inst_PRR_lag44[itr,j]-precip_44_ctrl_d7_inst_PR2_lag44[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag44[itr,j]=precip_44_ctrl_d8_inst_PRR_lag44[itr,j]-precip_44_ctrl_d8_inst_PR2_lag44[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag44[itr,j]=precip_44_ctrl_d9_inst_PRR_lag44[itr,j]-precip_44_ctrl_d9_inst_PR2_lag44[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag44[itr,j]=precip_44_ctrl_d10_inst_PRR_lag44[itr,j]-precip_44_ctrl_d10_inst_PR2_lag44[itr,j]

          if ( j >=45 ) :
               precip_44_ctrl_d2_inst_PR2_lag45[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-45]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag45[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-45]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag45[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-45]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag45[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-45]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag45[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-45]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag45[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-45]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag45[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-45]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag45[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-45]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag45[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-45]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag45[itr,j]=precip_44_ctrl_d2_inst_PRR_lag45[itr,j]-precip_44_ctrl_d2_inst_PR2_lag45[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag45[itr,j]=precip_44_ctrl_d3_inst_PRR_lag45[itr,j]-precip_44_ctrl_d3_inst_PR2_lag45[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag45[itr,j]=precip_44_ctrl_d4_inst_PRR_lag45[itr,j]-precip_44_ctrl_d4_inst_PR2_lag45[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag45[itr,j]=precip_44_ctrl_d5_inst_PRR_lag45[itr,j]-precip_44_ctrl_d5_inst_PR2_lag45[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag45[itr,j]=precip_44_ctrl_d6_inst_PRR_lag45[itr,j]-precip_44_ctrl_d6_inst_PR2_lag45[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag45[itr,j]=precip_44_ctrl_d7_inst_PRR_lag45[itr,j]-precip_44_ctrl_d7_inst_PR2_lag45[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag45[itr,j]=precip_44_ctrl_d8_inst_PRR_lag45[itr,j]-precip_44_ctrl_d8_inst_PR2_lag45[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag45[itr,j]=precip_44_ctrl_d9_inst_PRR_lag45[itr,j]-precip_44_ctrl_d9_inst_PR2_lag45[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag45[itr,j]=precip_44_ctrl_d10_inst_PRR_lag45[itr,j]-precip_44_ctrl_d10_inst_PR2_lag45[itr,j]

          if ( j >=46 ) :
               precip_44_ctrl_d2_inst_PR2_lag46[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-46]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag46[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-46]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag46[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-46]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag46[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-46]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag46[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-46]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag46[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-46]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag46[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-46]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag46[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-46]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag46[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-46]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag46[itr,j]=precip_44_ctrl_d2_inst_PRR_lag46[itr,j]-precip_44_ctrl_d2_inst_PR2_lag46[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag46[itr,j]=precip_44_ctrl_d3_inst_PRR_lag46[itr,j]-precip_44_ctrl_d3_inst_PR2_lag46[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag46[itr,j]=precip_44_ctrl_d4_inst_PRR_lag46[itr,j]-precip_44_ctrl_d4_inst_PR2_lag46[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag46[itr,j]=precip_44_ctrl_d5_inst_PRR_lag46[itr,j]-precip_44_ctrl_d5_inst_PR2_lag46[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag46[itr,j]=precip_44_ctrl_d6_inst_PRR_lag46[itr,j]-precip_44_ctrl_d6_inst_PR2_lag46[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag46[itr,j]=precip_44_ctrl_d7_inst_PRR_lag46[itr,j]-precip_44_ctrl_d7_inst_PR2_lag46[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag46[itr,j]=precip_44_ctrl_d8_inst_PRR_lag46[itr,j]-precip_44_ctrl_d8_inst_PR2_lag46[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag46[itr,j]=precip_44_ctrl_d9_inst_PRR_lag46[itr,j]-precip_44_ctrl_d9_inst_PR2_lag46[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag46[itr,j]=precip_44_ctrl_d10_inst_PRR_lag46[itr,j]-precip_44_ctrl_d10_inst_PR2_lag46[itr,j]

          if ( j >=47 ) :
               precip_44_ctrl_d2_inst_PR2_lag47[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-47]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag47[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-47]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag47[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-47]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag47[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-47]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag47[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-47]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag47[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-47]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag47[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-47]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag47[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-47]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag47[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-47]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag47[itr,j]=precip_44_ctrl_d2_inst_PRR_lag47[itr,j]-precip_44_ctrl_d2_inst_PR2_lag47[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag47[itr,j]=precip_44_ctrl_d3_inst_PRR_lag47[itr,j]-precip_44_ctrl_d3_inst_PR2_lag47[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag47[itr,j]=precip_44_ctrl_d4_inst_PRR_lag47[itr,j]-precip_44_ctrl_d4_inst_PR2_lag47[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag47[itr,j]=precip_44_ctrl_d5_inst_PRR_lag47[itr,j]-precip_44_ctrl_d5_inst_PR2_lag47[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag47[itr,j]=precip_44_ctrl_d6_inst_PRR_lag47[itr,j]-precip_44_ctrl_d6_inst_PR2_lag47[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag47[itr,j]=precip_44_ctrl_d7_inst_PRR_lag47[itr,j]-precip_44_ctrl_d7_inst_PR2_lag47[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag47[itr,j]=precip_44_ctrl_d8_inst_PRR_lag47[itr,j]-precip_44_ctrl_d8_inst_PR2_lag47[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag47[itr,j]=precip_44_ctrl_d9_inst_PRR_lag47[itr,j]-precip_44_ctrl_d9_inst_PR2_lag47[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag47[itr,j]=precip_44_ctrl_d10_inst_PRR_lag47[itr,j]-precip_44_ctrl_d10_inst_PR2_lag47[itr,j]

          if ( j >=48 ) :
               precip_44_ctrl_d2_inst_PR2_lag48[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-48]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag48[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-48]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag48[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-48]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag48[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-48]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag48[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-48]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag48[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-48]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag48[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-48]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag48[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-48]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag48[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-48]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag48[itr,j]=precip_44_ctrl_d2_inst_PRR_lag48[itr,j]-precip_44_ctrl_d2_inst_PR2_lag48[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag48[itr,j]=precip_44_ctrl_d3_inst_PRR_lag48[itr,j]-precip_44_ctrl_d3_inst_PR2_lag48[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag48[itr,j]=precip_44_ctrl_d4_inst_PRR_lag48[itr,j]-precip_44_ctrl_d4_inst_PR2_lag48[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag48[itr,j]=precip_44_ctrl_d5_inst_PRR_lag48[itr,j]-precip_44_ctrl_d5_inst_PR2_lag48[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag48[itr,j]=precip_44_ctrl_d6_inst_PRR_lag48[itr,j]-precip_44_ctrl_d6_inst_PR2_lag48[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag48[itr,j]=precip_44_ctrl_d7_inst_PRR_lag48[itr,j]-precip_44_ctrl_d7_inst_PR2_lag48[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag48[itr,j]=precip_44_ctrl_d8_inst_PRR_lag48[itr,j]-precip_44_ctrl_d8_inst_PR2_lag48[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag48[itr,j]=precip_44_ctrl_d9_inst_PRR_lag48[itr,j]-precip_44_ctrl_d9_inst_PR2_lag48[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag48[itr,j]=precip_44_ctrl_d10_inst_PRR_lag48[itr,j]-precip_44_ctrl_d10_inst_PR2_lag48[itr,j]

          if ( j >=49 ) :
               precip_44_ctrl_d2_inst_PR2_lag49[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-49]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag49[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-49]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag49[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-49]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag49[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-49]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag49[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-49]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag49[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-49]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag49[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-49]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag49[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-49]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag49[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-49]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag49[itr,j]=precip_44_ctrl_d2_inst_PRR_lag49[itr,j]-precip_44_ctrl_d2_inst_PR2_lag49[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag49[itr,j]=precip_44_ctrl_d3_inst_PRR_lag49[itr,j]-precip_44_ctrl_d3_inst_PR2_lag49[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag49[itr,j]=precip_44_ctrl_d4_inst_PRR_lag49[itr,j]-precip_44_ctrl_d4_inst_PR2_lag49[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag49[itr,j]=precip_44_ctrl_d5_inst_PRR_lag49[itr,j]-precip_44_ctrl_d5_inst_PR2_lag49[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag49[itr,j]=precip_44_ctrl_d6_inst_PRR_lag49[itr,j]-precip_44_ctrl_d6_inst_PR2_lag49[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag49[itr,j]=precip_44_ctrl_d7_inst_PRR_lag49[itr,j]-precip_44_ctrl_d7_inst_PR2_lag49[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag49[itr,j]=precip_44_ctrl_d8_inst_PRR_lag49[itr,j]-precip_44_ctrl_d8_inst_PR2_lag49[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag49[itr,j]=precip_44_ctrl_d9_inst_PRR_lag49[itr,j]-precip_44_ctrl_d9_inst_PR2_lag49[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag49[itr,j]=precip_44_ctrl_d10_inst_PRR_lag49[itr,j]-precip_44_ctrl_d10_inst_PR2_lag49[itr,j]

          proR_lag_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR_lag[itr,j]+precip_44_ctrl_d3_inst_PR_lag[itr,j]+precip_44_ctrl_d4_inst_PR_lag[itr,j]+precip_44_ctrl_d5_inst_PR_lag[itr,j]+precip_44_ctrl_d6_inst_PR_lag[itr,j]+precip_44_ctrl_d7_inst_PR_lag[itr,j]+precip_44_ctrl_d8_inst_PR_lag[itr,j]+precip_44_ctrl_d9_inst_PR_lag[itr,j]+precip_44_ctrl_d10_inst_PR_lag[itr,j])/9.

          proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag1[itr,j]+precip_44_ctrl_d3_inst_PRR_lag1[itr,j]+precip_44_ctrl_d4_inst_PRR_lag1[itr,j]+precip_44_ctrl_d5_inst_PRR_lag1[itr,j]+precip_44_ctrl_d6_inst_PRR_lag1[itr,j]+precip_44_ctrl_d7_inst_PRR_lag1[itr,j]+precip_44_ctrl_d8_inst_PRR_lag1[itr,j]+precip_44_ctrl_d9_inst_PRR_lag1[itr,j]+precip_44_ctrl_d10_inst_PRR_lag1[itr,j])/9.

          proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag2[itr,j]+precip_44_ctrl_d3_inst_PRR_lag2[itr,j]+precip_44_ctrl_d4_inst_PRR_lag2[itr,j]+precip_44_ctrl_d5_inst_PRR_lag2[itr,j]+precip_44_ctrl_d6_inst_PRR_lag2[itr,j]+precip_44_ctrl_d7_inst_PRR_lag2[itr,j]+precip_44_ctrl_d8_inst_PRR_lag2[itr,j]+precip_44_ctrl_d9_inst_PRR_lag2[itr,j]+precip_44_ctrl_d10_inst_PRR_lag2[itr,j])/9.

          proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag3[itr,j]+precip_44_ctrl_d3_inst_PRR_lag3[itr,j]+precip_44_ctrl_d4_inst_PRR_lag3[itr,j]+precip_44_ctrl_d5_inst_PRR_lag3[itr,j]+precip_44_ctrl_d6_inst_PRR_lag3[itr,j]+precip_44_ctrl_d7_inst_PRR_lag3[itr,j]+precip_44_ctrl_d8_inst_PRR_lag3[itr,j]+precip_44_ctrl_d9_inst_PRR_lag3[itr,j]+precip_44_ctrl_d10_inst_PRR_lag3[itr,j])/9.

          proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag4[itr,j]+precip_44_ctrl_d3_inst_PRR_lag4[itr,j]+precip_44_ctrl_d4_inst_PRR_lag4[itr,j]+precip_44_ctrl_d5_inst_PRR_lag4[itr,j]+precip_44_ctrl_d6_inst_PRR_lag4[itr,j]+precip_44_ctrl_d7_inst_PRR_lag4[itr,j]+precip_44_ctrl_d8_inst_PRR_lag4[itr,j]+precip_44_ctrl_d9_inst_PRR_lag4[itr,j]+precip_44_ctrl_d10_inst_PRR_lag4[itr,j])/9.

          proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag5[itr,j]+precip_44_ctrl_d3_inst_PRR_lag5[itr,j]+precip_44_ctrl_d4_inst_PRR_lag5[itr,j]+precip_44_ctrl_d5_inst_PRR_lag5[itr,j]+precip_44_ctrl_d6_inst_PRR_lag5[itr,j]+precip_44_ctrl_d7_inst_PRR_lag5[itr,j]+precip_44_ctrl_d8_inst_PRR_lag5[itr,j]+precip_44_ctrl_d9_inst_PRR_lag5[itr,j]+precip_44_ctrl_d10_inst_PRR_lag5[itr,j])/9.

          proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag6[itr,j]+precip_44_ctrl_d3_inst_PRR_lag6[itr,j]+precip_44_ctrl_d4_inst_PRR_lag6[itr,j]+precip_44_ctrl_d5_inst_PRR_lag6[itr,j]+precip_44_ctrl_d6_inst_PRR_lag6[itr,j]+precip_44_ctrl_d7_inst_PRR_lag6[itr,j]+precip_44_ctrl_d8_inst_PRR_lag6[itr,j]+precip_44_ctrl_d9_inst_PRR_lag6[itr,j]+precip_44_ctrl_d10_inst_PRR_lag6[itr,j])/9.

          proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag7[itr,j]+precip_44_ctrl_d3_inst_PRR_lag7[itr,j]+precip_44_ctrl_d4_inst_PRR_lag7[itr,j]+precip_44_ctrl_d5_inst_PRR_lag7[itr,j]+precip_44_ctrl_d6_inst_PRR_lag7[itr,j]+precip_44_ctrl_d7_inst_PRR_lag7[itr,j]+precip_44_ctrl_d8_inst_PRR_lag7[itr,j]+precip_44_ctrl_d9_inst_PRR_lag7[itr,j]+precip_44_ctrl_d10_inst_PRR_lag7[itr,j])/9.

          proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag8[itr,j]+precip_44_ctrl_d3_inst_PRR_lag8[itr,j]+precip_44_ctrl_d4_inst_PRR_lag8[itr,j]+precip_44_ctrl_d5_inst_PRR_lag8[itr,j]+precip_44_ctrl_d6_inst_PRR_lag8[itr,j]+precip_44_ctrl_d7_inst_PRR_lag8[itr,j]+precip_44_ctrl_d8_inst_PRR_lag8[itr,j]+precip_44_ctrl_d9_inst_PRR_lag8[itr,j]+precip_44_ctrl_d10_inst_PRR_lag8[itr,j])/9.

          proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag9[itr,j]+precip_44_ctrl_d3_inst_PRR_lag9[itr,j]+precip_44_ctrl_d4_inst_PRR_lag9[itr,j]+precip_44_ctrl_d5_inst_PRR_lag9[itr,j]+precip_44_ctrl_d6_inst_PRR_lag9[itr,j]+precip_44_ctrl_d7_inst_PRR_lag9[itr,j]+precip_44_ctrl_d8_inst_PRR_lag9[itr,j]+precip_44_ctrl_d9_inst_PRR_lag9[itr,j]+precip_44_ctrl_d10_inst_PRR_lag9[itr,j])/9.

          proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag10[itr,j]+precip_44_ctrl_d3_inst_PRR_lag10[itr,j]+precip_44_ctrl_d4_inst_PRR_lag10[itr,j]+precip_44_ctrl_d5_inst_PRR_lag10[itr,j]+precip_44_ctrl_d6_inst_PRR_lag10[itr,j]+precip_44_ctrl_d7_inst_PRR_lag10[itr,j]+precip_44_ctrl_d8_inst_PRR_lag10[itr,j]+precip_44_ctrl_d9_inst_PRR_lag10[itr,j]+precip_44_ctrl_d10_inst_PRR_lag10[itr,j])/9.

          proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag11[itr,j]+precip_44_ctrl_d3_inst_PRR_lag11[itr,j]+precip_44_ctrl_d4_inst_PRR_lag11[itr,j]+precip_44_ctrl_d5_inst_PRR_lag11[itr,j]+precip_44_ctrl_d6_inst_PRR_lag11[itr,j]+precip_44_ctrl_d7_inst_PRR_lag11[itr,j]+precip_44_ctrl_d8_inst_PRR_lag11[itr,j]+precip_44_ctrl_d9_inst_PRR_lag11[itr,j]+precip_44_ctrl_d10_inst_PRR_lag11[itr,j])/9.

          proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag12[itr,j]+precip_44_ctrl_d3_inst_PRR_lag12[itr,j]+precip_44_ctrl_d4_inst_PRR_lag12[itr,j]+precip_44_ctrl_d5_inst_PRR_lag12[itr,j]+precip_44_ctrl_d6_inst_PRR_lag12[itr,j]+precip_44_ctrl_d7_inst_PRR_lag12[itr,j]+precip_44_ctrl_d8_inst_PRR_lag12[itr,j]+precip_44_ctrl_d9_inst_PRR_lag12[itr,j]+precip_44_ctrl_d10_inst_PRR_lag12[itr,j])/9.

          proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag13[itr,j]+precip_44_ctrl_d3_inst_PRR_lag13[itr,j]+precip_44_ctrl_d4_inst_PRR_lag13[itr,j]+precip_44_ctrl_d5_inst_PRR_lag13[itr,j]+precip_44_ctrl_d6_inst_PRR_lag13[itr,j]+precip_44_ctrl_d7_inst_PRR_lag13[itr,j]+precip_44_ctrl_d8_inst_PRR_lag13[itr,j]+precip_44_ctrl_d9_inst_PRR_lag13[itr,j]+precip_44_ctrl_d10_inst_PRR_lag13[itr,j])/9.

          proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag14[itr,j]+precip_44_ctrl_d3_inst_PRR_lag14[itr,j]+precip_44_ctrl_d4_inst_PRR_lag14[itr,j]+precip_44_ctrl_d5_inst_PRR_lag14[itr,j]+precip_44_ctrl_d6_inst_PRR_lag14[itr,j]+precip_44_ctrl_d7_inst_PRR_lag14[itr,j]+precip_44_ctrl_d8_inst_PRR_lag14[itr,j]+precip_44_ctrl_d9_inst_PRR_lag14[itr,j]+precip_44_ctrl_d10_inst_PRR_lag14[itr,j])/9.

          proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag15[itr,j]+precip_44_ctrl_d3_inst_PRR_lag15[itr,j]+precip_44_ctrl_d4_inst_PRR_lag15[itr,j]+precip_44_ctrl_d5_inst_PRR_lag15[itr,j]+precip_44_ctrl_d6_inst_PRR_lag15[itr,j]+precip_44_ctrl_d7_inst_PRR_lag15[itr,j]+precip_44_ctrl_d8_inst_PRR_lag15[itr,j]+precip_44_ctrl_d9_inst_PRR_lag15[itr,j]+precip_44_ctrl_d10_inst_PRR_lag15[itr,j])/9.

          proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag16[itr,j]+precip_44_ctrl_d3_inst_PRR_lag16[itr,j]+precip_44_ctrl_d4_inst_PRR_lag16[itr,j]+precip_44_ctrl_d5_inst_PRR_lag16[itr,j]+precip_44_ctrl_d6_inst_PRR_lag16[itr,j]+precip_44_ctrl_d7_inst_PRR_lag16[itr,j]+precip_44_ctrl_d8_inst_PRR_lag16[itr,j]+precip_44_ctrl_d9_inst_PRR_lag16[itr,j]+precip_44_ctrl_d10_inst_PRR_lag16[itr,j])/9.

          proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag17[itr,j]+precip_44_ctrl_d3_inst_PRR_lag17[itr,j]+precip_44_ctrl_d4_inst_PRR_lag17[itr,j]+precip_44_ctrl_d5_inst_PRR_lag17[itr,j]+precip_44_ctrl_d6_inst_PRR_lag17[itr,j]+precip_44_ctrl_d7_inst_PRR_lag17[itr,j]+precip_44_ctrl_d8_inst_PRR_lag17[itr,j]+precip_44_ctrl_d9_inst_PRR_lag17[itr,j]+precip_44_ctrl_d10_inst_PRR_lag17[itr,j])/9.

          proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag18[itr,j]+precip_44_ctrl_d3_inst_PRR_lag18[itr,j]+precip_44_ctrl_d4_inst_PRR_lag18[itr,j]+precip_44_ctrl_d5_inst_PRR_lag18[itr,j]+precip_44_ctrl_d6_inst_PRR_lag18[itr,j]+precip_44_ctrl_d7_inst_PRR_lag18[itr,j]+precip_44_ctrl_d8_inst_PRR_lag18[itr,j]+precip_44_ctrl_d9_inst_PRR_lag18[itr,j]+precip_44_ctrl_d10_inst_PRR_lag18[itr,j])/9.

          proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag19[itr,j]+precip_44_ctrl_d3_inst_PRR_lag19[itr,j]+precip_44_ctrl_d4_inst_PRR_lag19[itr,j]+precip_44_ctrl_d5_inst_PRR_lag19[itr,j]+precip_44_ctrl_d6_inst_PRR_lag19[itr,j]+precip_44_ctrl_d7_inst_PRR_lag19[itr,j]+precip_44_ctrl_d8_inst_PRR_lag19[itr,j]+precip_44_ctrl_d9_inst_PRR_lag19[itr,j]+precip_44_ctrl_d10_inst_PRR_lag19[itr,j])/9.

          proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag20[itr,j]+precip_44_ctrl_d3_inst_PRR_lag20[itr,j]+precip_44_ctrl_d4_inst_PRR_lag20[itr,j]+precip_44_ctrl_d5_inst_PRR_lag20[itr,j]+precip_44_ctrl_d6_inst_PRR_lag20[itr,j]+precip_44_ctrl_d7_inst_PRR_lag20[itr,j]+precip_44_ctrl_d8_inst_PRR_lag20[itr,j]+precip_44_ctrl_d9_inst_PRR_lag20[itr,j]+precip_44_ctrl_d10_inst_PRR_lag20[itr,j])/9.

          proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag21[itr,j]+precip_44_ctrl_d3_inst_PRR_lag21[itr,j]+precip_44_ctrl_d4_inst_PRR_lag21[itr,j]+precip_44_ctrl_d5_inst_PRR_lag21[itr,j]+precip_44_ctrl_d6_inst_PRR_lag21[itr,j]+precip_44_ctrl_d7_inst_PRR_lag21[itr,j]+precip_44_ctrl_d8_inst_PRR_lag21[itr,j]+precip_44_ctrl_d9_inst_PRR_lag21[itr,j]+precip_44_ctrl_d10_inst_PRR_lag21[itr,j])/9.

          proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag22[itr,j]+precip_44_ctrl_d3_inst_PRR_lag22[itr,j]+precip_44_ctrl_d4_inst_PRR_lag22[itr,j]+precip_44_ctrl_d5_inst_PRR_lag22[itr,j]+precip_44_ctrl_d6_inst_PRR_lag22[itr,j]+precip_44_ctrl_d7_inst_PRR_lag22[itr,j]+precip_44_ctrl_d8_inst_PRR_lag22[itr,j]+precip_44_ctrl_d9_inst_PRR_lag22[itr,j]+precip_44_ctrl_d10_inst_PRR_lag22[itr,j])/9.

          proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag23[itr,j]+precip_44_ctrl_d3_inst_PRR_lag23[itr,j]+precip_44_ctrl_d4_inst_PRR_lag23[itr,j]+precip_44_ctrl_d5_inst_PRR_lag23[itr,j]+precip_44_ctrl_d6_inst_PRR_lag23[itr,j]+precip_44_ctrl_d7_inst_PRR_lag23[itr,j]+precip_44_ctrl_d8_inst_PRR_lag23[itr,j]+precip_44_ctrl_d9_inst_PRR_lag23[itr,j]+precip_44_ctrl_d10_inst_PRR_lag23[itr,j])/9.

          proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag24[itr,j]+precip_44_ctrl_d3_inst_PRR_lag24[itr,j]+precip_44_ctrl_d4_inst_PRR_lag24[itr,j]+precip_44_ctrl_d5_inst_PRR_lag24[itr,j]+precip_44_ctrl_d6_inst_PRR_lag24[itr,j]+precip_44_ctrl_d7_inst_PRR_lag24[itr,j]+precip_44_ctrl_d8_inst_PRR_lag24[itr,j]+precip_44_ctrl_d9_inst_PRR_lag24[itr,j]+precip_44_ctrl_d10_inst_PRR_lag24[itr,j])/9.

          proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag25[itr,j]+precip_44_ctrl_d3_inst_PRR_lag25[itr,j]+precip_44_ctrl_d4_inst_PRR_lag25[itr,j]+precip_44_ctrl_d5_inst_PRR_lag25[itr,j]+precip_44_ctrl_d6_inst_PRR_lag25[itr,j]+precip_44_ctrl_d7_inst_PRR_lag25[itr,j]+precip_44_ctrl_d8_inst_PRR_lag25[itr,j]+precip_44_ctrl_d9_inst_PRR_lag25[itr,j]+precip_44_ctrl_d10_inst_PRR_lag25[itr,j])/9.

          proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag26[itr,j]+precip_44_ctrl_d3_inst_PRR_lag26[itr,j]+precip_44_ctrl_d4_inst_PRR_lag26[itr,j]+precip_44_ctrl_d5_inst_PRR_lag26[itr,j]+precip_44_ctrl_d6_inst_PRR_lag26[itr,j]+precip_44_ctrl_d7_inst_PRR_lag26[itr,j]+precip_44_ctrl_d8_inst_PRR_lag26[itr,j]+precip_44_ctrl_d9_inst_PRR_lag26[itr,j]+precip_44_ctrl_d10_inst_PRR_lag26[itr,j])/9.

          proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag27[itr,j]+precip_44_ctrl_d3_inst_PRR_lag27[itr,j]+precip_44_ctrl_d4_inst_PRR_lag27[itr,j]+precip_44_ctrl_d5_inst_PRR_lag27[itr,j]+precip_44_ctrl_d6_inst_PRR_lag27[itr,j]+precip_44_ctrl_d7_inst_PRR_lag27[itr,j]+precip_44_ctrl_d8_inst_PRR_lag27[itr,j]+precip_44_ctrl_d9_inst_PRR_lag27[itr,j]+precip_44_ctrl_d10_inst_PRR_lag27[itr,j])/9.

          proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag28[itr,j]+precip_44_ctrl_d3_inst_PRR_lag28[itr,j]+precip_44_ctrl_d4_inst_PRR_lag28[itr,j]+precip_44_ctrl_d5_inst_PRR_lag28[itr,j]+precip_44_ctrl_d6_inst_PRR_lag28[itr,j]+precip_44_ctrl_d7_inst_PRR_lag28[itr,j]+precip_44_ctrl_d8_inst_PRR_lag28[itr,j]+precip_44_ctrl_d9_inst_PRR_lag28[itr,j]+precip_44_ctrl_d10_inst_PRR_lag28[itr,j])/9.

          proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag29[itr,j]+precip_44_ctrl_d3_inst_PRR_lag29[itr,j]+precip_44_ctrl_d4_inst_PRR_lag29[itr,j]+precip_44_ctrl_d5_inst_PRR_lag29[itr,j]+precip_44_ctrl_d6_inst_PRR_lag29[itr,j]+precip_44_ctrl_d7_inst_PRR_lag29[itr,j]+precip_44_ctrl_d8_inst_PRR_lag29[itr,j]+precip_44_ctrl_d9_inst_PRR_lag29[itr,j]+precip_44_ctrl_d10_inst_PRR_lag29[itr,j])/9.

          proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag30[itr,j]+precip_44_ctrl_d3_inst_PRR_lag30[itr,j]+precip_44_ctrl_d4_inst_PRR_lag30[itr,j]+precip_44_ctrl_d5_inst_PRR_lag30[itr,j]+precip_44_ctrl_d6_inst_PRR_lag30[itr,j]+precip_44_ctrl_d7_inst_PRR_lag30[itr,j]+precip_44_ctrl_d8_inst_PRR_lag30[itr,j]+precip_44_ctrl_d9_inst_PRR_lag30[itr,j]+precip_44_ctrl_d10_inst_PRR_lag30[itr,j])/9.

          proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag31[itr,j]+precip_44_ctrl_d3_inst_PRR_lag31[itr,j]+precip_44_ctrl_d4_inst_PRR_lag31[itr,j]+precip_44_ctrl_d5_inst_PRR_lag31[itr,j]+precip_44_ctrl_d6_inst_PRR_lag31[itr,j]+precip_44_ctrl_d7_inst_PRR_lag31[itr,j]+precip_44_ctrl_d8_inst_PRR_lag31[itr,j]+precip_44_ctrl_d9_inst_PRR_lag31[itr,j]+precip_44_ctrl_d10_inst_PRR_lag31[itr,j])/9.

          proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_lag32[itr,j]+precip_44_ctrl_d3_inst_PRR_lag32[itr,j]+precip_44_ctrl_d4_inst_PRR_lag32[itr,j]+precip_44_ctrl_d5_inst_PRR_lag32[itr,j]+precip_44_ctrl_d6_inst_PRR_lag32[itr,j]+precip_44_ctrl_d7_inst_PRR_lag32[itr,j]+precip_44_ctrl_d8_inst_PRR_lag32[itr,j]+precip_44_ctrl_d9_inst_PRR_lag32[itr,j]+precip_44_ctrl_d10_inst_PRR_lag32[itr,j])/9.

          proRNR_lag1_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRNR_lag1[itr,j]+precip_44_ctrl_d3_inst_PRNR_lag1[itr,j]+precip_44_ctrl_d4_inst_PRNR_lag1[itr,j]+precip_44_ctrl_d5_inst_PRNR_lag1[itr,j]+precip_44_ctrl_d6_inst_PRNR_lag1[itr,j]+precip_44_ctrl_d7_inst_PRNR_lag1[itr,j]+precip_44_ctrl_d8_inst_PRNR_lag1[itr,j]+precip_44_ctrl_d9_inst_PRNR_lag1[itr,j]+precip_44_ctrl_d10_inst_PRNR_lag1[itr,j])/9.

          proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag1[itr,j]+precip_44_ctrl_d3_inst_PR2_lag1[itr,j]+precip_44_ctrl_d4_inst_PR2_lag1[itr,j]+precip_44_ctrl_d5_inst_PR2_lag1[itr,j]+precip_44_ctrl_d6_inst_PR2_lag1[itr,j]+precip_44_ctrl_d7_inst_PR2_lag1[itr,j]+precip_44_ctrl_d8_inst_PR2_lag1[itr,j]+precip_44_ctrl_d9_inst_PR2_lag1[itr,j]+precip_44_ctrl_d10_inst_PR2_lag1[itr,j])/9.

          proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag2[itr,j]+precip_44_ctrl_d3_inst_PR2_lag2[itr,j]+precip_44_ctrl_d4_inst_PR2_lag2[itr,j]+precip_44_ctrl_d5_inst_PR2_lag2[itr,j]+precip_44_ctrl_d6_inst_PR2_lag2[itr,j]+precip_44_ctrl_d7_inst_PR2_lag2[itr,j]+precip_44_ctrl_d8_inst_PR2_lag2[itr,j]+precip_44_ctrl_d9_inst_PR2_lag2[itr,j]+precip_44_ctrl_d10_inst_PR2_lag2[itr,j])/9.

          proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag3[itr,j]+precip_44_ctrl_d3_inst_PR2_lag3[itr,j]+precip_44_ctrl_d4_inst_PR2_lag3[itr,j]+precip_44_ctrl_d5_inst_PR2_lag3[itr,j]+precip_44_ctrl_d6_inst_PR2_lag3[itr,j]+precip_44_ctrl_d7_inst_PR2_lag3[itr,j]+precip_44_ctrl_d8_inst_PR2_lag3[itr,j]+precip_44_ctrl_d9_inst_PR2_lag3[itr,j]+precip_44_ctrl_d10_inst_PR2_lag3[itr,j])/9.

          proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag4[itr,j]+precip_44_ctrl_d3_inst_PR2_lag4[itr,j]+precip_44_ctrl_d4_inst_PR2_lag4[itr,j]+precip_44_ctrl_d5_inst_PR2_lag4[itr,j]+precip_44_ctrl_d6_inst_PR2_lag4[itr,j]+precip_44_ctrl_d7_inst_PR2_lag4[itr,j]+precip_44_ctrl_d8_inst_PR2_lag4[itr,j]+precip_44_ctrl_d9_inst_PR2_lag4[itr,j]+precip_44_ctrl_d10_inst_PR2_lag4[itr,j])/9.

          proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag5[itr,j]+precip_44_ctrl_d3_inst_PR2_lag5[itr,j]+precip_44_ctrl_d4_inst_PR2_lag5[itr,j]+precip_44_ctrl_d5_inst_PR2_lag5[itr,j]+precip_44_ctrl_d6_inst_PR2_lag5[itr,j]+precip_44_ctrl_d7_inst_PR2_lag5[itr,j]+precip_44_ctrl_d8_inst_PR2_lag5[itr,j]+precip_44_ctrl_d9_inst_PR2_lag5[itr,j]+precip_44_ctrl_d10_inst_PR2_lag5[itr,j])/9.

          proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag6[itr,j]+precip_44_ctrl_d3_inst_PR2_lag6[itr,j]+precip_44_ctrl_d4_inst_PR2_lag6[itr,j]+precip_44_ctrl_d5_inst_PR2_lag6[itr,j]+precip_44_ctrl_d6_inst_PR2_lag6[itr,j]+precip_44_ctrl_d7_inst_PR2_lag6[itr,j]+precip_44_ctrl_d8_inst_PR2_lag6[itr,j]+precip_44_ctrl_d9_inst_PR2_lag6[itr,j]+precip_44_ctrl_d10_inst_PR2_lag6[itr,j])/9.

          proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag7[itr,j]+precip_44_ctrl_d3_inst_PR2_lag7[itr,j]+precip_44_ctrl_d4_inst_PR2_lag7[itr,j]+precip_44_ctrl_d5_inst_PR2_lag7[itr,j]+precip_44_ctrl_d6_inst_PR2_lag7[itr,j]+precip_44_ctrl_d7_inst_PR2_lag7[itr,j]+precip_44_ctrl_d8_inst_PR2_lag7[itr,j]+precip_44_ctrl_d9_inst_PR2_lag7[itr,j]+precip_44_ctrl_d10_inst_PR2_lag7[itr,j])/9.

          proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag8[itr,j]+precip_44_ctrl_d3_inst_PR2_lag8[itr,j]+precip_44_ctrl_d4_inst_PR2_lag8[itr,j]+precip_44_ctrl_d5_inst_PR2_lag8[itr,j]+precip_44_ctrl_d6_inst_PR2_lag8[itr,j]+precip_44_ctrl_d7_inst_PR2_lag8[itr,j]+precip_44_ctrl_d8_inst_PR2_lag8[itr,j]+precip_44_ctrl_d9_inst_PR2_lag8[itr,j]+precip_44_ctrl_d10_inst_PR2_lag8[itr,j])/9.

          proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag9[itr,j]+precip_44_ctrl_d3_inst_PR2_lag9[itr,j]+precip_44_ctrl_d4_inst_PR2_lag9[itr,j]+precip_44_ctrl_d5_inst_PR2_lag9[itr,j]+precip_44_ctrl_d6_inst_PR2_lag9[itr,j]+precip_44_ctrl_d7_inst_PR2_lag9[itr,j]+precip_44_ctrl_d8_inst_PR2_lag9[itr,j]+precip_44_ctrl_d9_inst_PR2_lag9[itr,j]+precip_44_ctrl_d10_inst_PR2_lag9[itr,j])/9.

          proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag10[itr,j]+precip_44_ctrl_d3_inst_PR2_lag10[itr,j]+precip_44_ctrl_d4_inst_PR2_lag10[itr,j]+precip_44_ctrl_d5_inst_PR2_lag10[itr,j]+precip_44_ctrl_d6_inst_PR2_lag10[itr,j]+precip_44_ctrl_d7_inst_PR2_lag10[itr,j]+precip_44_ctrl_d8_inst_PR2_lag10[itr,j]+precip_44_ctrl_d9_inst_PR2_lag10[itr,j]+precip_44_ctrl_d10_inst_PR2_lag10[itr,j])/9.

          proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag11[itr,j]+precip_44_ctrl_d3_inst_PR2_lag11[itr,j]+precip_44_ctrl_d4_inst_PR2_lag11[itr,j]+precip_44_ctrl_d5_inst_PR2_lag11[itr,j]+precip_44_ctrl_d6_inst_PR2_lag11[itr,j]+precip_44_ctrl_d7_inst_PR2_lag11[itr,j]+precip_44_ctrl_d8_inst_PR2_lag11[itr,j]+precip_44_ctrl_d9_inst_PR2_lag11[itr,j]+precip_44_ctrl_d10_inst_PR2_lag11[itr,j])/9.

          proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag12[itr,j]+precip_44_ctrl_d3_inst_PR2_lag12[itr,j]+precip_44_ctrl_d4_inst_PR2_lag12[itr,j]+precip_44_ctrl_d5_inst_PR2_lag12[itr,j]+precip_44_ctrl_d6_inst_PR2_lag12[itr,j]+precip_44_ctrl_d7_inst_PR2_lag12[itr,j]+precip_44_ctrl_d8_inst_PR2_lag12[itr,j]+precip_44_ctrl_d9_inst_PR2_lag12[itr,j]+precip_44_ctrl_d10_inst_PR2_lag12[itr,j])/9.

          proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag13[itr,j]+precip_44_ctrl_d3_inst_PR2_lag13[itr,j]+precip_44_ctrl_d4_inst_PR2_lag13[itr,j]+precip_44_ctrl_d5_inst_PR2_lag13[itr,j]+precip_44_ctrl_d6_inst_PR2_lag13[itr,j]+precip_44_ctrl_d7_inst_PR2_lag13[itr,j]+precip_44_ctrl_d8_inst_PR2_lag13[itr,j]+precip_44_ctrl_d9_inst_PR2_lag13[itr,j]+precip_44_ctrl_d10_inst_PR2_lag13[itr,j])/9.

          proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag14[itr,j]+precip_44_ctrl_d3_inst_PR2_lag14[itr,j]+precip_44_ctrl_d4_inst_PR2_lag14[itr,j]+precip_44_ctrl_d5_inst_PR2_lag14[itr,j]+precip_44_ctrl_d6_inst_PR2_lag14[itr,j]+precip_44_ctrl_d7_inst_PR2_lag14[itr,j]+precip_44_ctrl_d8_inst_PR2_lag14[itr,j]+precip_44_ctrl_d9_inst_PR2_lag14[itr,j]+precip_44_ctrl_d10_inst_PR2_lag14[itr,j])/9.

          proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag15[itr,j]+precip_44_ctrl_d3_inst_PR2_lag15[itr,j]+precip_44_ctrl_d4_inst_PR2_lag15[itr,j]+precip_44_ctrl_d5_inst_PR2_lag15[itr,j]+precip_44_ctrl_d6_inst_PR2_lag15[itr,j]+precip_44_ctrl_d7_inst_PR2_lag15[itr,j]+precip_44_ctrl_d8_inst_PR2_lag15[itr,j]+precip_44_ctrl_d9_inst_PR2_lag15[itr,j]+precip_44_ctrl_d10_inst_PR2_lag15[itr,j])/9.

          proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag16[itr,j]+precip_44_ctrl_d3_inst_PR2_lag16[itr,j]+precip_44_ctrl_d4_inst_PR2_lag16[itr,j]+precip_44_ctrl_d5_inst_PR2_lag16[itr,j]+precip_44_ctrl_d6_inst_PR2_lag16[itr,j]+precip_44_ctrl_d7_inst_PR2_lag16[itr,j]+precip_44_ctrl_d8_inst_PR2_lag16[itr,j]+precip_44_ctrl_d9_inst_PR2_lag16[itr,j]+precip_44_ctrl_d10_inst_PR2_lag16[itr,j])/9.

          proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag17[itr,j]+precip_44_ctrl_d3_inst_PR2_lag17[itr,j]+precip_44_ctrl_d4_inst_PR2_lag17[itr,j]+precip_44_ctrl_d5_inst_PR2_lag17[itr,j]+precip_44_ctrl_d6_inst_PR2_lag17[itr,j]+precip_44_ctrl_d7_inst_PR2_lag17[itr,j]+precip_44_ctrl_d8_inst_PR2_lag17[itr,j]+precip_44_ctrl_d9_inst_PR2_lag17[itr,j]+precip_44_ctrl_d10_inst_PR2_lag17[itr,j])/9.

          proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag18[itr,j]+precip_44_ctrl_d3_inst_PR2_lag18[itr,j]+precip_44_ctrl_d4_inst_PR2_lag18[itr,j]+precip_44_ctrl_d5_inst_PR2_lag18[itr,j]+precip_44_ctrl_d6_inst_PR2_lag18[itr,j]+precip_44_ctrl_d7_inst_PR2_lag18[itr,j]+precip_44_ctrl_d8_inst_PR2_lag18[itr,j]+precip_44_ctrl_d9_inst_PR2_lag18[itr,j]+precip_44_ctrl_d10_inst_PR2_lag18[itr,j])/9.

          proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag19[itr,j]+precip_44_ctrl_d3_inst_PR2_lag19[itr,j]+precip_44_ctrl_d4_inst_PR2_lag19[itr,j]+precip_44_ctrl_d5_inst_PR2_lag19[itr,j]+precip_44_ctrl_d6_inst_PR2_lag19[itr,j]+precip_44_ctrl_d7_inst_PR2_lag19[itr,j]+precip_44_ctrl_d8_inst_PR2_lag19[itr,j]+precip_44_ctrl_d9_inst_PR2_lag19[itr,j]+precip_44_ctrl_d10_inst_PR2_lag19[itr,j])/9.

          proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag20[itr,j]+precip_44_ctrl_d3_inst_PR2_lag20[itr,j]+precip_44_ctrl_d4_inst_PR2_lag20[itr,j]+precip_44_ctrl_d5_inst_PR2_lag20[itr,j]+precip_44_ctrl_d6_inst_PR2_lag20[itr,j]+precip_44_ctrl_d7_inst_PR2_lag20[itr,j]+precip_44_ctrl_d8_inst_PR2_lag20[itr,j]+precip_44_ctrl_d9_inst_PR2_lag20[itr,j]+precip_44_ctrl_d10_inst_PR2_lag20[itr,j])/9.

          proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag21[itr,j]+precip_44_ctrl_d3_inst_PR2_lag21[itr,j]+precip_44_ctrl_d4_inst_PR2_lag21[itr,j]+precip_44_ctrl_d5_inst_PR2_lag21[itr,j]+precip_44_ctrl_d6_inst_PR2_lag21[itr,j]+precip_44_ctrl_d7_inst_PR2_lag21[itr,j]+precip_44_ctrl_d8_inst_PR2_lag21[itr,j]+precip_44_ctrl_d9_inst_PR2_lag21[itr,j]+precip_44_ctrl_d10_inst_PR2_lag21[itr,j])/9.

          proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag22[itr,j]+precip_44_ctrl_d3_inst_PR2_lag22[itr,j]+precip_44_ctrl_d4_inst_PR2_lag22[itr,j]+precip_44_ctrl_d5_inst_PR2_lag22[itr,j]+precip_44_ctrl_d6_inst_PR2_lag22[itr,j]+precip_44_ctrl_d7_inst_PR2_lag22[itr,j]+precip_44_ctrl_d8_inst_PR2_lag22[itr,j]+precip_44_ctrl_d9_inst_PR2_lag22[itr,j]+precip_44_ctrl_d10_inst_PR2_lag22[itr,j])/9.

          proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag23[itr,j]+precip_44_ctrl_d3_inst_PR2_lag23[itr,j]+precip_44_ctrl_d4_inst_PR2_lag23[itr,j]+precip_44_ctrl_d5_inst_PR2_lag23[itr,j]+precip_44_ctrl_d6_inst_PR2_lag23[itr,j]+precip_44_ctrl_d7_inst_PR2_lag23[itr,j]+precip_44_ctrl_d8_inst_PR2_lag23[itr,j]+precip_44_ctrl_d9_inst_PR2_lag23[itr,j]+precip_44_ctrl_d10_inst_PR2_lag23[itr,j])/9.

          proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag24[itr,j]+precip_44_ctrl_d3_inst_PR2_lag24[itr,j]+precip_44_ctrl_d4_inst_PR2_lag24[itr,j]+precip_44_ctrl_d5_inst_PR2_lag24[itr,j]+precip_44_ctrl_d6_inst_PR2_lag24[itr,j]+precip_44_ctrl_d7_inst_PR2_lag24[itr,j]+precip_44_ctrl_d8_inst_PR2_lag24[itr,j]+precip_44_ctrl_d9_inst_PR2_lag24[itr,j]+precip_44_ctrl_d10_inst_PR2_lag24[itr,j])/9.

          proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag25[itr,j]+precip_44_ctrl_d3_inst_PR2_lag25[itr,j]+precip_44_ctrl_d4_inst_PR2_lag25[itr,j]+precip_44_ctrl_d5_inst_PR2_lag25[itr,j]+precip_44_ctrl_d6_inst_PR2_lag25[itr,j]+precip_44_ctrl_d7_inst_PR2_lag25[itr,j]+precip_44_ctrl_d8_inst_PR2_lag25[itr,j]+precip_44_ctrl_d9_inst_PR2_lag25[itr,j]+precip_44_ctrl_d10_inst_PR2_lag25[itr,j])/9.

          proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag26[itr,j]+precip_44_ctrl_d3_inst_PR2_lag26[itr,j]+precip_44_ctrl_d4_inst_PR2_lag26[itr,j]+precip_44_ctrl_d5_inst_PR2_lag26[itr,j]+precip_44_ctrl_d6_inst_PR2_lag26[itr,j]+precip_44_ctrl_d7_inst_PR2_lag26[itr,j]+precip_44_ctrl_d8_inst_PR2_lag26[itr,j]+precip_44_ctrl_d9_inst_PR2_lag26[itr,j]+precip_44_ctrl_d10_inst_PR2_lag26[itr,j])/9.

          proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag27[itr,j]+precip_44_ctrl_d3_inst_PR2_lag27[itr,j]+precip_44_ctrl_d4_inst_PR2_lag27[itr,j]+precip_44_ctrl_d5_inst_PR2_lag27[itr,j]+precip_44_ctrl_d6_inst_PR2_lag27[itr,j]+precip_44_ctrl_d7_inst_PR2_lag27[itr,j]+precip_44_ctrl_d8_inst_PR2_lag27[itr,j]+precip_44_ctrl_d9_inst_PR2_lag27[itr,j]+precip_44_ctrl_d10_inst_PR2_lag27[itr,j])/9.

          proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag28[itr,j]+precip_44_ctrl_d3_inst_PR2_lag28[itr,j]+precip_44_ctrl_d4_inst_PR2_lag28[itr,j]+precip_44_ctrl_d5_inst_PR2_lag28[itr,j]+precip_44_ctrl_d6_inst_PR2_lag28[itr,j]+precip_44_ctrl_d7_inst_PR2_lag28[itr,j]+precip_44_ctrl_d8_inst_PR2_lag28[itr,j]+precip_44_ctrl_d9_inst_PR2_lag28[itr,j]+precip_44_ctrl_d10_inst_PR2_lag28[itr,j])/9.

          proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag29[itr,j]+precip_44_ctrl_d3_inst_PR2_lag29[itr,j]+precip_44_ctrl_d4_inst_PR2_lag29[itr,j]+precip_44_ctrl_d5_inst_PR2_lag29[itr,j]+precip_44_ctrl_d6_inst_PR2_lag29[itr,j]+precip_44_ctrl_d7_inst_PR2_lag29[itr,j]+precip_44_ctrl_d8_inst_PR2_lag29[itr,j]+precip_44_ctrl_d9_inst_PR2_lag29[itr,j]+precip_44_ctrl_d10_inst_PR2_lag29[itr,j])/9.

          proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag30[itr,j]+precip_44_ctrl_d3_inst_PR2_lag30[itr,j]+precip_44_ctrl_d4_inst_PR2_lag30[itr,j]+precip_44_ctrl_d5_inst_PR2_lag30[itr,j]+precip_44_ctrl_d6_inst_PR2_lag30[itr,j]+precip_44_ctrl_d7_inst_PR2_lag30[itr,j]+precip_44_ctrl_d8_inst_PR2_lag30[itr,j]+precip_44_ctrl_d9_inst_PR2_lag30[itr,j]+precip_44_ctrl_d10_inst_PR2_lag30[itr,j])/9.

          proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag31[itr,j]+precip_44_ctrl_d3_inst_PR2_lag31[itr,j]+precip_44_ctrl_d4_inst_PR2_lag31[itr,j]+precip_44_ctrl_d5_inst_PR2_lag31[itr,j]+precip_44_ctrl_d6_inst_PR2_lag31[itr,j]+precip_44_ctrl_d7_inst_PR2_lag31[itr,j]+precip_44_ctrl_d8_inst_PR2_lag31[itr,j]+precip_44_ctrl_d9_inst_PR2_lag31[itr,j]+precip_44_ctrl_d10_inst_PR2_lag31[itr,j])/9.

          proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PR2_lag32[itr,j]+precip_44_ctrl_d3_inst_PR2_lag32[itr,j]+precip_44_ctrl_d4_inst_PR2_lag32[itr,j]+precip_44_ctrl_d5_inst_PR2_lag32[itr,j]+precip_44_ctrl_d6_inst_PR2_lag32[itr,j]+precip_44_ctrl_d7_inst_PR2_lag32[itr,j]+precip_44_ctrl_d8_inst_PR2_lag32[itr,j]+precip_44_ctrl_d9_inst_PR2_lag32[itr,j]+precip_44_ctrl_d10_inst_PR2_lag32[itr,j])/9.


          proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag1[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag1[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag1[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag1[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag1[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag1[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag1[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag1[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag1[itr,j])/9.

          proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag2[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag2[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag2[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag2[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag2[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag2[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag2[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag2[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag2[itr,j])/9.
          proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag3[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag3[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag3[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag3[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag3[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag3[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag3[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag3[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag3[itr,j])/9.
          proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag4[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag4[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag4[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag4[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag4[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag4[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag4[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag4[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag4[itr,j])/9.
          proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag5[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag5[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag5[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag5[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag5[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag5[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag5[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag5[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag5[itr,j])/9.
          proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag6[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag6[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag6[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag6[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag6[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag6[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag6[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag6[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag6[itr,j])/9.
          proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag7[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag7[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag7[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag7[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag7[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag7[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag7[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag7[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag7[itr,j])/9.
          proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag8[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag8[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag8[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag8[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag8[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag8[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag8[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag8[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag8[itr,j])/9.
          proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag9[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag9[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag9[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag9[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag9[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag9[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag9[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag9[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag9[itr,j])/9.
          proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag10[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag10[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag10[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag10[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag10[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag10[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag10[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag10[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag10[itr,j])/9.
          proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag11[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag11[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag11[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag11[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag11[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag11[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag11[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag11[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag11[itr,j])/9.
          proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag12[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag12[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag12[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag12[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag12[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag12[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag12[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag12[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag12[itr,j])/9.
          proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag13[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag13[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag13[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag13[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag13[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag13[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag13[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag13[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag13[itr,j])/9.
          proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag14[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag14[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag14[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag14[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag14[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag14[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag14[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag14[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag14[itr,j])/9.
          proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag15[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag15[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag15[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag15[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag15[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag15[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag15[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag15[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag15[itr,j])/9.
          proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag16[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag16[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag16[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag16[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag16[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag16[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag16[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag16[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag16[itr,j])/9.
          proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag17[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag17[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag17[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag17[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag17[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag17[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag17[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag17[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag17[itr,j])/9.
          proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag18[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag18[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag18[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag18[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag18[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag18[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag18[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag18[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag18[itr,j])/9.
          proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag19[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag19[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag19[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag19[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag19[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag19[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag19[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag19[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag19[itr,j])/9.
          proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag20[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag20[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag20[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag20[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag20[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag20[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag20[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag20[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag20[itr,j])/9.
          proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag21[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag21[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag21[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag21[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag21[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag21[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag21[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag21[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag21[itr,j])/9.
          proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag22[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag22[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag22[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag22[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag22[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag22[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag22[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag22[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag22[itr,j])/9.
          proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag23[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag23[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag23[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag23[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag23[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag23[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag23[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag23[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag23[itr,j])/9.
          proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag24[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag24[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag24[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag24[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag24[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag24[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag24[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag24[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag24[itr,j])/9.
          proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag25[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag25[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag25[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag25[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag25[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag25[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag25[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag25[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag25[itr,j])/9.
          proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag26[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag26[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag26[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag26[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag26[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag26[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag26[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag26[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag26[itr,j])/9.
          proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag27[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag27[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag27[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag27[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag27[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag27[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag27[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag27[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag27[itr,j])/9.
          proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag28[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag28[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag28[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag28[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag28[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag28[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag28[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag28[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag28[itr,j])/9.
          proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag29[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag29[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag29[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag29[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag29[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag29[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag29[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag29[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag29[itr,j])/9.
          proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag30[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag30[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag30[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag30[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag30[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag30[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag30[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag30[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag30[itr,j])/9.
          proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag31[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag31[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag31[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag31[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag31[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag31[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag31[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag31[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag31[itr,j])/9.
          proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag32[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag32[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag32[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag32[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag32[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag32[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag32[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag32[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag32[itr,j])/9.
          proRR_proR2_lag33_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag33[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag33[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag33[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag33[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag33[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag33[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag33[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag33[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag33[itr,j])/9.
          proRR_proR2_lag34_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag34[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag34[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag34[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag34[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag34[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag34[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag34[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag34[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag34[itr,j])/9.
          proRR_proR2_lag35_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag35[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag35[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag35[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag35[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag35[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag35[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag35[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag35[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag35[itr,j])/9.
          proRR_proR2_lag36_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag36[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag36[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag36[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag36[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag36[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag36[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag36[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag36[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag36[itr,j])/9.
          proRR_proR2_lag37_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag37[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag37[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag37[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag37[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag37[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag37[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag37[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag37[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag37[itr,j])/9.
          proRR_proR2_lag38_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag38[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag38[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag38[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag38[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag38[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag38[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag38[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag38[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag38[itr,j])/9.
          proRR_proR2_lag39_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag39[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag39[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag39[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag39[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag39[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag39[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag39[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag39[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag39[itr,j])/9.
          proRR_proR2_lag40_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag40[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag40[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag40[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag40[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag40[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag40[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag40[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag40[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag40[itr,j])/9.
          proRR_proR2_lag41_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag41[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag41[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag41[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag41[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag41[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag41[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag41[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag41[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag41[itr,j])/9.
          proRR_proR2_lag42_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag42[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag42[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag42[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag42[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag42[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag42[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag42[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag42[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag42[itr,j])/9.
          proRR_proR2_lag43_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag43[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag43[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag43[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag43[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag43[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag43[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag43[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag43[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag43[itr,j])/9.
          proRR_proR2_lag44_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag44[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag44[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag44[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag44[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag44[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag44[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag44[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag44[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag44[itr,j])/9.
          proRR_proR2_lag45_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag45[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag45[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag45[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag45[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag45[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag45[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag45[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag45[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag45[itr,j])/9.
          proRR_proR2_lag46_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag46[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag46[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag46[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag46[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag46[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag46[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag46[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag46[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag46[itr,j])/9.
          proRR_proR2_lag47_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag47[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag47[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag47[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag47[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag47[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag47[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag47[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag47[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag47[itr,j])/9.
          proRR_proR2_lag48_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag48[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag48[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag48[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag48[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag48[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag48[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag48[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag48[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag48[itr,j])/9.
          proRR_proR2_lag49_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag49[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag49[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag49[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag49[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag49[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag49[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag49[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag49[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag49[itr,j])/9.
          proRR_proR2_lag50_ens_mean_inst_44_vs_AveRR_thr[itr,j]=(precip_44_ctrl_d2_inst_PRR_PR2_lag50[itr,j]+precip_44_ctrl_d3_inst_PRR_PR2_lag50[itr,j]+precip_44_ctrl_d4_inst_PRR_PR2_lag50[itr,j]+precip_44_ctrl_d5_inst_PRR_PR2_lag50[itr,j]+precip_44_ctrl_d6_inst_PRR_PR2_lag50[itr,j]+precip_44_ctrl_d7_inst_PRR_PR2_lag50[itr,j]+precip_44_ctrl_d8_inst_PRR_PR2_lag50[itr,j]+precip_44_ctrl_d9_inst_PRR_PR2_lag50[itr,j]+precip_44_ctrl_d10_inst_PRR_PR2_lag50[itr,j])/9.

  
nlag=55
lag=N.zeros((nlag))
zero_lag=N.zeros((nlag))
lag[0]=0
for j in N.arange(nlag-1):
     lag[j+1]=lag[j]-0.25
zero_lag[:]=0.0


proR2t5_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t6_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t7_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t8_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t9_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t10_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t11_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t12_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t13_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t14_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t15_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t16_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t17_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t18_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t19_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t20_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t21_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t22_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t23_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t24_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t25_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t26_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t27_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t28_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t29_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t30_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t31_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t32_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t33_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t34_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t35_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t36_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t37_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t38_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t39_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t40_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t41_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t42_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t43_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t44_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t45_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t46_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t47_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t48_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t49_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t50_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t51_lag_ens_mean_inst_44=N.zeros((nlag))
proR2t52_lag_ens_mean_inst_44=N.zeros((nlag))


proR2t5_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,4]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t6_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,5]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t7_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,6]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t8_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,7]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t9_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,8]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t10_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,9]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t11_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,10]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t12_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,11]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t13_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,12]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t14_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,13]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t15_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,14]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t16_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,15]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t17_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,16]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t18_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,17]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t19_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,18]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t20_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,19]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t21_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,20]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t22_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,21]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t23_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,22]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t24_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,23]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t25_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,24]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t26_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,25]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t27_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,26]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t28_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,27]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t29_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,28]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t30_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,29]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t31_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,30]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t32_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,31]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t33_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,32]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t34_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,33]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t35_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,34]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t36_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,35]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t37_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,36]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t38_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,37]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t39_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,38]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t40_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,39]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t41_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,10]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t42_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,41]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t43_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,42]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t44_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,43]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t45_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,44]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t46_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,45]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t47_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,46]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t48_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,47]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t49_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,48]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t50_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,49]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t51_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,50]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t52_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,51]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,51]

proR2t5_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,4]
proR2t5_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,4]

proR2t6_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,5]

proR2t6_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,5]
proR2t6_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,5]

proR2t7_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,6]
proR2t7_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,6]

proR2t8_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,7]
proR2t8_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,7]

proR2t9_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,8]
proR2t9_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,8]





proR2t10_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,9]

proR2t10_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,9]
proR2t10_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,9]


proR2t11_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,10]
proR2t11_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,10]


proR2t12_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,11]
proR2t12_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,11]


proR2t13_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,12]

proR2t14_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,13]

proR2t15_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,14]

proR2t16_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,15]

proR2t17_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,16]

proR2t18_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,17]

proR2t19_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,18]

proR2t20_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,19]

proR2t21_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,20]

proR2t22_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,21]

proR2t23_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,22]

proR2t24_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,23]

proR2t25_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,24]

proR2t26_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,25]

proR2t27_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,26]

proR2t28_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,27]

proR2t13_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,12]
proR2t13_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,12]


proR2t14_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,13]
proR2t14_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,13]

proR2t15_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,14]
proR2t15_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,14]


proR2t16_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,15]
proR2t16_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,15]

proR2t17_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,16]
proR2t17_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,16]

proR2t18_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,17]
proR2t18_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,17]


proR2t19_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,18]
proR2t19_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,18]

proR2t20_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,19]
proR2t20_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,19]


proR2t21_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,20]
proR2t21_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,20]

proR2t22_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,21]
proR2t22_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,21]

proR2t23_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,22]
proR2t23_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,22]

proR2t24_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,23]
proR2t24_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,23]

proR2t25_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,24]
proR2t25_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,24]


proR2t26_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,25]
proR2t26_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,25]

proR2t27_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,26]
proR2t27_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,26]


proR2t28_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,27]
proR2t28_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,27]


proR2t29_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[28]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,28]
proR2t29_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,28]

proR2t30_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,29]
proR2t30_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,29]

proR2t31_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,30]
proR2t31_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,30]






proR2t32_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,31]
proR2t32_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,31]


proR2t33_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,32]
proR2t33_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,32]


proR2t34_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,33]
proR2t34_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,33]







proR2t35_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,34]
proR2t35_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,34]



proR2t36_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,35]
proR2t36_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,35]

proR2t37_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,36]
proR2t37_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,36]


proR2t38_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,37]
proR2t38_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,37]



proR2t39_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,38]
proR2t39_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,38]

proR2t40_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,39]
proR2t40_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,39]




proR2t41_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,40]
proR2t41_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,40]

proR2t42_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,41]
proR2t42_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,41]

proR2t43_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,42]
proR2t43_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,42]


proR2t44_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,43]
proR2t44_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,43]



proR2t45_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,44]
proR2t45_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,44]


proR2t46_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,45]
proR2t46_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,45]



proR2t47_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,46]
proR2t47_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,46]


proR2t48_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,47]
proR2t48_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,47]

proR2t49_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,48]
proR2t49_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,48]


proR2t50_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,49]
proR2t50_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,49]


proR2t51_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,50]
proR2t51_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,50]


proR2t52_lag_ens_mean_inst_44[1]=proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[2]=proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[3]=proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[4]=proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[5]=proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[6]=proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[7]=proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[8]=proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[9]=proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[10]=proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[11]=proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[12]=proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[13]=proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[14]=proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[15]=proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[16]=proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[17]=proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[18]=proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[19]=proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[20]=proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[21]=proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[22]=proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[23]=proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[24]=proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[25]=proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[26]=proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[27]=proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[28]=proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[29]=proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[30]=proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[31]=proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,51]
proR2t52_lag_ens_mean_inst_44[32]=proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,51]

proRRt5_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt6_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt7_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt8_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt9_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt10_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt11_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt12_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt13_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt14_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt15_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt16_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt17_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt18_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt19_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt20_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt21_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt22_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt23_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt24_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt25_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt26_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt27_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt28_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt29_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt30_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt31_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt32_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt33_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt34_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt35_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt36_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt37_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt38_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt39_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt40_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt41_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt42_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt43_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt44_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt45_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt46_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt47_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt48_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt49_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt50_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt51_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt52_lag_ens_mean_inst_44=N.zeros((nlag))


proRRt5_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,4]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt6_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,5]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt7_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,6]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt8_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,7]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt9_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,8]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt10_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,9]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt11_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,10]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt12_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,11]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt13_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,12]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt14_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,13]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt15_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,14]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt16_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,15]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt17_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,16]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt18_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,17]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt19_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,18]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt20_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,19]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt21_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,20]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt22_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,21]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt23_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,22]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt24_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,23]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt25_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,24]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt26_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,25]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt27_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,26]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt28_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,27]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt29_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,28]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt30_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,29]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt31_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,30]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt32_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,31]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt33_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,32]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt34_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,33]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt35_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,34]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt36_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,35]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt37_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,36]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt38_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,37]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt39_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,38]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt40_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,39]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt41_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,10]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt42_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,41]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt43_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,42]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt44_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,43]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt45_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,44]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt46_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,45]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt47_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,46]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt48_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,47]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt49_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,48]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt50_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,49]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt51_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,50]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt52_lag_ens_mean_inst_44[0]=proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,51]*proR_lag_ens_mean_inst_44_vs_AveRR_thr[0,51]

proRRt5_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,4]

proRRt6_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,5]

proRRt6_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,5]

proRRt7_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,6]

proRRt8_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,7]

proRRt9_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,8]





proRRt10_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,9]

proRRt10_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,9]


proRRt11_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,10]


proRRt12_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,11]


proRRt13_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,12]

proRRt14_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,13]

proRRt15_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,14]

proRRt16_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,15]

proRRt17_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,16]

proRRt18_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,17]

proRRt19_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,18]

proRRt20_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,19]

proRRt21_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,20]

proRRt22_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,21]

proRRt23_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,22]

proRRt24_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,23]

proRRt25_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,24]

proRRt26_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,25]

proRRt27_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,26]

proRRt28_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,27]

proRRt13_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,12]


proRRt14_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,13]

proRRt15_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,14]


proRRt16_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,15]

proRRt17_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,16]

proRRt18_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,17]


proRRt19_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,18]

proRRt20_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,19]


proRRt21_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,20]

proRRt22_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,21]

proRRt23_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,22]

proRRt24_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,23]

proRRt25_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,24]


proRRt26_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,25]

proRRt27_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,26]


proRRt28_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,27]


proRRt29_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[28]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,28]

proRRt30_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,29]

proRRt31_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,30]






proRRt32_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,31]


proRRt33_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,32]







proRRt34_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,33]







proRRt35_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,34]



proRRt36_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,35]

proRRt37_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,36]


proRRt38_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,37]



proRRt39_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,38]

proRRt40_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,39]




proRRt41_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,40]

proRRt42_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,41]

proRRt43_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,42]


proRRt44_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,43]



proRRt45_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,44]


proRRt46_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,45]



proRRt47_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,46]


proRRt48_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,47]

proRRt49_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,48]


proRRt50_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,49]


proRRt51_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,50]


proRRt52_lag_ens_mean_inst_44[1]=proRR_lag1_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[2]=proRR_lag2_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[3]=proRR_lag3_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[4]=proRR_lag4_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[5]=proRR_lag5_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[6]=proRR_lag6_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[7]=proRR_lag7_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[8]=proRR_lag8_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[9]=proRR_lag9_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[10]=proRR_lag10_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[11]=proRR_lag11_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[12]=proRR_lag12_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[13]=proRR_lag13_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[14]=proRR_lag14_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[15]=proRR_lag15_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[16]=proRR_lag16_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[17]=proRR_lag17_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[18]=proRR_lag18_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[19]=proRR_lag19_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[20]=proRR_lag20_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[21]=proRR_lag21_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[22]=proRR_lag22_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[23]=proRR_lag23_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[24]=proRR_lag24_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[25]=proRR_lag25_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[26]=proRR_lag26_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[27]=proRR_lag27_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[28]=proRR_lag28_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[29]=proRR_lag29_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[30]=proRR_lag30_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[31]=proRR_lag31_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_lag_ens_mean_inst_44[32]=proRR_lag32_ens_mean_inst_44_vs_AveRR_thr[0,51]


proRRt5_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt6_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt7_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt8_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt9_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt10_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt11_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt12_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt13_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt14_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt15_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt16_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt17_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt18_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt19_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt20_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt21_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt22_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt23_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt24_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt25_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt26_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt27_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt28_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt29_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt30_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt31_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt32_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt33_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt34_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt35_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt36_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt37_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt38_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt39_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt40_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt41_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt42_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt43_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt44_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt45_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt46_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt47_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt48_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt49_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt50_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt51_proR2_lag_ens_mean_inst_44=N.zeros((nlag))
proRRt52_proR2_lag_ens_mean_inst_44=N.zeros((nlag))


proRRt5_proR2_lag_ens_mean_inst_44[0]=proRRt5_lag_ens_mean_inst_44[0]-proR2t5_lag_ens_mean_inst_44[0]
proRRt6_proR2_lag_ens_mean_inst_44[0]=proRRt6_lag_ens_mean_inst_44[0]-proR2t6_lag_ens_mean_inst_44[0]
proRRt7_proR2_lag_ens_mean_inst_44[0]=proRRt7_lag_ens_mean_inst_44[0]-proR2t7_lag_ens_mean_inst_44[0]
proRRt8_proR2_lag_ens_mean_inst_44[0]=proRRt8_lag_ens_mean_inst_44[0]-proR2t8_lag_ens_mean_inst_44[0]
proRRt9_proR2_lag_ens_mean_inst_44[0]=proRRt9_lag_ens_mean_inst_44[0]-proR2t9_lag_ens_mean_inst_44[0]
proRRt10_proR2_lag_ens_mean_inst_44[0]=proRRt10_lag_ens_mean_inst_44[0]-proR2t10_lag_ens_mean_inst_44[0]
proRRt11_proR2_lag_ens_mean_inst_44[0]=proRRt11_lag_ens_mean_inst_44[0]-proR2t11_lag_ens_mean_inst_44[0]
proRRt12_proR2_lag_ens_mean_inst_44[0]=proRRt12_lag_ens_mean_inst_44[0]-proR2t12_lag_ens_mean_inst_44[0]
proRRt13_proR2_lag_ens_mean_inst_44[0]=proRRt13_lag_ens_mean_inst_44[0]-proR2t13_lag_ens_mean_inst_44[0]
proRRt14_proR2_lag_ens_mean_inst_44[0]=proRRt14_lag_ens_mean_inst_44[0]-proR2t14_lag_ens_mean_inst_44[0]
proRRt15_proR2_lag_ens_mean_inst_44[0]=proRRt15_lag_ens_mean_inst_44[0]-proR2t15_lag_ens_mean_inst_44[0]
proRRt16_proR2_lag_ens_mean_inst_44[0]=proRRt16_lag_ens_mean_inst_44[0]-proR2t16_lag_ens_mean_inst_44[0]
proRRt17_proR2_lag_ens_mean_inst_44[0]=proRRt17_lag_ens_mean_inst_44[0]-proR2t17_lag_ens_mean_inst_44[0]
proRRt18_proR2_lag_ens_mean_inst_44[0]=proRRt18_lag_ens_mean_inst_44[0]-proR2t18_lag_ens_mean_inst_44[0]
proRRt19_proR2_lag_ens_mean_inst_44[0]=proRRt19_lag_ens_mean_inst_44[0]-proR2t19_lag_ens_mean_inst_44[0]
proRRt20_proR2_lag_ens_mean_inst_44[0]=proRRt20_lag_ens_mean_inst_44[0]-proR2t20_lag_ens_mean_inst_44[0]
proRRt21_proR2_lag_ens_mean_inst_44[0]=proRRt21_lag_ens_mean_inst_44[0]-proR2t21_lag_ens_mean_inst_44[0]
proRRt22_proR2_lag_ens_mean_inst_44[0]=proRRt22_lag_ens_mean_inst_44[0]-proR2t22_lag_ens_mean_inst_44[0]
proRRt23_proR2_lag_ens_mean_inst_44[0]=proRRt23_lag_ens_mean_inst_44[0]-proR2t23_lag_ens_mean_inst_44[0]
proRRt24_proR2_lag_ens_mean_inst_44[0]=proRRt24_lag_ens_mean_inst_44[0]-proR2t24_lag_ens_mean_inst_44[0]
proRRt25_proR2_lag_ens_mean_inst_44[0]=proRRt25_lag_ens_mean_inst_44[0]-proR2t25_lag_ens_mean_inst_44[0]
proRRt26_proR2_lag_ens_mean_inst_44[0]=proRRt26_lag_ens_mean_inst_44[0]-proR2t26_lag_ens_mean_inst_44[0]
proRRt27_proR2_lag_ens_mean_inst_44[0]=proRRt27_lag_ens_mean_inst_44[0]-proR2t27_lag_ens_mean_inst_44[0]
proRRt28_proR2_lag_ens_mean_inst_44[0]=proRRt28_lag_ens_mean_inst_44[0]-proR2t28_lag_ens_mean_inst_44[0]
proRRt29_proR2_lag_ens_mean_inst_44[0]=proRRt29_lag_ens_mean_inst_44[0]-proR2t29_lag_ens_mean_inst_44[0]
proRRt30_proR2_lag_ens_mean_inst_44[0]=proRRt30_lag_ens_mean_inst_44[0]-proR2t30_lag_ens_mean_inst_44[0]
proRRt31_proR2_lag_ens_mean_inst_44[0]=proRRt31_lag_ens_mean_inst_44[0]-proR2t31_lag_ens_mean_inst_44[0]
proRRt32_proR2_lag_ens_mean_inst_44[0]=proRRt32_lag_ens_mean_inst_44[0]-proR2t32_lag_ens_mean_inst_44[0]
proRRt33_proR2_lag_ens_mean_inst_44[0]=proRRt33_lag_ens_mean_inst_44[0]-proR2t33_lag_ens_mean_inst_44[0]
proRRt34_proR2_lag_ens_mean_inst_44[0]=proRRt34_lag_ens_mean_inst_44[0]-proR2t34_lag_ens_mean_inst_44[0]
proRRt35_proR2_lag_ens_mean_inst_44[0]=proRRt35_lag_ens_mean_inst_44[0]-proR2t35_lag_ens_mean_inst_44[0]
proRRt36_proR2_lag_ens_mean_inst_44[0]=proRRt36_lag_ens_mean_inst_44[0]-proR2t36_lag_ens_mean_inst_44[0]
proRRt37_proR2_lag_ens_mean_inst_44[0]=proRRt37_lag_ens_mean_inst_44[0]-proR2t37_lag_ens_mean_inst_44[0]
proRRt38_proR2_lag_ens_mean_inst_44[0]=proRRt38_lag_ens_mean_inst_44[0]-proR2t38_lag_ens_mean_inst_44[0]
proRRt39_proR2_lag_ens_mean_inst_44[0]=proRRt39_lag_ens_mean_inst_44[0]-proR2t39_lag_ens_mean_inst_44[0]
proRRt40_proR2_lag_ens_mean_inst_44[0]=proRRt40_lag_ens_mean_inst_44[0]-proR2t40_lag_ens_mean_inst_44[0]
proRRt41_proR2_lag_ens_mean_inst_44[0]=proRRt41_lag_ens_mean_inst_44[0]-proR2t41_lag_ens_mean_inst_44[0]
proRRt42_proR2_lag_ens_mean_inst_44[0]=proRRt42_lag_ens_mean_inst_44[0]-proR2t42_lag_ens_mean_inst_44[0]
proRRt43_proR2_lag_ens_mean_inst_44[0]=proRRt43_lag_ens_mean_inst_44[0]-proR2t43_lag_ens_mean_inst_44[0]
proRRt44_proR2_lag_ens_mean_inst_44[0]=proRRt44_lag_ens_mean_inst_44[0]-proR2t44_lag_ens_mean_inst_44[0]
proRRt45_proR2_lag_ens_mean_inst_44[0]=proRRt45_lag_ens_mean_inst_44[0]-proR2t45_lag_ens_mean_inst_44[0]
proRRt46_proR2_lag_ens_mean_inst_44[0]=proRRt46_lag_ens_mean_inst_44[0]-proR2t46_lag_ens_mean_inst_44[0]
proRRt47_proR2_lag_ens_mean_inst_44[0]=proRRt47_lag_ens_mean_inst_44[0]-proR2t47_lag_ens_mean_inst_44[0]
proRRt48_proR2_lag_ens_mean_inst_44[0]=proRRt48_lag_ens_mean_inst_44[0]-proR2t48_lag_ens_mean_inst_44[0]
proRRt49_proR2_lag_ens_mean_inst_44[0]=proRRt49_lag_ens_mean_inst_44[0]-proR2t49_lag_ens_mean_inst_44[0]
proRRt50_proR2_lag_ens_mean_inst_44[0]=proRRt50_lag_ens_mean_inst_44[0]-proR2t50_lag_ens_mean_inst_44[0]
proRRt51_proR2_lag_ens_mean_inst_44[0]=proRRt51_lag_ens_mean_inst_44[0]-proR2t51_lag_ens_mean_inst_44[0]
proRRt52_proR2_lag_ens_mean_inst_44[0]=proRRt52_lag_ens_mean_inst_44[0]-proR2t52_lag_ens_mean_inst_44[0]

proRRt5_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,4]
proRRt5_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,4]

proRRt6_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,5]

proRRt6_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,5]
proRRt6_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,5]

proRRt7_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,6]
proRRt7_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,6]

proRRt8_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,7]
proRRt8_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,7]

proRRt9_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,8]
proRRt9_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,8]





proRRt10_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,9]

proRRt10_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,9]
proRRt10_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,9]


proRRt11_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,10]
proRRt11_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,10]


proRRt12_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,11]
proRRt12_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,11]


proRRt13_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,12]

proRRt14_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,13]

proRRt15_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,14]

proRRt16_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,15]

proRRt17_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,16]

proRRt18_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,17]

proRRt19_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,18]

proRRt20_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,19]

proRRt21_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,20]

proRRt22_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,21]

proRRt23_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,22]

proRRt24_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,23]

proRRt25_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,24]

proRRt26_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,25]

proRRt27_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,26]

proRRt28_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,27]

proRRt13_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,12]
proRRt13_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,12]


proRRt14_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,13]
proRRt14_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,13]

proRRt15_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,14]
proRRt15_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,14]


proRRt16_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,15]
proRRt16_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,15]

proRRt17_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,16]
proRRt17_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,16]

proRRt18_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,17]
proRRt18_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,17]


proRRt19_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,18]
proRRt19_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,18]

proRRt20_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,19]
proRRt20_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,19]


proRRt21_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,20]
proRRt21_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,20]

proRRt22_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,21]
proRRt22_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,21]

proRRt23_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,22]
proRRt23_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,22]

proRRt24_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,23]
proRRt24_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,23]

proRRt25_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,24]
proRRt25_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,24]


proRRt26_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,25]
proRRt26_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,25]

proRRt27_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,26]
proRRt27_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,26]


proRRt28_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,27]
proRRt28_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,27]


proRRt29_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,28]
proRRt29_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,28]

proRRt30_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,29]
proRRt30_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,29]

proRRt31_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,30]
proRRt31_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,30]






proRRt32_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,31]
proRRt32_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,31]


proRRt33_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,32]
proRRt33_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,32]







proRRt34_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,33]
proRRt34_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,33]







proRRt35_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,34]
proRRt35_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,34]



proRRt36_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,35]
proRRt36_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,35]

proRRt37_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,36]
proRRt37_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,36]


proRRt38_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,37]
proRRt38_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,37]



proRRt39_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,38]
proRRt39_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,38]

proRRt40_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,39]
proRRt40_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,39]




proRRt41_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,40]
proRRt41_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,40]

proRRt42_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,41]
proRRt42_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,41]

proRRt43_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,42]
proRRt43_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,42]


proRRt44_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,43]
proRRt44_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,43]



proRRt45_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,44]
proRRt45_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,44]


proRRt46_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,45]
proRRt46_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,45]



proRRt47_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,46]
proRRt47_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,46]


proRRt48_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,47]
proRRt48_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,47]

proRRt49_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,48]
proRRt49_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,48]


proRRt50_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,49]
proRRt50_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,49]


proRRt51_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,50]
proRRt51_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,50]


proRRt52_proR2_lag_ens_mean_inst_44[1]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[2]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[3]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[4]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[5]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[6]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[7]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[8]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[9]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[10]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[11]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[12]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[13]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[14]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[15]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[16]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[17]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[18]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[19]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[20]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[21]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[22]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[23]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[24]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[25]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[26]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[27]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[28]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[29]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[30]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[31]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,51]
proRRt52_proR2_lag_ens_mean_inst_44[32]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,51]


proRR_ProR2t_lag_dt=N.zeros((nlag, n_30min))
proRR_ProR2t_lag_dt[0,:]=float('nan')
for j in N.arange(48):
     proRR_ProR2t_lag_dt[1,j]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[2,j]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[3,j]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[4,j]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[5,j]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[6,j]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[7,j]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[8,j]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[9,j]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[10,j]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[11,j]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[12,j]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[13,j]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[14,j]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[15,j]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[16,j]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[17,j]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[18,j]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[19,j]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[20,j]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[21,j]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[22,j]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[23,j]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[24,j]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[25,j]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[26,j]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[27,j]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[28,j]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[29,j]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[30,j]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[31,j]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[32,j]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,j]
     if ( j <=11):
          proRR_ProR2t_lag_dt[:,j]=float('nan')
     if ( j == 12):
          proRR_ProR2t_lag_dt[5:nlag,j]=float('nan')
     if ( j == 13):
          proRR_ProR2t_lag_dt[6:nlag,j]=float('nan')
     if ( j == 14):
          proRR_ProR2t_lag_dt[7:nlag,j]=float('nan')
     if ( j == 15):
          proRR_ProR2t_lag_dt[8:nlag,j]=float('nan')
     if ( j == 16):
          proRR_ProR2t_lag_dt[9:nlag,j]=float('nan')
     if ( j == 17):
          proRR_ProR2t_lag_dt[10:nlag,j]=float('nan')
     if ( j == 18):
          proRR_ProR2t_lag_dt[11:nlag,j]=float('nan')
     if ( j == 19):
          proRR_ProR2t_lag_dt[12:nlag,j]=float('nan')
     if ( j == 20):
          proRR_ProR2t_lag_dt[13:nlag,j]=float('nan')
     if ( j == 21):
          proRR_ProR2t_lag_dt[14:nlag,j]=float('nan')
     if ( j == 22):
          proRR_ProR2t_lag_dt[15:nlag,j]=float('nan')
     if ( j == 23):
          proRR_ProR2t_lag_dt[16:nlag,j]=float('nan')
     if ( j == 24):
          proRR_ProR2t_lag_dt[17:nlag,j]=float('nan')
     if ( j == 25):
          proRR_ProR2t_lag_dt[18:nlag,j]=float('nan')
     if ( j == 26):
          proRR_ProR2t_lag_dt[19:nlag,j]=float('nan')
     if ( j == 27):
          proRR_ProR2t_lag_dt[20:nlag,j]=float('nan')
     if ( j == 28):
          proRR_ProR2t_lag_dt[21:nlag,j]=float('nan')
     if ( j == 29):
          proRR_ProR2t_lag_dt[22:nlag,j]=float('nan')
     if ( j == 30):
          proRR_ProR2t_lag_dt[23:nlag,j]=float('nan')
     if ( j == 31):
          proRR_ProR2t_lag_dt[24:nlag,j]=float('nan')
     if ( j == 32):
          proRR_ProR2t_lag_dt[25:nlag,j]=float('nan')
     if ( j == 33):
          proRR_ProR2t_lag_dt[26:nlag,j]=float('nan')
     if ( j == 34):
          proRR_ProR2t_lag_dt[27:nlag,j]=float('nan')
     if ( j == 35):
          proRR_ProR2t_lag_dt[28:nlag,j]=float('nan')
     if ( j == 36):
          proRR_ProR2t_lag_dt[29:nlag,j]=float('nan')
     if ( j == 37):
          proRR_ProR2t_lag_dt[30:nlag,j]=float('nan')
     if ( j == 38):
          proRR_ProR2t_lag_dt[31:nlag,j]=float('nan')
     if ( j == 39):
          proRR_ProR2t_lag_dt[32:nlag,j]=float('nan')
     if ( j == 40):
          proRR_ProR2t_lag_dt[33:nlag,j]=float('nan')
     if ( j> 40):
          proRR_ProR2t_lag_dt[34:nlag,j]=float('nan')


np.save('proRR_ProR2t_lag_dt_CTRL',proRR_ProR2t_lag_dt)
