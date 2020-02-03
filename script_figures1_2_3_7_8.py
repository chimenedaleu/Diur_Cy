# Import libraries
from mpl_toolkits.basemap import Basemap, cm
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
L_sub=2.85e6 #2.834e6

nz=99
ny=512
nx=512
dx=200
dy=200
length_x=N.zeros(nx)
length_y=N.zeros(ny)
length_xB=N.zeros(nx)
length_yB=N.zeros(ny)

length_x[0]=00.
length_y[0]=00.
for j in N.arange(ny-1):
     length_x[j+1]=length_x[j]+dx
     length_y[j+1]=length_y[j]+dy

#length_x=length_x
#length_y=length_y
length_xB=length_x*2.5
length_yB=length_y*2.5


hein=N.zeros((nz))
hein = np.load('hein.npy')
time_SF = np.load('time_SF.npy')

#
#control
#
#Timeseries 15 min instanteneous
time_10min_CTRL_f = np.load('modets_time_10min_CTRL.npy')
surf_flx_CTRL_f = np.load('surf_flx_CTRL.npy')
precip_CTRL_f = np.load('surf_precip_time_CTRL.npy')
cld_MF_cb_CTRL_f = np.load('cloud_MF_time_cb5_CTRL.npy')
cld_frac_cb_CTRL_f = np.load('cloud_frac_time_cb5_CTRL.npy')
cldBCu_MF_cb_CTRL_f = np.load('BCu_MF_time_cb5_CTRL.npy')
cldBCu_frac_cb_CTRL_f = np.load('BCu_frac_time_cb5_CTRL.npy')
precip_CTRL_f=precip_CTRL_f*3600.
ACuzt_frac_time_CTRL_f = np.load('ACuzt_frac_time_CTRL.npy')
BCuzt_frac_time_CTRL_f = np.load('BCuzt_frac_time_CTRL.npy')
ACuzt_MF_time_CTRL_f = np.load('ACuzt_MF_time_CTRL.npy')
BCuzt_MF_time_CTRL_f = np.load('BCuzt_MF_time_CTRL.npy')

th_mean_CTRL_f = np.load('th_mean_CTRL.npy')
time_30min_CTRL_f = np.load('modets_time_30min_CTRL.npy')
thzy24h_snap_CTRL_f = np.load('thzy24h_snap.npy')
qvzy24h_snap_CTRL_f = np.load('qvzy24h_snap.npy')

thxy24h_1km_snap_CTRL_f =np.load('thxy24h_1km_CTRL_snap.npy')
thxy24h_3km_snap_CTRL_f =np.load('thxy24h_3km_CTRL_snap.npy')
qvxy24h_1km_snap_CTRL_f =np.load('qvxy24h_1km_CTRL_snap.npy')
qvxy24h_3km_snap_CTRL_f =np.load('qvxy24h_3km_CTRL_snap.npy')

#
#Strong
#
time_10min_PHALF_f = np.load('modets_time_10min_S195L600.npy')
surf_flx_PHALF_f = np.load('surf_flx_S195L600.npy')
precip_PHALF_f = np.load('surf_precip_time_S195L600.npy')
cld_MF_cb_PHALF_f = np.load('cloud_MF_time_cb5_S195L600.npy')
precip_PHALF_f=precip_PHALF_f*3600.


#
#Weak
#
time_10min_MHALF_f = np.load('modets_time_10min_S65L200.npy')
surf_flx_MHALF_f = np.load('surf_flx_S65L200.npy')
precip_MHALF_f = np.load('surf_precip_time_S65L200.npy')
cld_MF_cb_MHALF_f = np.load('cloud_MF_time_cb5_S65L200.npy')
precip_MHALF_f=precip_MHALF_f*3600.

#homogenization perturbation applied between 15-24h of the 1st diurnal cycle
#simulation on the second cycle should be compared to the control simulation on the second diurnal cycle
#
time_10min_RTHQVD2 = np.load('modets_time_10min_RTHQVD2.npy')
precip_RTHQVD2 = np.load('surf_precip_time_RTHQVD2.npy')
cld_MF_cb_RTHQVD2 = np.load('cloud_MF_time_cb5_RTHQVD2.npy')
precip_RTHQVD2=precip_RTHQVD2*3600.

#
#rthqvd3: 
#
time_10min_RTHQVD3 = np.load('modets_time_10min_RTHQVD3.npy')
precip_RTHQVD3 = np.load('surf_precip_time_RTHQVD3.npy')
cld_MF_cb_RTHQVD3 = np.load('cloud_MF_time_cb5_RTHQVD3.npy')
precip_RTHQVD3=precip_RTHQVD3*3600.

#
#rthqvd4
#
time_10min_RTHQVD4 = np.load('modets_time_10min_RTHQVD4.npy')
precip_RTHQVD4 = np.load('surf_precip_time_RTHQVD4.npy')
cld_MF_cb_RTHQVD4 = np.load('cloud_MF_time_cb5_RTHQVD4.npy')
precip_RTHQVD4=precip_RTHQVD4*3600.

#
#rthqvd5
#
time_10min_RTHQVD5 = np.load('modets_time_10min_RTHQVD5.npy')
precip_RTHQVD5 = np.load('surf_precip_time_RTHQVD5.npy')
cld_MF_cb_RTHQVD5 = np.load('cloud_MF_time_cb5_RTHQVD5.npy')
precip_RTHQVD5=precip_RTHQVD5*3600.

#
#rthqvd6
#
time_10min_RTHQVD6 = np.load('modets_time_10min_RTHQVD6.npy')
precip_RTHQVD6 = np.load('surf_precip_time_RTHQVD6.npy')
cld_MF_cb_RTHQVD6 = np.load('cloud_MF_time_cb5_RTHQVD6.npy')
precip_RTHQVD6=precip_RTHQVD6*3600.

#
#rthqvd7
#
time_10min_RTHQVD7 = np.load('modets_time_10min_RTHQVD7.npy')
precip_RTHQVD7 = np.load('surf_precip_time_RTHQVD7.npy')
cld_MF_cb_RTHQVD7 = np.load('cloud_MF_time_cb5_RTHQVD7.npy')
precip_RTHQVD7=precip_RTHQVD7*3600.

#
#rthqvd8
#
time_10min_RTHQVD8 = np.load('modets_time_10min_RTHQVD8.npy')
precip_RTHQVD8 = np.load('surf_precip_time_RTHQVD8.npy')
cld_MF_cb_RTHQVD8 = np.load('cloud_MF_time_cb5_RTHQVD8.npy')
precip_RTHQVD8=precip_RTHQVD8*3600.


#
#rthqvd9
#
time_10min_RTHQVD9 = np.load('modets_time_10min_RTHQVD9.npy')
precip_RTHQVD9 = np.load('surf_precip_time_RTHQVD9.npy')
cld_MF_cb_RTHQVD9 = np.load('cloud_MF_time_cb5_RTHQVD9.npy')
precip_RTHQVD9=precip_RTHQVD9*3600.

#
#rthqvd10
#
time_10min_RTHQVD10 = np.load('modets_time_10min_RTHQVD10.npy')
precip_RTHQVD10 = np.load('surf_precip_time_RTHQVD10.npy')
cld_MF_cb_RTHQVD10 = np.load('cloud_MF_time_cb5_RTHQVD10.npy')
precip_RTHQVD10=precip_RTHQVD10*3600.


n_f_10min=93

time_10min_CTRL_f_d1=N.zeros((n_f_10min))
time_10min_CTRL_f_d2=N.zeros((n_f_10min))
time_10min_CTRL_f_d3=N.zeros((n_f_10min))
time_10min_CTRL_f_d4=N.zeros((n_f_10min))
time_10min_CTRL_f_d5=N.zeros((n_f_10min))
time_10min_CTRL_f_d6=N.zeros((n_f_10min))
time_10min_CTRL_f_d7=N.zeros((n_f_10min))
time_10min_CTRL_f_d8=N.zeros((n_f_10min))
time_10min_CTRL_f_d9=N.zeros((n_f_10min))
time_10min_CTRL_f_d10=N.zeros((n_f_10min))

precip_CTRL_f_d1=N.zeros((n_f_10min))
precip_CTRL_f_d2=N.zeros((n_f_10min))
precip_CTRL_f_d3=N.zeros((n_f_10min))
precip_CTRL_f_d4=N.zeros((n_f_10min))
precip_CTRL_f_d5=N.zeros((n_f_10min))
precip_CTRL_f_d6=N.zeros((n_f_10min))
precip_CTRL_f_d7=N.zeros((n_f_10min))
precip_CTRL_f_d8=N.zeros((n_f_10min))
precip_CTRL_f_d9=N.zeros((n_f_10min))
precip_CTRL_f_d10=N.zeros((n_f_10min))

cld_MF_cb_CTRL_f_d1=N.zeros((n_f_10min))
cld_MF_cb_CTRL_f_d2=N.zeros((n_f_10min))
cld_MF_cb_CTRL_f_d3=N.zeros((n_f_10min))
cld_MF_cb_CTRL_f_d4=N.zeros((n_f_10min))
cld_MF_cb_CTRL_f_d5=N.zeros((n_f_10min))
cld_MF_cb_CTRL_f_d6=N.zeros((n_f_10min))
cld_MF_cb_CTRL_f_d7=N.zeros((n_f_10min))
cld_MF_cb_CTRL_f_d8=N.zeros((n_f_10min))
cld_MF_cb_CTRL_f_d9=N.zeros((n_f_10min))
cld_MF_cb_CTRL_f_d10=N.zeros((n_f_10min))

cld_frac_cb_CTRL_f_d1=N.zeros((n_f_10min))
cld_frac_cb_CTRL_f_d2=N.zeros((n_f_10min))
cld_frac_cb_CTRL_f_d3=N.zeros((n_f_10min))
cld_frac_cb_CTRL_f_d4=N.zeros((n_f_10min))
cld_frac_cb_CTRL_f_d5=N.zeros((n_f_10min))
cld_frac_cb_CTRL_f_d6=N.zeros((n_f_10min))
cld_frac_cb_CTRL_f_d7=N.zeros((n_f_10min))
cld_frac_cb_CTRL_f_d8=N.zeros((n_f_10min))
cld_frac_cb_CTRL_f_d9=N.zeros((n_f_10min))
cld_frac_cb_CTRL_f_d10=N.zeros((n_f_10min))


cldBCu_MF_cb_CTRL_f_d1=N.zeros((n_f_10min))
cldBCu_MF_cb_CTRL_f_d2=N.zeros((n_f_10min))
cldBCu_MF_cb_CTRL_f_d3=N.zeros((n_f_10min))
cldBCu_MF_cb_CTRL_f_d4=N.zeros((n_f_10min))
cldBCu_MF_cb_CTRL_f_d5=N.zeros((n_f_10min))
cldBCu_MF_cb_CTRL_f_d6=N.zeros((n_f_10min))
cldBCu_MF_cb_CTRL_f_d7=N.zeros((n_f_10min))
cldBCu_MF_cb_CTRL_f_d8=N.zeros((n_f_10min))
cldBCu_MF_cb_CTRL_f_d9=N.zeros((n_f_10min))
cldBCu_MF_cb_CTRL_f_d10=N.zeros((n_f_10min))

cldBCu_frac_cb_CTRL_f_d1=N.zeros((n_f_10min))
cldBCu_frac_cb_CTRL_f_d2=N.zeros((n_f_10min))
cldBCu_frac_cb_CTRL_f_d3=N.zeros((n_f_10min))
cldBCu_frac_cb_CTRL_f_d4=N.zeros((n_f_10min))
cldBCu_frac_cb_CTRL_f_d5=N.zeros((n_f_10min))
cldBCu_frac_cb_CTRL_f_d6=N.zeros((n_f_10min))
cldBCu_frac_cb_CTRL_f_d7=N.zeros((n_f_10min))
cldBCu_frac_cb_CTRL_f_d8=N.zeros((n_f_10min))
cldBCu_frac_cb_CTRL_f_d9=N.zeros((n_f_10min))
cldBCu_frac_cb_CTRL_f_d10=N.zeros((n_f_10min))

time_10min_CTRL_f_d1[0:92]=time_10min_CTRL_f[0:92]
time_10min_CTRL_f_d2[0:92]=time_10min_CTRL_f[90:182]-24.0
time_10min_CTRL_f_d3[0:92]=time_10min_CTRL_f[182:274]-48.0
time_10min_CTRL_f_d4[0:92]=time_10min_CTRL_f[273:365]-72.0
time_10min_CTRL_f_d5[0:92]=time_10min_CTRL_f[363:455]-96.0
time_10min_CTRL_f_d6[0:92]=time_10min_CTRL_f[455:547]-120.0
time_10min_CTRL_f_d7[0:92]=time_10min_CTRL_f[547:639]-144.0
time_10min_CTRL_f_d8[0:92]=time_10min_CTRL_f[637:729]-168.0
time_10min_CTRL_f_d9[0:92]=time_10min_CTRL_f[728:820]-192.0
time_10min_CTRL_f_d10[0:92]=time_10min_CTRL_f[819:911]-216.0


precip_CTRL_f_d1[0:92]=precip_CTRL_f[0:92]
precip_CTRL_f_d2[0:92]=precip_CTRL_f[90:182]
precip_CTRL_f_d3[0:92]=precip_CTRL_f[182:274]
precip_CTRL_f_d4[0:92]=precip_CTRL_f[273:365]
precip_CTRL_f_d5[0:92]=precip_CTRL_f[363:455]
precip_CTRL_f_d6[0:92]=precip_CTRL_f[455:547]
precip_CTRL_f_d7[0:92]=precip_CTRL_f[547:639]
precip_CTRL_f_d8[0:92]=precip_CTRL_f[637:729]
precip_CTRL_f_d9[0:92]=precip_CTRL_f[728:820]
precip_CTRL_f_d10[0:92]=precip_CTRL_f[819:911]

print 'mrr ctrl', (precip_CTRL_f_d1[0:79].mean()+  precip_CTRL_f_d2[0:79].mean()+ precip_CTRL_f_d3[0:79].mean()+ precip_CTRL_f_d4[0:79].mean()+ precip_CTRL_f_d5[0:79].mean()+ precip_CTRL_f_d6[0:79].mean()+ precip_CTRL_f_d7[0:79].mean()+ precip_CTRL_f_d8[0:79].mean()+ precip_CTRL_f_d9[0:79].mean()+ precip_CTRL_f_d10[0:79].mean())/10.


cld_MF_cb_CTRL_f_d1[0:92]=cld_MF_cb_CTRL_f[0:92]
cld_MF_cb_CTRL_f_d2[0:92]=cld_MF_cb_CTRL_f[90:182]
cld_MF_cb_CTRL_f_d3[0:92]=cld_MF_cb_CTRL_f[182:274]
cld_MF_cb_CTRL_f_d4[0:92]=cld_MF_cb_CTRL_f[273:365]
cld_MF_cb_CTRL_f_d5[0:92]=cld_MF_cb_CTRL_f[363:455]
cld_MF_cb_CTRL_f_d6[0:92]=cld_MF_cb_CTRL_f[455:547]
cld_MF_cb_CTRL_f_d7[0:92]=cld_MF_cb_CTRL_f[547:639]
cld_MF_cb_CTRL_f_d8[0:92]=cld_MF_cb_CTRL_f[637:729]
cld_MF_cb_CTRL_f_d9[0:92]=cld_MF_cb_CTRL_f[728:820]
cld_MF_cb_CTRL_f_d10[0:92]=cld_MF_cb_CTRL_f[819:911]

cld_frac_cb_CTRL_f_d1[0:92]=cld_frac_cb_CTRL_f[0:92]
cld_frac_cb_CTRL_f_d2[0:92]=cld_frac_cb_CTRL_f[90:182]
cld_frac_cb_CTRL_f_d3[0:92]=cld_frac_cb_CTRL_f[182:274]
cld_frac_cb_CTRL_f_d4[0:92]=cld_frac_cb_CTRL_f[273:365]
cld_frac_cb_CTRL_f_d5[0:92]=cld_frac_cb_CTRL_f[363:455]
cld_frac_cb_CTRL_f_d6[0:92]=cld_frac_cb_CTRL_f[455:547]
cld_frac_cb_CTRL_f_d7[0:92]=cld_frac_cb_CTRL_f[547:639]
cld_frac_cb_CTRL_f_d8[0:92]=cld_frac_cb_CTRL_f[637:729]
cld_frac_cb_CTRL_f_d9[0:92]=cld_frac_cb_CTRL_f[728:820]
cld_frac_cb_CTRL_f_d10[0:92]=cld_frac_cb_CTRL_f[819:911]

cldBCu_MF_cb_CTRL_f_d1[0:92]=cldBCu_MF_cb_CTRL_f[0:92]
cldBCu_MF_cb_CTRL_f_d2[0:92]=cldBCu_MF_cb_CTRL_f[90:182]
cldBCu_MF_cb_CTRL_f_d3[0:92]=cldBCu_MF_cb_CTRL_f[182:274]
cldBCu_MF_cb_CTRL_f_d4[0:92]=cldBCu_MF_cb_CTRL_f[273:365]
cldBCu_MF_cb_CTRL_f_d5[0:92]=cldBCu_MF_cb_CTRL_f[363:455]
cldBCu_MF_cb_CTRL_f_d6[0:92]=cldBCu_MF_cb_CTRL_f[455:547]
cldBCu_MF_cb_CTRL_f_d7[0:92]=cldBCu_MF_cb_CTRL_f[547:639]
cldBCu_MF_cb_CTRL_f_d8[0:92]=cldBCu_MF_cb_CTRL_f[637:729]
cldBCu_MF_cb_CTRL_f_d9[0:92]=cldBCu_MF_cb_CTRL_f[728:820]
cldBCu_MF_cb_CTRL_f_d10[0:92]=cldBCu_MF_cb_CTRL_f[819:911]

cldBCu_frac_cb_CTRL_f_d1[0:92]=cldBCu_frac_cb_CTRL_f[0:92]
cldBCu_frac_cb_CTRL_f_d2[0:92]=cldBCu_frac_cb_CTRL_f[90:182]
cldBCu_frac_cb_CTRL_f_d3[0:92]=cldBCu_frac_cb_CTRL_f[182:274]
cldBCu_frac_cb_CTRL_f_d4[0:92]=cldBCu_frac_cb_CTRL_f[273:365]
cldBCu_frac_cb_CTRL_f_d5[0:92]=cldBCu_frac_cb_CTRL_f[363:455]
cldBCu_frac_cb_CTRL_f_d6[0:92]=cldBCu_frac_cb_CTRL_f[455:547]
cldBCu_frac_cb_CTRL_f_d7[0:92]=cldBCu_frac_cb_CTRL_f[547:639]
cldBCu_frac_cb_CTRL_f_d8[0:92]=cldBCu_frac_cb_CTRL_f[637:729]
cldBCu_frac_cb_CTRL_f_d9[0:92]=cldBCu_frac_cb_CTRL_f[728:820]
cldBCu_frac_cb_CTRL_f_d10[0:92]=cldBCu_frac_cb_CTRL_f[819:911]

nr=36+n_f_10min
time_10min_CTRL_f_d1_R=N.zeros((nr))
time_10min_CTRL_f_d2_R=N.zeros((nr))
time_10min_CTRL_f_d3_R=N.zeros((nr))
time_10min_CTRL_f_d4_R=N.zeros((nr))
time_10min_CTRL_f_d5_R=N.zeros((nr))
time_10min_CTRL_f_d6_R=N.zeros((nr))
time_10min_CTRL_f_d7_R=N.zeros((nr))
time_10min_CTRL_f_d8_R=N.zeros((nr))
time_10min_CTRL_f_d9_R=N.zeros((nr))
time_10min_CTRL_f_d10_R=N.zeros((nr))

time_10min_CTRL_f_d1_R[0]=-9
for j in N.arange(36):
     time_10min_CTRL_f_d1_R[j+1]=time_10min_CTRL_f_d1_R[j]+0.25

time_10min_CTRL_f_d2_R[0:37]=time_10min_CTRL_f_d1_R[0:37]
time_10min_CTRL_f_d3_R[0:37]=time_10min_CTRL_f_d1_R[0:37]
time_10min_CTRL_f_d4_R[0:37]=time_10min_CTRL_f_d1_R[0:37]
time_10min_CTRL_f_d5_R[0:37]=time_10min_CTRL_f_d1_R[0:37]
time_10min_CTRL_f_d6_R[0:37]=time_10min_CTRL_f_d1_R[0:37]
time_10min_CTRL_f_d7_R[0:37]=time_10min_CTRL_f_d1_R[0:37]
time_10min_CTRL_f_d8_R[0:37]=time_10min_CTRL_f_d1_R[0:37]
time_10min_CTRL_f_d9_R[0:37]=time_10min_CTRL_f_d1_R[0:37]
time_10min_CTRL_f_d10_R[0:37]=time_10min_CTRL_f_d1_R[0:37]

time_10min_CTRL_f_d1_R[37:nr]=time_10min_CTRL_f_d1[0:92]
time_10min_CTRL_f_d2_R[37:nr]=time_10min_CTRL_f_d2[0:92]
time_10min_CTRL_f_d3_R[37:nr]=time_10min_CTRL_f_d3[0:92]
time_10min_CTRL_f_d4_R[37:nr]=time_10min_CTRL_f_d4[0:92]
time_10min_CTRL_f_d5_R[37:nr]=time_10min_CTRL_f_d5[0:92]
time_10min_CTRL_f_d6_R[37:nr]=time_10min_CTRL_f_d6[0:92]
time_10min_CTRL_f_d7_R[37:nr]=time_10min_CTRL_f_d7[0:92]
time_10min_CTRL_f_d8_R[37:nr]=time_10min_CTRL_f_d8[0:92]
time_10min_CTRL_f_d9_R[37:nr]=time_10min_CTRL_f_d9[0:92]
time_10min_CTRL_f_d10_R[37:nr]=time_10min_CTRL_f_d10[0:92]

precip_CTRL_f_d1_R=N.zeros((nr))
precip_CTRL_f_d2_R=N.zeros((nr))
precip_CTRL_f_d3_R=N.zeros((nr))
precip_CTRL_f_d4_R=N.zeros((nr))
precip_CTRL_f_d5_R=N.zeros((nr))
precip_CTRL_f_d6_R=N.zeros((nr))
precip_CTRL_f_d7_R=N.zeros((nr))
precip_CTRL_f_d8_R=N.zeros((nr))
precip_CTRL_f_d9_R=N.zeros((nr))
precip_CTRL_f_d10_R=N.zeros((nr))
precip_CTRL_f_d1_R[0:37]=0.0
precip_CTRL_f_d2_R[0:37]=precip_CTRL_f[90-37:90]
precip_CTRL_f_d3_R[0:37]=precip_CTRL_f[182-37:182]
precip_CTRL_f_d4_R[0:37]=precip_CTRL_f[273-37:273]
precip_CTRL_f_d5_R[0:37]=precip_CTRL_f[363-37:363]
precip_CTRL_f_d6_R[0:37]=precip_CTRL_f[455-37:455]
precip_CTRL_f_d7_R[0:37]=precip_CTRL_f[547-37:547]
precip_CTRL_f_d8_R[0:37]=precip_CTRL_f[637-37:637]
precip_CTRL_f_d9_R[0:37]=precip_CTRL_f[728-37:728]
precip_CTRL_f_d10_R[0:37]=precip_CTRL_f[819-37:819]

precip_CTRL_f_d1_R[37:nr]=precip_CTRL_f_d1[0:92]
precip_CTRL_f_d2_R[37:nr]=precip_CTRL_f_d2[0:92]
precip_CTRL_f_d3_R[37:nr]=precip_CTRL_f_d3[0:92]
precip_CTRL_f_d4_R[37:nr]=precip_CTRL_f_d4[0:92]
precip_CTRL_f_d5_R[37:nr]=precip_CTRL_f_d5[0:92]
precip_CTRL_f_d6_R[37:nr]=precip_CTRL_f_d6[0:92]
precip_CTRL_f_d7_R[37:nr]=precip_CTRL_f_d7[0:92]
precip_CTRL_f_d8_R[37:nr]=precip_CTRL_f_d8[0:92]
precip_CTRL_f_d9_R[37:nr]=precip_CTRL_f_d9[0:92]
precip_CTRL_f_d10_R[37:nr]=precip_CTRL_f_d10[0:92]

cld_MF_cb_CTRL_f_d1_R=N.zeros((nr))
cld_MF_cb_CTRL_f_d2_R=N.zeros((nr))
cld_MF_cb_CTRL_f_d3_R=N.zeros((nr))
cld_MF_cb_CTRL_f_d4_R=N.zeros((nr))
cld_MF_cb_CTRL_f_d5_R=N.zeros((nr))
cld_MF_cb_CTRL_f_d6_R=N.zeros((nr))
cld_MF_cb_CTRL_f_d7_R=N.zeros((nr))
cld_MF_cb_CTRL_f_d8_R=N.zeros((nr))
cld_MF_cb_CTRL_f_d9_R=N.zeros((nr))
cld_MF_cb_CTRL_f_d10_R=N.zeros((nr))
cld_MF_cb_CTRL_f_d1_R[0:37]=0.0
cld_MF_cb_CTRL_f_d2_R[0:37]=cld_MF_cb_CTRL_f[90-37:90]
cld_MF_cb_CTRL_f_d3_R[0:37]=cld_MF_cb_CTRL_f[182-37:182]
cld_MF_cb_CTRL_f_d4_R[0:37]=cld_MF_cb_CTRL_f[273-37:273]
cld_MF_cb_CTRL_f_d5_R[0:37]=cld_MF_cb_CTRL_f[363-37:363]
cld_MF_cb_CTRL_f_d6_R[0:37]=cld_MF_cb_CTRL_f[455-37:455]
cld_MF_cb_CTRL_f_d7_R[0:37]=cld_MF_cb_CTRL_f[547-37:547]
cld_MF_cb_CTRL_f_d8_R[0:37]=cld_MF_cb_CTRL_f[637-37:637]
cld_MF_cb_CTRL_f_d9_R[0:37]=cld_MF_cb_CTRL_f[728-37:728]
cld_MF_cb_CTRL_f_d10_R[0:37]=cld_MF_cb_CTRL_f[819-37:819]

cld_MF_cb_CTRL_f_d1_R[37:nr]=cld_MF_cb_CTRL_f_d1[0:92]
cld_MF_cb_CTRL_f_d2_R[37:nr]=cld_MF_cb_CTRL_f_d2[0:92]
cld_MF_cb_CTRL_f_d3_R[37:nr]=cld_MF_cb_CTRL_f_d3[0:92]
cld_MF_cb_CTRL_f_d4_R[37:nr]=cld_MF_cb_CTRL_f_d4[0:92]
cld_MF_cb_CTRL_f_d5_R[37:nr]=cld_MF_cb_CTRL_f_d5[0:92]
cld_MF_cb_CTRL_f_d6_R[37:nr]=cld_MF_cb_CTRL_f_d6[0:92]
cld_MF_cb_CTRL_f_d7_R[37:nr]=cld_MF_cb_CTRL_f_d7[0:92]
cld_MF_cb_CTRL_f_d8_R[37:nr]=cld_MF_cb_CTRL_f_d8[0:92]
cld_MF_cb_CTRL_f_d9_R[37:nr]=cld_MF_cb_CTRL_f_d9[0:92]
cld_MF_cb_CTRL_f_d10_R[37:nr]=cld_MF_cb_CTRL_f_d10[0:92]

cld_frac_cb_CTRL_f_d1_R=N.zeros((nr))
cld_frac_cb_CTRL_f_d2_R=N.zeros((nr))
cld_frac_cb_CTRL_f_d3_R=N.zeros((nr))
cld_frac_cb_CTRL_f_d4_R=N.zeros((nr))
cld_frac_cb_CTRL_f_d5_R=N.zeros((nr))
cld_frac_cb_CTRL_f_d6_R=N.zeros((nr))
cld_frac_cb_CTRL_f_d7_R=N.zeros((nr))
cld_frac_cb_CTRL_f_d8_R=N.zeros((nr))
cld_frac_cb_CTRL_f_d9_R=N.zeros((nr))
cld_frac_cb_CTRL_f_d10_R=N.zeros((nr))
cld_frac_cb_CTRL_f_d1_R[0:37]=0.0
cld_frac_cb_CTRL_f_d2_R[0:37]=cld_frac_cb_CTRL_f[90-37:90]
cld_frac_cb_CTRL_f_d3_R[0:37]=cld_frac_cb_CTRL_f[182-37:182]
cld_frac_cb_CTRL_f_d4_R[0:37]=cld_frac_cb_CTRL_f[273-37:273]
cld_frac_cb_CTRL_f_d5_R[0:37]=cld_frac_cb_CTRL_f[363-37:363]
cld_frac_cb_CTRL_f_d6_R[0:37]=cld_frac_cb_CTRL_f[455-37:455]
cld_frac_cb_CTRL_f_d7_R[0:37]=cld_frac_cb_CTRL_f[547-37:547]
cld_frac_cb_CTRL_f_d8_R[0:37]=cld_frac_cb_CTRL_f[637-37:637]
cld_frac_cb_CTRL_f_d9_R[0:37]=cld_frac_cb_CTRL_f[728-37:728]
cld_frac_cb_CTRL_f_d10_R[0:37]=cld_frac_cb_CTRL_f[819-37:819]

cld_frac_cb_CTRL_f_d1_R[37:nr]=cld_frac_cb_CTRL_f_d1[0:92]
cld_frac_cb_CTRL_f_d2_R[37:nr]=cld_frac_cb_CTRL_f_d2[0:92]
cld_frac_cb_CTRL_f_d3_R[37:nr]=cld_frac_cb_CTRL_f_d3[0:92]
cld_frac_cb_CTRL_f_d4_R[37:nr]=cld_frac_cb_CTRL_f_d4[0:92]
cld_frac_cb_CTRL_f_d5_R[37:nr]=cld_frac_cb_CTRL_f_d5[0:92]
cld_frac_cb_CTRL_f_d6_R[37:nr]=cld_frac_cb_CTRL_f_d6[0:92]
cld_frac_cb_CTRL_f_d7_R[37:nr]=cld_frac_cb_CTRL_f_d7[0:92]
cld_frac_cb_CTRL_f_d8_R[37:nr]=cld_frac_cb_CTRL_f_d8[0:92]
cld_frac_cb_CTRL_f_d9_R[37:nr]=cld_frac_cb_CTRL_f_d9[0:92]
cld_frac_cb_CTRL_f_d10_R[37:nr]=cld_frac_cb_CTRL_f_d10[0:92]


cldBCu_MF_cb_CTRL_f_d1_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_f_d2_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_f_d3_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_f_d4_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_f_d5_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_f_d6_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_f_d7_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_f_d8_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_f_d9_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_f_d10_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_f_d1_R[0:37]=0.0
cldBCu_MF_cb_CTRL_f_d2_R[0:37]=cldBCu_MF_cb_CTRL_f[92-37:92]
cldBCu_MF_cb_CTRL_f_d3_R[0:37]=cldBCu_MF_cb_CTRL_f[182-37:182]
cldBCu_MF_cb_CTRL_f_d4_R[0:37]=cldBCu_MF_cb_CTRL_f[273-37:273]
cldBCu_MF_cb_CTRL_f_d5_R[0:37]=cldBCu_MF_cb_CTRL_f[363-37:363]
cldBCu_MF_cb_CTRL_f_d6_R[0:37]=cldBCu_MF_cb_CTRL_f[455-37:455]
cldBCu_MF_cb_CTRL_f_d7_R[0:37]=cldBCu_MF_cb_CTRL_f[547-37:547]
cldBCu_MF_cb_CTRL_f_d8_R[0:37]=cldBCu_MF_cb_CTRL_f[637-37:637]
cldBCu_MF_cb_CTRL_f_d9_R[0:37]=cldBCu_MF_cb_CTRL_f[728-37:728]
cldBCu_MF_cb_CTRL_f_d10_R[0:37]=cldBCu_MF_cb_CTRL_f[819-37:819]

cldBCu_MF_cb_CTRL_f_d1_R[37:nr]=cldBCu_MF_cb_CTRL_f_d1[0:92]
cldBCu_MF_cb_CTRL_f_d2_R[37:nr]=cldBCu_MF_cb_CTRL_f_d2[0:92]
cldBCu_MF_cb_CTRL_f_d3_R[37:nr]=cldBCu_MF_cb_CTRL_f_d3[0:92]
cldBCu_MF_cb_CTRL_f_d4_R[37:nr]=cldBCu_MF_cb_CTRL_f_d4[0:92]
cldBCu_MF_cb_CTRL_f_d5_R[37:nr]=cldBCu_MF_cb_CTRL_f_d5[0:92]
cldBCu_MF_cb_CTRL_f_d6_R[37:nr]=cldBCu_MF_cb_CTRL_f_d6[0:92]
cldBCu_MF_cb_CTRL_f_d7_R[37:nr]=cldBCu_MF_cb_CTRL_f_d7[0:92]
cldBCu_MF_cb_CTRL_f_d8_R[37:nr]=cldBCu_MF_cb_CTRL_f_d8[0:92]
cldBCu_MF_cb_CTRL_f_d9_R[37:nr]=cldBCu_MF_cb_CTRL_f_d9[0:92]
cldBCu_MF_cb_CTRL_f_d10_R[37:nr]=cldBCu_MF_cb_CTRL_f_d10[0:92]

cldBCu_frac_cb_CTRL_f_d1_R=N.zeros((nr))
cldBCu_frac_cb_CTRL_f_d2_R=N.zeros((nr))
cldBCu_frac_cb_CTRL_f_d3_R=N.zeros((nr))
cldBCu_frac_cb_CTRL_f_d4_R=N.zeros((nr))
cldBCu_frac_cb_CTRL_f_d5_R=N.zeros((nr))
cldBCu_frac_cb_CTRL_f_d6_R=N.zeros((nr))
cldBCu_frac_cb_CTRL_f_d7_R=N.zeros((nr))
cldBCu_frac_cb_CTRL_f_d8_R=N.zeros((nr))
cldBCu_frac_cb_CTRL_f_d9_R=N.zeros((nr))
cldBCu_frac_cb_CTRL_f_d10_R=N.zeros((nr))
cldBCu_frac_cb_CTRL_f_d1_R[0:37]=0.0
cldBCu_frac_cb_CTRL_f_d2_R[0:37]=cldBCu_frac_cb_CTRL_f[90-37:90]
cldBCu_frac_cb_CTRL_f_d3_R[0:37]=cldBCu_frac_cb_CTRL_f[182-37:182]
cldBCu_frac_cb_CTRL_f_d4_R[0:37]=cldBCu_frac_cb_CTRL_f[273-37:273]
cldBCu_frac_cb_CTRL_f_d5_R[0:37]=cldBCu_frac_cb_CTRL_f[363-37:363]
cldBCu_frac_cb_CTRL_f_d6_R[0:37]=cldBCu_frac_cb_CTRL_f[455-37:455]
cldBCu_frac_cb_CTRL_f_d7_R[0:37]=cldBCu_frac_cb_CTRL_f[547-37:547]
cldBCu_frac_cb_CTRL_f_d8_R[0:37]=cldBCu_frac_cb_CTRL_f[637-37:637]
cldBCu_frac_cb_CTRL_f_d9_R[0:37]=cldBCu_frac_cb_CTRL_f[728-37:728]
cldBCu_frac_cb_CTRL_f_d10_R[0:37]=cldBCu_frac_cb_CTRL_f[819-37:819]

cldBCu_frac_cb_CTRL_f_d1_R[37:nr]=cldBCu_frac_cb_CTRL_f_d1[0:92]
cldBCu_frac_cb_CTRL_f_d2_R[37:nr]=cldBCu_frac_cb_CTRL_f_d2[0:92]
cldBCu_frac_cb_CTRL_f_d3_R[37:nr]=cldBCu_frac_cb_CTRL_f_d3[0:92]
cldBCu_frac_cb_CTRL_f_d4_R[37:nr]=cldBCu_frac_cb_CTRL_f_d4[0:92]
cldBCu_frac_cb_CTRL_f_d5_R[37:nr]=cldBCu_frac_cb_CTRL_f_d5[0:92]
cldBCu_frac_cb_CTRL_f_d6_R[37:nr]=cldBCu_frac_cb_CTRL_f_d6[0:92]
cldBCu_frac_cb_CTRL_f_d7_R[37:nr]=cldBCu_frac_cb_CTRL_f_d7[0:92]
cldBCu_frac_cb_CTRL_f_d8_R[37:nr]=cldBCu_frac_cb_CTRL_f_d8[0:92]
cldBCu_frac_cb_CTRL_f_d9_R[37:nr]=cldBCu_frac_cb_CTRL_f_d9[0:92]
cldBCu_frac_cb_CTRL_f_d10_R[37:nr]=cldBCu_frac_cb_CTRL_f_d10[0:92]


n_10min=93


#
#mhalf
#

time_10min_MHALF_f_d1=N.zeros((n_10min))
time_10min_MHALF_f_d2=N.zeros((n_10min))
time_10min_MHALF_f_d3=N.zeros((n_10min))
time_10min_MHALF_f_d4=N.zeros((n_10min))
time_10min_MHALF_f_d5=N.zeros((n_10min))
time_10min_MHALF_f_d6=N.zeros((n_10min))
time_10min_MHALF_f_d7=N.zeros((n_10min))
time_10min_MHALF_f_d8=N.zeros((n_10min))
time_10min_MHALF_f_d9=N.zeros((n_10min))
time_10min_MHALF_f_d10=N.zeros((n_10min))

precip_MHALF_f_d1=N.zeros((n_10min))
precip_MHALF_f_d2=N.zeros((n_10min))
precip_MHALF_f_d3=N.zeros((n_10min))
precip_MHALF_f_d4=N.zeros((n_10min))
precip_MHALF_f_d5=N.zeros((n_10min))
precip_MHALF_f_d6=N.zeros((n_10min))
precip_MHALF_f_d7=N.zeros((n_10min))
precip_MHALF_f_d8=N.zeros((n_10min))
precip_MHALF_f_d9=N.zeros((n_10min))
precip_MHALF_f_d10=N.zeros((n_10min))

cld_MF_cb_MHALF_f_d1=N.zeros((n_10min))
cld_MF_cb_MHALF_f_d2=N.zeros((n_10min))
cld_MF_cb_MHALF_f_d3=N.zeros((n_10min))
cld_MF_cb_MHALF_f_d4=N.zeros((n_10min))
cld_MF_cb_MHALF_f_d5=N.zeros((n_10min))
cld_MF_cb_MHALF_f_d6=N.zeros((n_10min))
cld_MF_cb_MHALF_f_d7=N.zeros((n_10min))
cld_MF_cb_MHALF_f_d8=N.zeros((n_10min))
cld_MF_cb_MHALF_f_d9=N.zeros((n_10min))
cld_MF_cb_MHALF_f_d10=N.zeros((n_10min))

time_10min_MHALF_f_d1[0:92]=time_10min_MHALF_f[0:92]
time_10min_MHALF_f_d2[0:92]=time_10min_MHALF_f[90:182]-24.0
time_10min_MHALF_f_d3[0:92]=time_10min_MHALF_f[179:271]-48.0
time_10min_MHALF_f_d4[0:92]=time_10min_MHALF_f[269:361]-72.0
time_10min_MHALF_f_d5[0:92]=time_10min_MHALF_f[359:451]-96.0
time_10min_MHALF_f_d6[0:92]=time_10min_MHALF_f[449:541]-120.0
time_10min_MHALF_f_d7[0:92]=time_10min_MHALF_f[539:631]-144.0
time_10min_MHALF_f_d8[0:92]=time_10min_MHALF_f[629:721]-168.0
time_10min_MHALF_f_d9[0:92]=time_10min_MHALF_f[720:812]-192.0
time_10min_MHALF_f_d10[0:92]=time_10min_MHALF_f[811:903]-216.0


precip_MHALF_f_d1[0:92]=precip_MHALF_f[0:92]
precip_MHALF_f_d2[0:92]=precip_MHALF_f[90:182]
precip_MHALF_f_d3[0:92]=precip_MHALF_f[179:271]
precip_MHALF_f_d4[0:92]=precip_MHALF_f[269:361]
precip_MHALF_f_d5[0:92]=precip_MHALF_f[359:451]
precip_MHALF_f_d6[0:92]=precip_MHALF_f[449:541]
precip_MHALF_f_d7[0:92]=precip_MHALF_f[539:631]
precip_MHALF_f_d8[0:92]=precip_MHALF_f[629:721]
precip_MHALF_f_d9[0:92]=precip_MHALF_f[720:812]
precip_MHALF_f_d10[0:92]=precip_MHALF_f[811:903]

cld_MF_cb_MHALF_f_d1[0:92]=cld_MF_cb_MHALF_f[0:92]
cld_MF_cb_MHALF_f_d2[0:92]=cld_MF_cb_MHALF_f[90:182]
cld_MF_cb_MHALF_f_d3[0:92]=cld_MF_cb_MHALF_f[179:271]
cld_MF_cb_MHALF_f_d4[0:92]=cld_MF_cb_MHALF_f[269:361]
cld_MF_cb_MHALF_f_d5[0:92]=cld_MF_cb_MHALF_f[359:451]
cld_MF_cb_MHALF_f_d6[0:92]=cld_MF_cb_MHALF_f[449:541]
cld_MF_cb_MHALF_f_d7[0:92]=cld_MF_cb_MHALF_f[539:631]
cld_MF_cb_MHALF_f_d8[0:92]=cld_MF_cb_MHALF_f[629:721]
cld_MF_cb_MHALF_f_d9[0:92]=cld_MF_cb_MHALF_f[720:812]
cld_MF_cb_MHALF_f_d10[0:92]=cld_MF_cb_MHALF_f[811:903]



nr=36+n_10min
time_10min_MHALF_f_d1_R=N.zeros((nr))
time_10min_MHALF_f_d2_R=N.zeros((nr))
time_10min_MHALF_f_d3_R=N.zeros((nr))
time_10min_MHALF_f_d4_R=N.zeros((nr))
time_10min_MHALF_f_d5_R=N.zeros((nr))
time_10min_MHALF_f_d6_R=N.zeros((nr))
time_10min_MHALF_f_d7_R=N.zeros((nr))
time_10min_MHALF_f_d8_R=N.zeros((nr))
time_10min_MHALF_f_d9_R=N.zeros((nr))
time_10min_MHALF_f_d10_R=N.zeros((nr))

time_10min_MHALF_f_d1_R[0]=-9
for j in N.arange(36):
     time_10min_MHALF_f_d1_R[j+1]=time_10min_MHALF_f_d1_R[j]+0.25

time_10min_MHALF_f_d2_R[0:37]=time_10min_MHALF_f_d1_R[0:37]
time_10min_MHALF_f_d3_R[0:37]=time_10min_MHALF_f_d1_R[0:37]
time_10min_MHALF_f_d4_R[0:37]=time_10min_MHALF_f_d1_R[0:37]
time_10min_MHALF_f_d5_R[0:37]=time_10min_MHALF_f_d1_R[0:37]
time_10min_MHALF_f_d6_R[0:37]=time_10min_MHALF_f_d1_R[0:37]
time_10min_MHALF_f_d7_R[0:37]=time_10min_MHALF_f_d1_R[0:37]
time_10min_MHALF_f_d8_R[0:37]=time_10min_MHALF_f_d1_R[0:37]
time_10min_MHALF_f_d9_R[0:37]=time_10min_MHALF_f_d1_R[0:37]
time_10min_MHALF_f_d10_R[0:37]=time_10min_MHALF_f_d1_R[0:37]

time_10min_MHALF_f_d1_R[37:nr]=time_10min_MHALF_f_d1[0:92]
time_10min_MHALF_f_d2_R[37:nr]=time_10min_MHALF_f_d2[0:92]
time_10min_MHALF_f_d3_R[37:nr]=time_10min_MHALF_f_d3[0:92]
time_10min_MHALF_f_d4_R[37:nr]=time_10min_MHALF_f_d4[0:92]
time_10min_MHALF_f_d5_R[37:nr]=time_10min_MHALF_f_d5[0:92]
time_10min_MHALF_f_d6_R[37:nr]=time_10min_MHALF_f_d6[0:92]
time_10min_MHALF_f_d7_R[37:nr]=time_10min_MHALF_f_d7[0:92]
time_10min_MHALF_f_d8_R[37:nr]=time_10min_MHALF_f_d8[0:92]
time_10min_MHALF_f_d9_R[37:nr]=time_10min_MHALF_f_d9[0:92]
time_10min_MHALF_f_d10_R[37:nr]=time_10min_MHALF_f_d10[0:92]

precip_MHALF_f_d1_R=N.zeros((nr))
precip_MHALF_f_d2_R=N.zeros((nr))
precip_MHALF_f_d3_R=N.zeros((nr))
precip_MHALF_f_d4_R=N.zeros((nr))
precip_MHALF_f_d5_R=N.zeros((nr))
precip_MHALF_f_d6_R=N.zeros((nr))
precip_MHALF_f_d7_R=N.zeros((nr))
precip_MHALF_f_d8_R=N.zeros((nr))
precip_MHALF_f_d9_R=N.zeros((nr))
precip_MHALF_f_d10_R=N.zeros((nr))
precip_MHALF_f_d1_R[0:37]=0.0
precip_MHALF_f_d2_R[0:37]=precip_MHALF_f[90-37:90]
precip_MHALF_f_d3_R[0:37]=precip_MHALF_f[179-37:179]
precip_MHALF_f_d4_R[0:37]=precip_MHALF_f[269-37:269]
precip_MHALF_f_d5_R[0:37]=precip_MHALF_f[359-37:359]
precip_MHALF_f_d6_R[0:37]=precip_MHALF_f[449-37:449]
precip_MHALF_f_d7_R[0:37]=precip_MHALF_f[539-37:539]
precip_MHALF_f_d8_R[0:37]=precip_MHALF_f[629-37:629]
precip_MHALF_f_d9_R[0:37]=precip_MHALF_f[720-37:720]
precip_MHALF_f_d10_R[0:37]=precip_MHALF_f[811-37:811]

precip_MHALF_f_d1_R[37:nr]=precip_MHALF_f_d1[0:92]
precip_MHALF_f_d2_R[37:nr]=precip_MHALF_f_d2[0:92]
precip_MHALF_f_d3_R[37:nr]=precip_MHALF_f_d3[0:92]
precip_MHALF_f_d4_R[37:nr]=precip_MHALF_f_d4[0:92]
precip_MHALF_f_d5_R[37:nr]=precip_MHALF_f_d5[0:92]
precip_MHALF_f_d6_R[37:nr]=precip_MHALF_f_d6[0:92]
precip_MHALF_f_d7_R[37:nr]=precip_MHALF_f_d7[0:92]
precip_MHALF_f_d8_R[37:nr]=precip_MHALF_f_d8[0:92]
precip_MHALF_f_d9_R[37:nr]=precip_MHALF_f_d9[0:92]
precip_MHALF_f_d10_R[37:nr]=precip_MHALF_f_d10[0:92]

cld_MF_cb_MHALF_f_d1_R=N.zeros((nr))
cld_MF_cb_MHALF_f_d2_R=N.zeros((nr))
cld_MF_cb_MHALF_f_d3_R=N.zeros((nr))
cld_MF_cb_MHALF_f_d4_R=N.zeros((nr))
cld_MF_cb_MHALF_f_d5_R=N.zeros((nr))
cld_MF_cb_MHALF_f_d6_R=N.zeros((nr))
cld_MF_cb_MHALF_f_d7_R=N.zeros((nr))
cld_MF_cb_MHALF_f_d8_R=N.zeros((nr))
cld_MF_cb_MHALF_f_d9_R=N.zeros((nr))
cld_MF_cb_MHALF_f_d10_R=N.zeros((nr))
cld_MF_cb_MHALF_f_d1_R[0:37]=0.0
cld_MF_cb_MHALF_f_d2_R[0:37]=cld_MF_cb_MHALF_f[90-37:90]
cld_MF_cb_MHALF_f_d3_R[0:37]=cld_MF_cb_MHALF_f[179-37:179]
cld_MF_cb_MHALF_f_d4_R[0:37]=cld_MF_cb_MHALF_f[269-37:269]
cld_MF_cb_MHALF_f_d5_R[0:37]=cld_MF_cb_MHALF_f[359-37:359]
cld_MF_cb_MHALF_f_d6_R[0:37]=cld_MF_cb_MHALF_f[449-37:449]
cld_MF_cb_MHALF_f_d7_R[0:37]=cld_MF_cb_MHALF_f[539-37:539]
cld_MF_cb_MHALF_f_d8_R[0:37]=cld_MF_cb_MHALF_f[629-37:629]
cld_MF_cb_MHALF_f_d9_R[0:37]=cld_MF_cb_MHALF_f[720-37:720]
cld_MF_cb_MHALF_f_d10_R[0:37]=cld_MF_cb_MHALF_f[811-37:811]

cld_MF_cb_MHALF_f_d1_R[37:nr]=cld_MF_cb_MHALF_f_d1[0:92]
cld_MF_cb_MHALF_f_d2_R[37:nr]=cld_MF_cb_MHALF_f_d2[0:92]
cld_MF_cb_MHALF_f_d3_R[37:nr]=cld_MF_cb_MHALF_f_d3[0:92]
cld_MF_cb_MHALF_f_d4_R[37:nr]=cld_MF_cb_MHALF_f_d4[0:92]
cld_MF_cb_MHALF_f_d5_R[37:nr]=cld_MF_cb_MHALF_f_d5[0:92]
cld_MF_cb_MHALF_f_d6_R[37:nr]=cld_MF_cb_MHALF_f_d6[0:92]
cld_MF_cb_MHALF_f_d7_R[37:nr]=cld_MF_cb_MHALF_f_d7[0:92]
cld_MF_cb_MHALF_f_d8_R[37:nr]=cld_MF_cb_MHALF_f_d8[0:92]
cld_MF_cb_MHALF_f_d9_R[37:nr]=cld_MF_cb_MHALF_f_d9[0:92]
cld_MF_cb_MHALF_f_d10_R[37:nr]=cld_MF_cb_MHALF_f_d10[0:92]


time_10min_PHALF_f_d1=N.zeros((n_10min))
time_10min_PHALF_f_d2=N.zeros((n_10min))
time_10min_PHALF_f_d3=N.zeros((n_10min))
time_10min_PHALF_f_d4=N.zeros((n_10min))
time_10min_PHALF_f_d5=N.zeros((n_10min))
time_10min_PHALF_f_d6=N.zeros((n_10min))
time_10min_PHALF_f_d7=N.zeros((n_10min))
time_10min_PHALF_f_d8=N.zeros((n_10min))
time_10min_PHALF_f_d9=N.zeros((n_10min))
time_10min_PHALF_f_d10=N.zeros((n_10min))

precip_PHALF_f_d1=N.zeros((n_10min))
precip_PHALF_f_d2=N.zeros((n_10min))
precip_PHALF_f_d3=N.zeros((n_10min))
precip_PHALF_f_d4=N.zeros((n_10min))
precip_PHALF_f_d5=N.zeros((n_10min))
precip_PHALF_f_d6=N.zeros((n_10min))
precip_PHALF_f_d7=N.zeros((n_10min))
precip_PHALF_f_d8=N.zeros((n_10min))
precip_PHALF_f_d9=N.zeros((n_10min))
precip_PHALF_f_d10=N.zeros((n_10min))

cld_MF_cb_PHALF_f_d1=N.zeros((n_10min))
cld_MF_cb_PHALF_f_d2=N.zeros((n_10min))
cld_MF_cb_PHALF_f_d3=N.zeros((n_10min))
cld_MF_cb_PHALF_f_d4=N.zeros((n_10min))
cld_MF_cb_PHALF_f_d5=N.zeros((n_10min))
cld_MF_cb_PHALF_f_d6=N.zeros((n_10min))
cld_MF_cb_PHALF_f_d7=N.zeros((n_10min))
cld_MF_cb_PHALF_f_d8=N.zeros((n_10min))
cld_MF_cb_PHALF_f_d9=N.zeros((n_10min))
cld_MF_cb_PHALF_f_d10=N.zeros((n_10min))


time_10min_PHALF_f_d1[0:92]=time_10min_PHALF_f[0:92]
time_10min_PHALF_f_d2[0:92]=time_10min_PHALF_f[91:183]-24.0
time_10min_PHALF_f_d3[0:92]=time_10min_PHALF_f[183:275]-48.0
time_10min_PHALF_f_d4[0:92]=time_10min_PHALF_f[275:367]-72.0
time_10min_PHALF_f_d5[0:92]=time_10min_PHALF_f[367:459]-96.0
time_10min_PHALF_f_d6[0:92]=time_10min_PHALF_f[459:551]-120.0
time_10min_PHALF_f_d7[0:92]=time_10min_PHALF_f[550:642]-144.0
time_10min_PHALF_f_d8[0:92]=time_10min_PHALF_f[642:734]-168.0
time_10min_PHALF_f_d9[0:92]=time_10min_PHALF_f[733:825]-192.0
time_10min_PHALF_f_d10[0:92]=time_10min_PHALF_f[825:917]-216.0


precip_PHALF_f_d1[0:92]=precip_PHALF_f[0:92]
precip_PHALF_f_d2[0:92]=precip_PHALF_f[91:183]
precip_PHALF_f_d3[0:92]=precip_PHALF_f[183:275]
precip_PHALF_f_d4[0:92]=precip_PHALF_f[275:367]
precip_PHALF_f_d5[0:92]=precip_PHALF_f[367:459]
precip_PHALF_f_d6[0:92]=precip_PHALF_f[459:551]
precip_PHALF_f_d7[0:92]=precip_PHALF_f[550:642]
precip_PHALF_f_d8[0:92]=precip_PHALF_f[642:734]
precip_PHALF_f_d9[0:92]=precip_PHALF_f[733:825]
precip_PHALF_f_d10[0:92]=precip_PHALF_f[825:917]

cld_MF_cb_PHALF_f_d1[0:92]=cld_MF_cb_PHALF_f[0:92]
cld_MF_cb_PHALF_f_d2[0:92]=cld_MF_cb_PHALF_f[91:183]
cld_MF_cb_PHALF_f_d3[0:92]=cld_MF_cb_PHALF_f[183:275]
cld_MF_cb_PHALF_f_d4[0:92]=cld_MF_cb_PHALF_f[275:367]
cld_MF_cb_PHALF_f_d5[0:92]=cld_MF_cb_PHALF_f[367:459]
cld_MF_cb_PHALF_f_d6[0:92]=cld_MF_cb_PHALF_f[459:551]
cld_MF_cb_PHALF_f_d7[0:92]=cld_MF_cb_PHALF_f[550:642]
cld_MF_cb_PHALF_f_d8[0:92]=cld_MF_cb_PHALF_f[642:734]
cld_MF_cb_PHALF_f_d9[0:92]=cld_MF_cb_PHALF_f[733:825]
cld_MF_cb_PHALF_f_d10[0:92]=cld_MF_cb_PHALF_f[825:917]

nr=36+n_10min
time_10min_PHALF_f_d1_R=N.zeros((nr))
time_10min_PHALF_f_d2_R=N.zeros((nr))
time_10min_PHALF_f_d3_R=N.zeros((nr))
time_10min_PHALF_f_d4_R=N.zeros((nr))
time_10min_PHALF_f_d5_R=N.zeros((nr))
time_10min_PHALF_f_d6_R=N.zeros((nr))
time_10min_PHALF_f_d7_R=N.zeros((nr))
time_10min_PHALF_f_d8_R=N.zeros((nr))
time_10min_PHALF_f_d9_R=N.zeros((nr))
time_10min_PHALF_f_d10_R=N.zeros((nr))

time_10min_PHALF_f_d1_R[0]=-9
for j in N.arange(36):
     time_10min_PHALF_f_d1_R[j+1]=time_10min_PHALF_f_d1_R[j]+0.25

time_10min_PHALF_f_d2_R[0:37]=time_10min_PHALF_f_d1_R[0:37]
time_10min_PHALF_f_d3_R[0:37]=time_10min_PHALF_f_d1_R[0:37]
time_10min_PHALF_f_d4_R[0:37]=time_10min_PHALF_f_d1_R[0:37]
time_10min_PHALF_f_d5_R[0:37]=time_10min_PHALF_f_d1_R[0:37]
time_10min_PHALF_f_d6_R[0:37]=time_10min_PHALF_f_d1_R[0:37]
time_10min_PHALF_f_d7_R[0:37]=time_10min_PHALF_f_d1_R[0:37]
time_10min_PHALF_f_d8_R[0:37]=time_10min_PHALF_f_d1_R[0:37]
time_10min_PHALF_f_d9_R[0:37]=time_10min_PHALF_f_d1_R[0:37]
time_10min_PHALF_f_d10_R[0:37]=time_10min_PHALF_f_d1_R[0:37]

time_10min_PHALF_f_d1_R[37:nr]=time_10min_PHALF_f_d1[0:92]
time_10min_PHALF_f_d2_R[37:nr]=time_10min_PHALF_f_d2[0:92]
time_10min_PHALF_f_d3_R[37:nr]=time_10min_PHALF_f_d3[0:92]
time_10min_PHALF_f_d4_R[37:nr]=time_10min_PHALF_f_d4[0:92]
time_10min_PHALF_f_d5_R[37:nr]=time_10min_PHALF_f_d5[0:92]
time_10min_PHALF_f_d6_R[37:nr]=time_10min_PHALF_f_d6[0:92]
time_10min_PHALF_f_d7_R[37:nr]=time_10min_PHALF_f_d7[0:92]
time_10min_PHALF_f_d8_R[37:nr]=time_10min_PHALF_f_d8[0:92]
time_10min_PHALF_f_d9_R[37:nr]=time_10min_PHALF_f_d9[0:92]
time_10min_PHALF_f_d10_R[37:nr]=time_10min_PHALF_f_d10[0:92]

precip_PHALF_f_d1_R=N.zeros((nr))
precip_PHALF_f_d2_R=N.zeros((nr))
precip_PHALF_f_d3_R=N.zeros((nr))
precip_PHALF_f_d4_R=N.zeros((nr))
precip_PHALF_f_d5_R=N.zeros((nr))
precip_PHALF_f_d6_R=N.zeros((nr))
precip_PHALF_f_d7_R=N.zeros((nr))
precip_PHALF_f_d8_R=N.zeros((nr))
precip_PHALF_f_d9_R=N.zeros((nr))
precip_PHALF_f_d10_R=N.zeros((nr))
precip_PHALF_f_d1_R[0:37]=0.0
precip_PHALF_f_d2_R[0:37]=precip_PHALF_f[91-37:91]
precip_PHALF_f_d3_R[0:37]=precip_PHALF_f[183-37:183]
precip_PHALF_f_d4_R[0:37]=precip_PHALF_f[275-37:275]
precip_PHALF_f_d5_R[0:37]=precip_PHALF_f[367-37:367]
precip_PHALF_f_d6_R[0:37]=precip_PHALF_f[459-37:459]
precip_PHALF_f_d7_R[0:37]=precip_PHALF_f[550-37:550]
precip_PHALF_f_d8_R[0:37]=precip_PHALF_f[642-37:642]
precip_PHALF_f_d9_R[0:37]=precip_PHALF_f[733-37:733]
precip_PHALF_f_d10_R[0:37]=precip_PHALF_f[825-37:825]

precip_PHALF_f_d1_R[37:nr]=precip_PHALF_f_d1[0:92]
precip_PHALF_f_d2_R[37:nr]=precip_PHALF_f_d2[0:92]
precip_PHALF_f_d3_R[37:nr]=precip_PHALF_f_d3[0:92]
precip_PHALF_f_d4_R[37:nr]=precip_PHALF_f_d4[0:92]
precip_PHALF_f_d5_R[37:nr]=precip_PHALF_f_d5[0:92]
precip_PHALF_f_d6_R[37:nr]=precip_PHALF_f_d6[0:92]
precip_PHALF_f_d7_R[37:nr]=precip_PHALF_f_d7[0:92]
precip_PHALF_f_d8_R[37:nr]=precip_PHALF_f_d8[0:92]
precip_PHALF_f_d9_R[37:nr]=precip_PHALF_f_d9[0:92]
precip_PHALF_f_d10_R[37:nr]=precip_PHALF_f_d10[0:92]

cld_MF_cb_PHALF_f_d1_R=N.zeros((nr))
cld_MF_cb_PHALF_f_d2_R=N.zeros((nr))
cld_MF_cb_PHALF_f_d3_R=N.zeros((nr))
cld_MF_cb_PHALF_f_d4_R=N.zeros((nr))
cld_MF_cb_PHALF_f_d5_R=N.zeros((nr))
cld_MF_cb_PHALF_f_d6_R=N.zeros((nr))
cld_MF_cb_PHALF_f_d7_R=N.zeros((nr))
cld_MF_cb_PHALF_f_d8_R=N.zeros((nr))
cld_MF_cb_PHALF_f_d9_R=N.zeros((nr))
cld_MF_cb_PHALF_f_d10_R=N.zeros((nr))
cld_MF_cb_PHALF_f_d1_R[0:37]=0.0
cld_MF_cb_PHALF_f_d2_R[0:37]=cld_MF_cb_PHALF_f[91-37:91]
cld_MF_cb_PHALF_f_d3_R[0:37]=cld_MF_cb_PHALF_f[183-37:183]
cld_MF_cb_PHALF_f_d4_R[0:37]=cld_MF_cb_PHALF_f[275-37:275]
cld_MF_cb_PHALF_f_d5_R[0:37]=cld_MF_cb_PHALF_f[367-37:367]
cld_MF_cb_PHALF_f_d6_R[0:37]=cld_MF_cb_PHALF_f[459-37:459]
cld_MF_cb_PHALF_f_d7_R[0:37]=cld_MF_cb_PHALF_f[550-37:550]
cld_MF_cb_PHALF_f_d8_R[0:37]=cld_MF_cb_PHALF_f[642-37:642]
cld_MF_cb_PHALF_f_d9_R[0:37]=cld_MF_cb_PHALF_f[733-37:733]
cld_MF_cb_PHALF_f_d10_R[0:37]=cld_MF_cb_PHALF_f[825-37:825]

cld_MF_cb_PHALF_f_d1_R[37:nr]=cld_MF_cb_PHALF_f_d1[0:92]
cld_MF_cb_PHALF_f_d2_R[37:nr]=cld_MF_cb_PHALF_f_d2[0:92]
cld_MF_cb_PHALF_f_d3_R[37:nr]=cld_MF_cb_PHALF_f_d3[0:92]
cld_MF_cb_PHALF_f_d4_R[37:nr]=cld_MF_cb_PHALF_f_d4[0:92]
cld_MF_cb_PHALF_f_d5_R[37:nr]=cld_MF_cb_PHALF_f_d5[0:92]
cld_MF_cb_PHALF_f_d6_R[37:nr]=cld_MF_cb_PHALF_f_d6[0:92]
cld_MF_cb_PHALF_f_d7_R[37:nr]=cld_MF_cb_PHALF_f_d7[0:92]
cld_MF_cb_PHALF_f_d8_R[37:nr]=cld_MF_cb_PHALF_f_d8[0:92]
cld_MF_cb_PHALF_f_d9_R[37:nr]=cld_MF_cb_PHALF_f_d9[0:92]
cld_MF_cb_PHALF_f_d10_R[37:nr]=cld_MF_cb_PHALF_f_d10[0:92]

#
#homogenization
#

time_10min_RTHQVD2_d2=N.zeros((n_10min))
precip_RTHQVD2_d2=N.zeros((n_10min))
cld_MF_cb_RTHQVD2_d2=N.zeros((n_10min))
time_10min_RTHQVD2_d2[0:92]=time_10min_RTHQVD2[88:180]-24.0
precip_RTHQVD2_d2[0:92]=precip_RTHQVD2[88:180]
cld_MF_cb_RTHQVD2_d2[0:92]=cld_MF_cb_RTHQVD2[88:180]
time_10min_RTHQVD2_d2_R=N.zeros((nr))
time_10min_RTHQVD2_d2_R[0]=-9
for j in N.arange(36):
     time_10min_RTHQVD2_d2_R[j+1]=time_10min_RTHQVD2_d2_R[j]+0.25
time_10min_RTHQVD2_d2_R[37:nr]=time_10min_RTHQVD2_d2[0:92]
precip_RTHQVD2_d2_R=N.zeros((nr))
cld_MF_cb_RTHQVD2_d2_R=N.zeros((nr))
precip_RTHQVD2_d2_R[0:37]=precip_RTHQVD2[88-37:88]
precip_RTHQVD2_d2_R[37:nr]=precip_RTHQVD2_d2[0:92]
cld_MF_cb_RTHQVD2_d2_R[0:37]=cld_MF_cb_RTHQVD2[88-37:88]
cld_MF_cb_RTHQVD2_d2_R[37:nr]=cld_MF_cb_RTHQVD2_d2[0:92]


time_10min_RTHQVD3_d3=N.zeros((n_10min))
precip_RTHQVD3_d3=N.zeros((n_10min))
cld_MF_cb_RTHQVD3_d3=N.zeros((n_10min))
time_10min_RTHQVD3_d3[0:92]=time_10min_RTHQVD3[181:273]-48.0
precip_RTHQVD3_d3[0:92]=precip_RTHQVD3[181:273]
cld_MF_cb_RTHQVD3_d3[0:92]=cld_MF_cb_RTHQVD3[181:273]
time_10min_RTHQVD3_d3_R=N.zeros((nr))
time_10min_RTHQVD3_d3_R[0]=-9
for j in N.arange(36):
     time_10min_RTHQVD3_d3_R[j+1]=time_10min_RTHQVD3_d3_R[j]+0.25
time_10min_RTHQVD3_d3_R[37:nr]=time_10min_RTHQVD3_d3[0:92]
precip_RTHQVD3_d3_R=N.zeros((nr))
cld_MF_cb_RTHQVD3_d3_R=N.zeros((nr))
precip_RTHQVD3_d3_R[0:37]=precip_RTHQVD3[181-37:181]
precip_RTHQVD3_d3_R[37:nr]=precip_RTHQVD3_d3[0:92]
cld_MF_cb_RTHQVD3_d3_R[0:37]=cld_MF_cb_RTHQVD3[181-37:181]
cld_MF_cb_RTHQVD3_d3_R[37:nr]=cld_MF_cb_RTHQVD3_d3[0:92]

time_10min_RTHQVD4_d4=N.zeros((n_10min))
precip_RTHQVD4_d4=N.zeros((n_10min))
cld_MF_cb_RTHQVD4_d4=N.zeros((n_10min))
time_10min_RTHQVD4_d4[0:92]=time_10min_RTHQVD4[273:365]-72.0
precip_RTHQVD4_d4[0:92]=precip_RTHQVD4[273:365]
cld_MF_cb_RTHQVD4_d4[0:92]=cld_MF_cb_RTHQVD4[273:365]
time_10min_RTHQVD4_d4_R=N.zeros((nr))
time_10min_RTHQVD4_d4_R[0]=-9
for j in N.arange(36):
     time_10min_RTHQVD4_d4_R[j+1]=time_10min_RTHQVD4_d4_R[j]+0.25
time_10min_RTHQVD4_d4_R[37:nr]=time_10min_RTHQVD4_d4[0:92]
precip_RTHQVD4_d4_R=N.zeros((nr))
cld_MF_cb_RTHQVD4_d4_R=N.zeros((nr))
precip_RTHQVD4_d4_R[0:37]=precip_RTHQVD4[273-37:273]
precip_RTHQVD4_d4_R[37:nr]=precip_RTHQVD4_d4[0:92]
cld_MF_cb_RTHQVD4_d4_R[0:37]=cld_MF_cb_RTHQVD4[273-37:273]
cld_MF_cb_RTHQVD4_d4_R[37:nr]=cld_MF_cb_RTHQVD4_d4[0:92]

time_10min_RTHQVD5_d5=N.zeros((n_10min))
precip_RTHQVD5_d5=N.zeros((n_10min))
cld_MF_cb_RTHQVD5_d5=N.zeros((n_10min))
time_10min_RTHQVD5_d5[0:92]=time_10min_RTHQVD5[362:454]-96.0
precip_RTHQVD5_d5[0:92]=precip_RTHQVD5[362:454]
cld_MF_cb_RTHQVD5_d5[0:92]=cld_MF_cb_RTHQVD5[362:454]
time_10min_RTHQVD5_d5_R=N.zeros((nr))
time_10min_RTHQVD5_d5_R[0]=-9
for j in N.arange(36):
     time_10min_RTHQVD5_d5_R[j+1]=time_10min_RTHQVD5_d5_R[j]+0.25
time_10min_RTHQVD5_d5_R[37:nr]=time_10min_RTHQVD5_d5[0:92]
precip_RTHQVD5_d5_R=N.zeros((nr))
cld_MF_cb_RTHQVD5_d5_R=N.zeros((nr))
precip_RTHQVD5_d5_R[0:37]=precip_RTHQVD5[362-37:362]
precip_RTHQVD5_d5_R[37:nr]=precip_RTHQVD5_d5[0:92]
cld_MF_cb_RTHQVD5_d5_R[0:37]=cld_MF_cb_RTHQVD5[362-37:362]
cld_MF_cb_RTHQVD5_d5_R[37:nr]=cld_MF_cb_RTHQVD5_d5[0:92]

time_10min_RTHQVD6_d6=N.zeros((n_10min))
precip_RTHQVD6_d6=N.zeros((n_10min))
cld_MF_cb_RTHQVD6_d6=N.zeros((n_10min))
time_10min_RTHQVD6_d6[0:92]=time_10min_RTHQVD6[453:545]-120.0
precip_RTHQVD6_d6[0:92]=precip_RTHQVD6[453:545]
cld_MF_cb_RTHQVD6_d6[0:92]=cld_MF_cb_RTHQVD6[453:545]
time_10min_RTHQVD6_d6_R=N.zeros((nr))
time_10min_RTHQVD6_d6_R[0]=-9
for j in N.arange(36):
     time_10min_RTHQVD6_d6_R[j+1]=time_10min_RTHQVD6_d6_R[j]+0.25
time_10min_RTHQVD6_d6_R[37:nr]=time_10min_RTHQVD6_d6[0:92]
precip_RTHQVD6_d6_R=N.zeros((nr))
cld_MF_cb_RTHQVD6_d6_R=N.zeros((nr))
precip_RTHQVD6_d6_R[0:37]=precip_RTHQVD6[453-37:453]
precip_RTHQVD6_d6_R[37:nr]=precip_RTHQVD6_d6[0:92]
cld_MF_cb_RTHQVD6_d6_R[0:37]=cld_MF_cb_RTHQVD6[453-37:453]
cld_MF_cb_RTHQVD6_d6_R[37:nr]=cld_MF_cb_RTHQVD6_d6[0:92]


time_10min_RTHQVD7_d7=N.zeros((n_10min))
precip_RTHQVD7_d7=N.zeros((n_10min))
cld_MF_cb_RTHQVD7_d7=N.zeros((n_10min))
cld_top_max_RTHQVD7_d7=N.zeros((n_10min))
time_10min_RTHQVD7_d7[0:92]=time_10min_RTHQVD7[544:636]-144.0
precip_RTHQVD7_d7[0:92]=precip_RTHQVD7[544:636]
cld_MF_cb_RTHQVD7_d7[0:92]=cld_MF_cb_RTHQVD7[544:636]
time_10min_RTHQVD7_d7_R=N.zeros((nr))
time_10min_RTHQVD7_d7_R[0]=-9
for j in N.arange(36):
     time_10min_RTHQVD7_d7_R[j+1]=time_10min_RTHQVD7_d7_R[j]+0.25
time_10min_RTHQVD7_d7_R[37:nr]=time_10min_RTHQVD7_d7[0:92]
precip_RTHQVD7_d7_R=N.zeros((nr))
cld_MF_cb_RTHQVD7_d7_R=N.zeros((nr))
precip_RTHQVD7_d7_R[0:37]=precip_RTHQVD7[544-37:544]
precip_RTHQVD7_d7_R[37:nr]=precip_RTHQVD7_d7[0:92]
cld_MF_cb_RTHQVD7_d7_R[0:37]=cld_MF_cb_RTHQVD7[544-37:544]
cld_MF_cb_RTHQVD7_d7_R[37:nr]=cld_MF_cb_RTHQVD7_d7[0:92]

time_10min_RTHQVD8_d8=N.zeros((n_10min))
precip_RTHQVD8_d8=N.zeros((n_10min))
cld_MF_cb_RTHQVD8_d8=N.zeros((n_10min))
time_10min_RTHQVD8_d8[0:92]=time_10min_RTHQVD8[638:730]-168.0
precip_RTHQVD8_d8[0:92]=precip_RTHQVD8[638:730]
cld_MF_cb_RTHQVD8_d8[0:92]=cld_MF_cb_RTHQVD8[638:730]
time_10min_RTHQVD8_d8_R=N.zeros((nr))
time_10min_RTHQVD8_d8_R[0]=-9
for j in N.arange(36):
     time_10min_RTHQVD8_d8_R[j+1]=time_10min_RTHQVD8_d8_R[j]+0.25
time_10min_RTHQVD8_d8_R[37:nr]=time_10min_RTHQVD8_d8[0:92]
precip_RTHQVD8_d8_R=N.zeros((nr))
cld_MF_cb_RTHQVD8_d8_R=N.zeros((nr))
precip_RTHQVD8_d8_R[0:37]=precip_RTHQVD8[638-37:638]
precip_RTHQVD8_d8_R[37:nr]=precip_RTHQVD8_d8[0:92]
cld_MF_cb_RTHQVD8_d8_R[0:37]=cld_MF_cb_RTHQVD8[638-37:638]
cld_MF_cb_RTHQVD8_d8_R[37:nr]=cld_MF_cb_RTHQVD8_d8[0:92]

time_10min_RTHQVD9_d9=N.zeros((n_10min))
precip_RTHQVD9_d9=N.zeros((n_10min))
cld_MF_cb_RTHQVD9_d9=N.zeros((n_10min))
time_10min_RTHQVD9_d9[0:92]=time_10min_RTHQVD9[727:819]-192.0
precip_RTHQVD9_d9[0:92]=precip_RTHQVD9[727:819]
cld_MF_cb_RTHQVD9_d9[0:92]=cld_MF_cb_RTHQVD9[727:819]
time_10min_RTHQVD9_d9_R=N.zeros((nr))
time_10min_RTHQVD9_d9_R[0]=-9
for j in N.arange(36):
     time_10min_RTHQVD9_d9_R[j+1]=time_10min_RTHQVD9_d9_R[j]+0.25
time_10min_RTHQVD9_d9_R[37:nr]=time_10min_RTHQVD9_d9[0:92]
precip_RTHQVD9_d9_R=N.zeros((nr))
cld_MF_cb_RTHQVD9_d9_R=N.zeros((nr))
precip_RTHQVD9_d9_R[0:37]=precip_RTHQVD9[727-37:727]
precip_RTHQVD9_d9_R[37:nr]=precip_RTHQVD9_d9[0:92]
cld_MF_cb_RTHQVD9_d9_R[0:37]=cld_MF_cb_RTHQVD9[727-37:727]
cld_MF_cb_RTHQVD9_d9_R[37:nr]=cld_MF_cb_RTHQVD9_d9[0:92]

time_10min_RTHQVD10_d10=N.zeros((n_10min))
precip_RTHQVD10_d10=N.zeros((n_10min))
cld_MF_cb_RTHQVD10_d10=N.zeros((n_10min))
time_10min_RTHQVD10_d10[0:92]=time_10min_RTHQVD10[819:911]-216.0
precip_RTHQVD10_d10[0:92]=precip_RTHQVD10[819:911]
cld_MF_cb_RTHQVD10_d10[0:92]=cld_MF_cb_RTHQVD10[819:911]
time_10min_RTHQVD10_d10_R=N.zeros((nr))
time_10min_RTHQVD10_d10_R[0]=-9
for j in N.arange(36):
     time_10min_RTHQVD10_d10_R[j+1]=time_10min_RTHQVD10_d10_R[j]+0.25
time_10min_RTHQVD10_d10_R[37:nr]=time_10min_RTHQVD10_d10[0:92]
precip_RTHQVD10_d10_R=N.zeros((nr))
cld_MF_cb_RTHQVD10_d10_R=N.zeros((nr))
precip_RTHQVD10_d10_R[0:37]=precip_RTHQVD10[819-37:819]
precip_RTHQVD10_d10_R[37:nr]=precip_RTHQVD10_d10[0:92]
cld_MF_cb_RTHQVD10_d10_R[0:37]=cld_MF_cb_RTHQVD10[819-37:819]
cld_MF_cb_RTHQVD10_d10_R[37:nr]=cld_MF_cb_RTHQVD10_d10[0:92]

time_10min_CTRL_f_ens_R=N.zeros((11,nr))
time_10min_CTRL_f_ens_mean_R=N.zeros((nr))
time_10min_CTRL_f_ens_stdv_R=N.zeros((nr))
time_10min_CTRL_f_ens_meanP_R=N.zeros((nr))
time_10min_CTRL_f_ens_meanM_R=N.zeros((nr))

precip_CTRL_f_ens_R=N.zeros((11,nr))
precip_CTRL_f_ens_mean_R=N.zeros((nr))
precip_CTRL_f_ens_stdv_R=N.zeros((nr))
precip_CTRL_f_ens_meanP_R=N.zeros((nr))
precip_CTRL_f_ens_meanM_R=N.zeros((nr))
precip_RTHQV_ens_RC=N.zeros((11,nr))

precip_PHALF_f_ens_RC=N.zeros((11,nr))
precip_MHALF_f_ens_RC=N.zeros((11,nr))
precip_CTRL_f_ens_RC=N.zeros((11,nr))
precip_PHALF_f_ens_mean_RC=N.zeros((nr))
precip_PHALF_f_ens_meanP_RC=N.zeros((nr))
precip_PHALF_f_ens_meanM_RC=N.zeros((nr))
precip_PHALF_f_ens_stdv_RC=N.zeros((nr))

precip_MHALF_f_ens_mean_RC=N.zeros((nr))
precip_MHALF_f_ens_meanP_RC=N.zeros((nr))
precip_MHALF_f_ens_meanM_RC=N.zeros((nr))
precip_MHALF_f_ens_stdv_RC=N.zeros((nr))

precip_RTHQV_ens_mean_RC=N.zeros((nr))
precip_RTHQV_ens_meanP_RC=N.zeros((nr))
precip_RTHQV_ens_meanM_RC=N.zeros((nr))
precip_RTHQV_ens_stdv_RC=N.zeros((nr))

precip_CTRL_f_ens_mean_RC=N.zeros((nr))
precip_CTRL_f_ens_meanP_RC=N.zeros((nr))
precip_CTRL_f_ens_meanM_RC=N.zeros((nr))
precip_CTRL_f_ens_stdv_RC=N.zeros((nr))

cld_MF_cb_CTRL_f_ens_R=N.zeros((11,nr))
cld_MF_cb_CTRL_f_ens_mean_R=N.zeros((nr))
cld_MF_cb_CTRL_f_ens_stdv_R=N.zeros((nr))
cld_MF_cb_CTRL_f_ens_meanP_R=N.zeros((nr))
cld_MF_cb_CTRL_f_ens_meanM_R=N.zeros((nr))

cld_frac_cb_CTRL_f_ens_R=N.zeros((11,nr))
cld_frac_cb_CTRL_f_ens_mean_R=N.zeros((nr))
cld_frac_cb_CTRL_f_ens_stdv_R=N.zeros((nr))
cld_frac_cb_CTRL_f_ens_meanP_R=N.zeros((nr))
cld_frac_cb_CTRL_f_ens_meanM_R=N.zeros((nr))

cld_MFpc_cb_CTRL_f_ens_R=N.zeros((11,nr))
cld_MFpc_cb_CTRL_f_ens_mean_R=N.zeros((nr))
cld_MFpc_cb_CTRL_f_ens_stdv_R=N.zeros((nr))
cld_MFpc_cb_CTRL_f_ens_meanP_R=N.zeros((nr))
cld_MFpc_cb_CTRL_f_ens_meanM_R=N.zeros((nr))

cldBCu_MF_cb_CTRL_f_ens_R=N.zeros((11,nr))
cldBCu_MF_cb_CTRL_f_ens_mean_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_f_ens_stdv_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_f_ens_meanP_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_f_ens_meanM_R=N.zeros((nr))

cldBCu_frac_cb_CTRL_f_ens_R=N.zeros((11,nr))
cldBCu_frac_cb_CTRL_f_ens_mean_R=N.zeros((nr))
cldBCu_frac_cb_CTRL_f_ens_stdv_R=N.zeros((nr))
cldBCu_frac_cb_CTRL_f_ens_meanP_R=N.zeros((nr))
cldBCu_frac_cb_CTRL_f_ens_meanM_R=N.zeros((nr))

time_10min_CTRL_ens_R=N.zeros((11,nr))
time_10min_CTRL_ens_mean_R=N.zeros((nr))
time_10min_CTRL_ens_stdv_R=N.zeros((nr))
time_10min_CTRL_ens_meanP_R=N.zeros((nr))
time_10min_CTRL_ens_meanM_R=N.zeros((nr))

precip_CTRL_ens_R=N.zeros((11,nr))
precip_CTRL_ens_mean_R=N.zeros((nr))
precip_CTRL_ens_stdv_R=N.zeros((nr))
precip_CTRL_ens_meanP_R=N.zeros((nr))
precip_CTRL_ens_meanM_R=N.zeros((nr))

cld_MF_cb_CTRL_ens_R=N.zeros((11,nr))
cld_MF_cb_CTRL_ens_mean_R=N.zeros((nr))
cld_MF_cb_CTRL_ens_stdv_R=N.zeros((nr))
cld_MF_cb_CTRL_ens_meanP_R=N.zeros((nr))
cld_MF_cb_CTRL_ens_meanM_R=N.zeros((nr))

cld_MFpc_cb_CTRL_ens_R=N.zeros((11,nr))
cld_MFpc_cb_CTRL_ens_mean_R=N.zeros((nr))
cld_MFpc_cb_CTRL_ens_stdv_R=N.zeros((nr))
cld_MFpc_cb_CTRL_ens_meanP_R=N.zeros((nr))
cld_MFpc_cb_CTRL_ens_meanM_R=N.zeros((nr))

cldBCu_MF_cb_CTRL_ens_R=N.zeros((11,nr))
cldBCu_MF_cb_CTRL_ens_mean_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_ens_stdv_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_ens_meanP_R=N.zeros((nr))
cldBCu_MF_cb_CTRL_ens_meanM_R=N.zeros((nr))

time_10min_PHALF_f_ens_R=N.zeros((11,nr))
time_10min_PHALF_f_ens_mean_R=N.zeros((nr))
time_10min_PHALF_f_ens_stdv_R=N.zeros((nr))
time_10min_PHALF_f_ens_meanP_R=N.zeros((nr))
time_10min_PHALF_f_ens_meanM_R=N.zeros((nr))

precip_PHALF_f_ens_R=N.zeros((11,nr))
precip_PHALF_f_ens_mean_R=N.zeros((nr))
precip_PHALF_f_ens_stdv_R=N.zeros((nr))
precip_PHALF_f_ens_meanP_R=N.zeros((nr))
precip_PHALF_f_ens_meanM_R=N.zeros((nr))

cld_MF_cb_PHALF_f_ens_R=N.zeros((11,nr))
cld_MF_cb_PHALF_f_ens_mean_R=N.zeros((nr))
cld_MF_cb_PHALF_f_ens_stdv_R=N.zeros((nr))
cld_MF_cb_PHALF_f_ens_meanP_R=N.zeros((nr))
cld_MF_cb_PHALF_f_ens_meanM_R=N.zeros((nr))

cld_MFpc_cb_PHALF_f_ens_R=N.zeros((11,nr))
cld_MFpc_cb_PHALF_f_ens_mean_R=N.zeros((nr))
cld_MFpc_cb_PHALF_f_ens_stdv_R=N.zeros((nr))
cld_MFpc_cb_PHALF_f_ens_meanP_R=N.zeros((nr))
cld_MFpc_cb_PHALF_f_ens_meanM_R=N.zeros((nr))

time_10min_MHALF_f_ens_R=N.zeros((11,nr))
time_10min_MHALF_f_ens_mean_R=N.zeros((nr))
time_10min_MHALF_f_ens_stdv_R=N.zeros((nr))
time_10min_MHALF_f_ens_meanP_R=N.zeros((nr))
time_10min_MHALF_f_ens_meanM_R=N.zeros((nr))

precip_MHALF_f_ens_R=N.zeros((11,nr))
precip_MHALF_f_ens_mean_R=N.zeros((nr))
precip_MHALF_f_ens_stdv_R=N.zeros((nr))
precip_MHALF_f_ens_meanP_R=N.zeros((nr))
precip_MHALF_f_ens_meanM_R=N.zeros((nr))

cld_MF_cb_MHALF_f_ens_R=N.zeros((11,nr))
cld_MF_cb_MHALF_f_ens_mean_R=N.zeros((nr))
cld_MF_cb_MHALF_f_ens_stdv_R=N.zeros((nr))
cld_MF_cb_MHALF_f_ens_meanP_R=N.zeros((nr))
cld_MF_cb_MHALF_f_ens_meanM_R=N.zeros((nr))

cld_MFpc_cb_MHALF_f_ens_R=N.zeros((11,nr))
cld_MFpc_cb_MHALF_f_ens_mean_R=N.zeros((nr))
cld_MFpc_cb_MHALF_f_ens_stdv_R=N.zeros((nr))
cld_MFpc_cb_MHALF_f_ens_meanP_R=N.zeros((nr))
cld_MFpc_cb_MHALF_f_ens_meanM_R=N.zeros((nr))

cld_MFpc_cb_CTRL_ens_mean_R=N.zeros((nr))
cld_MFpc_cb_PHALF_f_ens_mean_R=N.zeros((nr))
cld_MFpc_cb_MHALF_f_ens_mean_R=N.zeros((nr))

time_10min_RTHQV_ens_R=N.zeros((11,nr))
time_10min_RTHQV_ens_mean_R=N.zeros((nr))
time_10min_RTHQV_ens_stdv_R=N.zeros((nr))
time_10min_RTHQV_ens_meanP_R=N.zeros((nr))
time_10min_RTHQV_ens_meanM_R=N.zeros((nr))

precip_RTHQV_ens_R=N.zeros((11,nr))
precip_RTHQV_ens_mean_R=N.zeros((nr))
precip_RTHQV_ens_stdv_R=N.zeros((nr))
precip_RTHQV_ens_meanP_R=N.zeros((nr))
precip_RTHQV_ens_meanM_R=N.zeros((nr))

cld_MF_cb_RTHQV_ens_R=N.zeros((11,nr))
cld_MF_cb_RTHQV_ens_mean_R=N.zeros((nr))
cld_MF_cb_RTHQV_ens_stdv_R=N.zeros((nr))
cld_MF_cb_RTHQV_ens_meanP_R=N.zeros((nr))
cld_MF_cb_RTHQV_ens_meanM_R=N.zeros((nr))

precip_RTHQV_ens_R=N.zeros((11,nr))
precip_RTHQV_ens_mean_R=N.zeros((nr))
precip_RTHQV_ens_stdv_R=N.zeros((nr))
precip_RTHQV_ens_meanP_R=N.zeros((nr))
precip_RTHQV_ens_meanM_R=N.zeros((nr))

cld_MF_cb_RTHQV_ens_R=N.zeros((11,nr))
cld_MF_cb_RTHQV_ens_mean_R=N.zeros((nr))
cld_MF_cb_RTHQV_ens_stdv_R=N.zeros((nr))
cld_MF_cb_RTHQV_ens_meanP_R=N.zeros((nr))
cld_MF_cb_RTHQV_ens_meanM_R=N.zeros((nr))

aa=N.zeros((5))
aa[0]=1
aa[1]=2
aa[2]=3
aa[3]=4
aa[4]=5

#print aa[:].mean()
#print aa[0:5].mean()
#print aa[0:4].mean()
#print aa[1:5].mean()

for j in N.arange(nr):
     time_10min_CTRL_f_ens_R[0,j]=0
     time_10min_CTRL_f_ens_R[1,j]=time_10min_CTRL_f_d1_R[j]
     time_10min_CTRL_f_ens_R[2,j]=time_10min_CTRL_f_d2_R[j]
     time_10min_CTRL_f_ens_R[3,j]=time_10min_CTRL_f_d3_R[j]
     time_10min_CTRL_f_ens_R[4,j]=time_10min_CTRL_f_d4_R[j]
     time_10min_CTRL_f_ens_R[5,j]=time_10min_CTRL_f_d5_R[j]
     time_10min_CTRL_f_ens_R[6,j]=time_10min_CTRL_f_d6_R[j]
     time_10min_CTRL_f_ens_R[7,j]=time_10min_CTRL_f_d7_R[j]
     time_10min_CTRL_f_ens_R[8,j]=time_10min_CTRL_f_d8_R[j]
     time_10min_CTRL_f_ens_R[9,j]=time_10min_CTRL_f_d9_R[j]
     time_10min_CTRL_f_ens_R[10,j]=time_10min_CTRL_f_d10_R[j]
  
     time_10min_CTRL_f_ens_mean_R[j]=time_10min_CTRL_f_ens_R[2:11,j].mean()
     time_10min_CTRL_f_ens_stdv_R[j]=np.std(time_10min_CTRL_f_ens_R[2:11,j])
     time_10min_CTRL_f_ens_meanP_R[j]=time_10min_CTRL_f_ens_mean_R[j]+time_10min_CTRL_f_ens_stdv_R[j]
     time_10min_CTRL_f_ens_meanM_R[j]=time_10min_CTRL_f_ens_mean_R[j]-time_10min_CTRL_f_ens_stdv_R[j]

     precip_CTRL_f_ens_R[0,j]=0
     precip_CTRL_f_ens_R[1,j]=precip_CTRL_f_d1_R[j]
     precip_CTRL_f_ens_R[2,j]=precip_CTRL_f_d2_R[j]
     precip_CTRL_f_ens_R[3,j]=precip_CTRL_f_d3_R[j]
     precip_CTRL_f_ens_R[4,j]=precip_CTRL_f_d4_R[j]
     precip_CTRL_f_ens_R[5,j]=precip_CTRL_f_d5_R[j]
     precip_CTRL_f_ens_R[6,j]=precip_CTRL_f_d6_R[j]
     precip_CTRL_f_ens_R[7,j]=precip_CTRL_f_d7_R[j]
     precip_CTRL_f_ens_R[8,j]=precip_CTRL_f_d8_R[j]
     precip_CTRL_f_ens_R[9,j]=precip_CTRL_f_d9_R[j]
     precip_CTRL_f_ens_R[10,j]=precip_CTRL_f_d10_R[j]

     precip_CTRL_f_ens_R[:,47]=7.93208267e-03
     precip_CTRL_f_ens_R[:,48]=1.58285317e-02


     precip_CTRL_f_ens_mean_R[j]=precip_CTRL_f_ens_R[2:11,j].mean()
     precip_CTRL_f_ens_stdv_R[j]=np.std(precip_CTRL_f_ens_R[2:11,j])
     precip_CTRL_f_ens_meanP_R[j]=precip_CTRL_f_ens_mean_R[j]+precip_CTRL_f_ens_stdv_R[j]
     precip_CTRL_f_ens_meanM_R[j]=precip_CTRL_f_ens_mean_R[j]-precip_CTRL_f_ens_stdv_R[j]

     cld_MF_cb_CTRL_f_ens_R[0,j]=0
     cld_MF_cb_CTRL_f_ens_R[1,j]=cld_MF_cb_CTRL_f_d1_R[j]
     cld_MF_cb_CTRL_f_ens_R[2,j]=cld_MF_cb_CTRL_f_d2_R[j]
     cld_MF_cb_CTRL_f_ens_R[3,j]=cld_MF_cb_CTRL_f_d3_R[j]
     cld_MF_cb_CTRL_f_ens_R[4,j]=cld_MF_cb_CTRL_f_d4_R[j]
     cld_MF_cb_CTRL_f_ens_R[5,j]=cld_MF_cb_CTRL_f_d5_R[j]
     cld_MF_cb_CTRL_f_ens_R[6,j]=cld_MF_cb_CTRL_f_d6_R[j]
     cld_MF_cb_CTRL_f_ens_R[7,j]=cld_MF_cb_CTRL_f_d7_R[j]
     cld_MF_cb_CTRL_f_ens_R[8,j]=cld_MF_cb_CTRL_f_d8_R[j]
     cld_MF_cb_CTRL_f_ens_R[9,j]=cld_MF_cb_CTRL_f_d9_R[j]
     cld_MF_cb_CTRL_f_ens_R[10,j]=cld_MF_cb_CTRL_f_d10_R[j]

     cld_MF_cb_CTRL_f_ens_mean_R[j]=cld_MF_cb_CTRL_f_ens_R[2:11,j].mean()
     cld_MF_cb_CTRL_f_ens_stdv_R[j]=np.std(cld_MF_cb_CTRL_f_ens_R[2:11,j])
     cld_MF_cb_CTRL_f_ens_meanP_R[j]=cld_MF_cb_CTRL_f_ens_mean_R[j]+cld_MF_cb_CTRL_f_ens_stdv_R[j]
     cld_MF_cb_CTRL_f_ens_meanM_R[j]=cld_MF_cb_CTRL_f_ens_mean_R[j]-cld_MF_cb_CTRL_f_ens_stdv_R[j]

     cld_frac_cb_CTRL_f_ens_R[0,j]=0
     cld_frac_cb_CTRL_f_ens_R[1,j]=cld_frac_cb_CTRL_f_d1_R[j]
     cld_frac_cb_CTRL_f_ens_R[2,j]=cld_frac_cb_CTRL_f_d2_R[j]
     cld_frac_cb_CTRL_f_ens_R[3,j]=cld_frac_cb_CTRL_f_d3_R[j]
     cld_frac_cb_CTRL_f_ens_R[4,j]=cld_frac_cb_CTRL_f_d4_R[j]
     cld_frac_cb_CTRL_f_ens_R[5,j]=cld_frac_cb_CTRL_f_d5_R[j]
     cld_frac_cb_CTRL_f_ens_R[6,j]=cld_frac_cb_CTRL_f_d6_R[j]
     cld_frac_cb_CTRL_f_ens_R[7,j]=cld_frac_cb_CTRL_f_d7_R[j]
     cld_frac_cb_CTRL_f_ens_R[8,j]=cld_frac_cb_CTRL_f_d8_R[j]
     cld_frac_cb_CTRL_f_ens_R[9,j]=cld_frac_cb_CTRL_f_d9_R[j]
     cld_frac_cb_CTRL_f_ens_R[10,j]=cld_frac_cb_CTRL_f_d10_R[j]

     cld_frac_cb_CTRL_f_ens_mean_R[j]=cld_frac_cb_CTRL_f_ens_R[2:11,j].mean()
     cld_frac_cb_CTRL_f_ens_stdv_R[j]=np.std(cld_frac_cb_CTRL_f_ens_R[2:11,j])
     cld_frac_cb_CTRL_f_ens_meanP_R[j]=cld_frac_cb_CTRL_f_ens_mean_R[j]+cld_frac_cb_CTRL_f_ens_stdv_R[j]
     cld_frac_cb_CTRL_f_ens_meanM_R[j]=cld_frac_cb_CTRL_f_ens_mean_R[j]-cld_frac_cb_CTRL_f_ens_stdv_R[j]

     cldBCu_MF_cb_CTRL_f_ens_R[0,j]=0
     cldBCu_MF_cb_CTRL_f_ens_R[1,j]=cldBCu_MF_cb_CTRL_f_d1_R[j]
     cldBCu_MF_cb_CTRL_f_ens_R[2,j]=cldBCu_MF_cb_CTRL_f_d2_R[j]
     cldBCu_MF_cb_CTRL_f_ens_R[3,j]=cldBCu_MF_cb_CTRL_f_d3_R[j]
     cldBCu_MF_cb_CTRL_f_ens_R[4,j]=cldBCu_MF_cb_CTRL_f_d4_R[j]
     cldBCu_MF_cb_CTRL_f_ens_R[5,j]=cldBCu_MF_cb_CTRL_f_d5_R[j]
     cldBCu_MF_cb_CTRL_f_ens_R[6,j]=cldBCu_MF_cb_CTRL_f_d6_R[j]
     cldBCu_MF_cb_CTRL_f_ens_R[7,j]=cldBCu_MF_cb_CTRL_f_d7_R[j]
     cldBCu_MF_cb_CTRL_f_ens_R[8,j]=cldBCu_MF_cb_CTRL_f_d8_R[j]
     cldBCu_MF_cb_CTRL_f_ens_R[9,j]=cldBCu_MF_cb_CTRL_f_d9_R[j]
     cldBCu_MF_cb_CTRL_f_ens_R[10,j]=cldBCu_MF_cb_CTRL_f_d10_R[j]

     cldBCu_MF_cb_CTRL_f_ens_mean_R[j]=cldBCu_MF_cb_CTRL_f_ens_R[2:11,j].mean()
     cldBCu_MF_cb_CTRL_f_ens_stdv_R[j]=np.std(cldBCu_MF_cb_CTRL_f_ens_R[2:11,j])
     cldBCu_MF_cb_CTRL_f_ens_meanP_R[j]=cldBCu_MF_cb_CTRL_f_ens_mean_R[j]+cldBCu_MF_cb_CTRL_f_ens_stdv_R[j]
     cldBCu_MF_cb_CTRL_f_ens_meanM_R[j]=cldBCu_MF_cb_CTRL_f_ens_mean_R[j]-cldBCu_MF_cb_CTRL_f_ens_stdv_R[j]

     cldBCu_frac_cb_CTRL_f_ens_R[0,j]=0
     cldBCu_frac_cb_CTRL_f_ens_R[1,j]=cldBCu_frac_cb_CTRL_f_d1_R[j]
     cldBCu_frac_cb_CTRL_f_ens_R[2,j]=cldBCu_frac_cb_CTRL_f_d2_R[j]
     cldBCu_frac_cb_CTRL_f_ens_R[3,j]=cldBCu_frac_cb_CTRL_f_d3_R[j]
     cldBCu_frac_cb_CTRL_f_ens_R[4,j]=cldBCu_frac_cb_CTRL_f_d4_R[j]
     cldBCu_frac_cb_CTRL_f_ens_R[5,j]=cldBCu_frac_cb_CTRL_f_d5_R[j]
     cldBCu_frac_cb_CTRL_f_ens_R[6,j]=cldBCu_frac_cb_CTRL_f_d6_R[j]
     cldBCu_frac_cb_CTRL_f_ens_R[7,j]=cldBCu_frac_cb_CTRL_f_d7_R[j]
     cldBCu_frac_cb_CTRL_f_ens_R[8,j]=cldBCu_frac_cb_CTRL_f_d8_R[j]
     cldBCu_frac_cb_CTRL_f_ens_R[9,j]=cldBCu_frac_cb_CTRL_f_d9_R[j]
     cldBCu_frac_cb_CTRL_f_ens_R[10,j]=cldBCu_frac_cb_CTRL_f_d10_R[j]

     cldBCu_frac_cb_CTRL_f_ens_mean_R[j]=cldBCu_frac_cb_CTRL_f_ens_R[2:11,j].mean()
     cldBCu_frac_cb_CTRL_f_ens_stdv_R[j]=np.std(cldBCu_frac_cb_CTRL_f_ens_R[2:11,j])
     cldBCu_frac_cb_CTRL_f_ens_meanP_R[j]=cldBCu_frac_cb_CTRL_f_ens_mean_R[j]+cldBCu_frac_cb_CTRL_f_ens_stdv_R[j]
     cldBCu_frac_cb_CTRL_f_ens_meanM_R[j]=cldBCu_frac_cb_CTRL_f_ens_mean_R[j]-cldBCu_frac_cb_CTRL_f_ens_stdv_R[j]

     time_10min_PHALF_f_ens_R[0,j]=0
     time_10min_PHALF_f_ens_R[1,j]=time_10min_PHALF_f_d1_R[j]
     time_10min_PHALF_f_ens_R[2,j]=time_10min_PHALF_f_d2_R[j]
     time_10min_PHALF_f_ens_R[3,j]=time_10min_PHALF_f_d3_R[j]
     time_10min_PHALF_f_ens_R[4,j]=time_10min_PHALF_f_d4_R[j]
     time_10min_PHALF_f_ens_R[5,j]=time_10min_PHALF_f_d5_R[j]
     time_10min_PHALF_f_ens_R[6,j]=time_10min_PHALF_f_d6_R[j]
     time_10min_PHALF_f_ens_R[7,j]=time_10min_PHALF_f_d7_R[j]
     time_10min_PHALF_f_ens_R[8,j]=time_10min_PHALF_f_d8_R[j]
     time_10min_PHALF_f_ens_R[9,j]=time_10min_PHALF_f_d9_R[j]
     time_10min_PHALF_f_ens_R[10,j]=time_10min_PHALF_f_d10_R[j]

     time_10min_PHALF_f_ens_mean_R[j]=time_10min_PHALF_f_ens_R[2:11,j].mean()
     time_10min_PHALF_f_ens_stdv_R[j]=np.std(time_10min_PHALF_f_ens_R[2:11,j])
     time_10min_PHALF_f_ens_meanP_R[j]=time_10min_PHALF_f_ens_mean_R[j]+time_10min_PHALF_f_ens_stdv_R[j]
     time_10min_PHALF_f_ens_meanM_R[j]=time_10min_PHALF_f_ens_mean_R[j]-time_10min_PHALF_f_ens_stdv_R[j]

     precip_PHALF_f_ens_R[0,j]=0
     precip_PHALF_f_ens_R[1,j]=precip_PHALF_f_d1_R[j]
     precip_PHALF_f_ens_R[2,j]=precip_PHALF_f_d2_R[j]
     precip_PHALF_f_ens_R[3,j]=precip_PHALF_f_d3_R[j]
     precip_PHALF_f_ens_R[4,j]=precip_PHALF_f_d4_R[j]
     precip_PHALF_f_ens_R[5,j]=precip_PHALF_f_d5_R[j]
     precip_PHALF_f_ens_R[6,j]=precip_PHALF_f_d6_R[j]
     precip_PHALF_f_ens_R[7,j]=precip_PHALF_f_d7_R[j]
     precip_PHALF_f_ens_R[8,j]=precip_PHALF_f_d8_R[j]
     precip_PHALF_f_ens_R[9,j]=precip_PHALF_f_d9_R[j]
     precip_PHALF_f_ens_R[10,j]=precip_PHALF_f_d10_R[j]

     precip_PHALF_f_ens_mean_R[j]=precip_PHALF_f_ens_R[2:11,j].mean()
     precip_PHALF_f_ens_stdv_R[j]=np.std(precip_PHALF_f_ens_R[2:11,j])
     precip_PHALF_f_ens_meanP_R[j]=precip_PHALF_f_ens_mean_R[j]+precip_PHALF_f_ens_stdv_R[j]
     precip_PHALF_f_ens_meanM_R[j]=precip_PHALF_f_ens_mean_R[j]-precip_PHALF_f_ens_stdv_R[j]

     cld_MF_cb_PHALF_f_ens_R[0,j]=0
     cld_MF_cb_PHALF_f_ens_R[1,j]=cld_MF_cb_PHALF_f_d1_R[j]
     cld_MF_cb_PHALF_f_ens_R[2,j]=cld_MF_cb_PHALF_f_d2_R[j]
     cld_MF_cb_PHALF_f_ens_R[3,j]=cld_MF_cb_PHALF_f_d3_R[j]
     cld_MF_cb_PHALF_f_ens_R[4,j]=cld_MF_cb_PHALF_f_d4_R[j]
     cld_MF_cb_PHALF_f_ens_R[5,j]=cld_MF_cb_PHALF_f_d5_R[j]
     cld_MF_cb_PHALF_f_ens_R[6,j]=cld_MF_cb_PHALF_f_d6_R[j]
     cld_MF_cb_PHALF_f_ens_R[7,j]=cld_MF_cb_PHALF_f_d7_R[j]
     cld_MF_cb_PHALF_f_ens_R[8,j]=cld_MF_cb_PHALF_f_d8_R[j]
     cld_MF_cb_PHALF_f_ens_R[9,j]=cld_MF_cb_PHALF_f_d9_R[j]
     cld_MF_cb_PHALF_f_ens_R[10,j]=cld_MF_cb_PHALF_f_d10_R[j]

     cld_MF_cb_PHALF_f_ens_mean_R[j]=cld_MF_cb_PHALF_f_ens_R[2:11,j].mean()
     cld_MF_cb_PHALF_f_ens_stdv_R[j]=np.std(cld_MF_cb_PHALF_f_ens_R[2:11,j])
     cld_MF_cb_PHALF_f_ens_meanP_R[j]=cld_MF_cb_PHALF_f_ens_mean_R[j]+cld_MF_cb_PHALF_f_ens_stdv_R[j]
     cld_MF_cb_PHALF_f_ens_meanM_R[j]=cld_MF_cb_PHALF_f_ens_mean_R[j]-cld_MF_cb_PHALF_f_ens_stdv_R[j]

     time_10min_MHALF_f_ens_R[0,j]=0
     time_10min_MHALF_f_ens_R[1,j]=time_10min_MHALF_f_d1_R[j]
     time_10min_MHALF_f_ens_R[2,j]=time_10min_MHALF_f_d2_R[j]
     time_10min_MHALF_f_ens_R[3,j]=time_10min_MHALF_f_d3_R[j]
     time_10min_MHALF_f_ens_R[4,j]=time_10min_MHALF_f_d4_R[j]
     time_10min_MHALF_f_ens_R[5,j]=time_10min_MHALF_f_d5_R[j]
     time_10min_MHALF_f_ens_R[6,j]=time_10min_MHALF_f_d6_R[j]
     time_10min_MHALF_f_ens_R[7,j]=time_10min_MHALF_f_d7_R[j]
     time_10min_MHALF_f_ens_R[8,j]=time_10min_MHALF_f_d8_R[j]
     time_10min_MHALF_f_ens_R[9,j]=time_10min_MHALF_f_d9_R[j]
     time_10min_MHALF_f_ens_R[10,j]=time_10min_MHALF_f_d10_R[j]
  
     time_10min_MHALF_f_ens_mean_R[j]=time_10min_MHALF_f_ens_R[2:11,j].mean()
     time_10min_MHALF_f_ens_stdv_R[j]=np.std(time_10min_MHALF_f_ens_R[2:11,j])
     time_10min_MHALF_f_ens_meanP_R[j]=time_10min_MHALF_f_ens_mean_R[j]+time_10min_MHALF_f_ens_stdv_R[j]
     time_10min_MHALF_f_ens_meanM_R[j]=time_10min_MHALF_f_ens_mean_R[j]-time_10min_MHALF_f_ens_stdv_R[j]

     precip_MHALF_f_ens_R[0,j]=0
     precip_MHALF_f_ens_R[1,j]=precip_MHALF_f_d1_R[j]
     precip_MHALF_f_ens_R[2,j]=precip_MHALF_f_d2_R[j]
     precip_MHALF_f_ens_R[3,j]=precip_MHALF_f_d3_R[j]
     precip_MHALF_f_ens_R[4,j]=precip_MHALF_f_d4_R[j]
     precip_MHALF_f_ens_R[5,j]=precip_MHALF_f_d5_R[j]
     precip_MHALF_f_ens_R[6,j]=precip_MHALF_f_d6_R[j]
     precip_MHALF_f_ens_R[7,j]=precip_MHALF_f_d7_R[j]
     precip_MHALF_f_ens_R[8,j]=precip_MHALF_f_d8_R[j]
     precip_MHALF_f_ens_R[9,j]=precip_MHALF_f_d9_R[j]
     precip_MHALF_f_ens_R[10,j]=precip_MHALF_f_d10_R[j]

     precip_MHALF_f_ens_mean_R[j]=precip_MHALF_f_ens_R[2:11,j].mean()
     precip_MHALF_f_ens_stdv_R[j]=np.std(precip_MHALF_f_ens_R[2:11,j])
     precip_MHALF_f_ens_meanP_R[j]=precip_MHALF_f_ens_mean_R[j]+precip_MHALF_f_ens_stdv_R[j]
     precip_MHALF_f_ens_meanM_R[j]=precip_MHALF_f_ens_mean_R[j]-precip_MHALF_f_ens_stdv_R[j]

     cld_MF_cb_MHALF_f_ens_R[0,j]=0
     cld_MF_cb_MHALF_f_ens_R[1,j]=cld_MF_cb_MHALF_f_d1_R[j]
     cld_MF_cb_MHALF_f_ens_R[2,j]=cld_MF_cb_MHALF_f_d2_R[j]
     cld_MF_cb_MHALF_f_ens_R[3,j]=cld_MF_cb_MHALF_f_d3_R[j]
     cld_MF_cb_MHALF_f_ens_R[4,j]=cld_MF_cb_MHALF_f_d4_R[j]
     cld_MF_cb_MHALF_f_ens_R[5,j]=cld_MF_cb_MHALF_f_d5_R[j]
     cld_MF_cb_MHALF_f_ens_R[6,j]=cld_MF_cb_MHALF_f_d6_R[j]
     cld_MF_cb_MHALF_f_ens_R[7,j]=cld_MF_cb_MHALF_f_d7_R[j]
     cld_MF_cb_MHALF_f_ens_R[8,j]=cld_MF_cb_MHALF_f_d8_R[j]
     cld_MF_cb_MHALF_f_ens_R[9,j]=cld_MF_cb_MHALF_f_d9_R[j]
     cld_MF_cb_MHALF_f_ens_R[10,j]=cld_MF_cb_MHALF_f_d10_R[j]

     cld_MF_cb_MHALF_f_ens_mean_R[j]=cld_MF_cb_MHALF_f_ens_R[2:11,j].mean()
     cld_MF_cb_MHALF_f_ens_stdv_R[j]=np.std(cld_MF_cb_MHALF_f_ens_R[2:11,j])
     cld_MF_cb_MHALF_f_ens_meanP_R[j]=cld_MF_cb_MHALF_f_ens_mean_R[j]+cld_MF_cb_MHALF_f_ens_stdv_R[j]
     cld_MF_cb_MHALF_f_ens_meanM_R[j]=cld_MF_cb_MHALF_f_ens_mean_R[j]-cld_MF_cb_MHALF_f_ens_stdv_R[j]

     precip_PHALF_f_ens_RC[:,j]=precip_PHALF_f_ens_R[:,j]/0.3
     precip_MHALF_f_ens_RC[:,j]=precip_MHALF_f_ens_R[:,j]/0.1
     precip_CTRL_f_ens_RC[:,j]=precip_CTRL_f_ens_R[:,j]/0.2

     precip_PHALF_f_ens_mean_RC[j]=precip_PHALF_f_ens_RC[2:11,j].mean()
     precip_PHALF_f_ens_stdv_RC[j]=np.std(precip_PHALF_f_ens_RC[2:11,j])
     precip_PHALF_f_ens_meanP_RC[j]=precip_PHALF_f_ens_mean_RC[j]+precip_PHALF_f_ens_stdv_RC[j]
     precip_PHALF_f_ens_meanM_RC[j]=precip_PHALF_f_ens_mean_RC[j]-precip_PHALF_f_ens_stdv_RC[j]
     precip_PHALF_f_ens_meanM_RC[48]=precip_PHALF_f_ens_mean_RC[48]-0.75*precip_PHALF_f_ens_stdv_RC[48]

     precip_MHALF_f_ens_mean_RC[j]=precip_MHALF_f_ens_RC[2:11,j].mean()
     precip_MHALF_f_ens_stdv_RC[j]=np.std(precip_MHALF_f_ens_RC[2:11,j])
     precip_MHALF_f_ens_meanP_RC[j]=precip_MHALF_f_ens_mean_RC[j]+precip_MHALF_f_ens_stdv_RC[j]
     precip_MHALF_f_ens_meanM_RC[j]=precip_MHALF_f_ens_mean_RC[j]-precip_MHALF_f_ens_stdv_RC[j]

     precip_CTRL_f_ens_mean_RC[j]=precip_CTRL_f_ens_RC[2:11,j].mean()
     precip_CTRL_f_ens_stdv_RC[j]=np.std(precip_CTRL_f_ens_RC[2:11,j])
     precip_CTRL_f_ens_meanP_RC[j]=precip_CTRL_f_ens_mean_RC[j]+precip_CTRL_f_ens_stdv_RC[j]
     precip_CTRL_f_ens_meanM_RC[j]=precip_CTRL_f_ens_mean_RC[j]-precip_CTRL_f_ens_stdv_RC[j]


     time_10min_RTHQV_ens_R[0,j]=0.0
     time_10min_RTHQV_ens_R[1,j]=0.0
     time_10min_RTHQV_ens_R[2,j]=time_10min_RTHQVD2_d2_R[j]
     time_10min_RTHQV_ens_R[3,j]=time_10min_RTHQVD3_d3_R[j]
     time_10min_RTHQV_ens_R[4,j]=time_10min_RTHQVD4_d4_R[j]
     time_10min_RTHQV_ens_R[5,j]=time_10min_RTHQVD5_d5_R[j]
     time_10min_RTHQV_ens_R[6,j]=time_10min_RTHQVD6_d6_R[j]
     time_10min_RTHQV_ens_R[7,j]=time_10min_RTHQVD7_d7_R[j]
     time_10min_RTHQV_ens_R[8,j]=time_10min_RTHQVD8_d8_R[j]
     time_10min_RTHQV_ens_R[9,j]=time_10min_RTHQVD9_d9_R[j]
     time_10min_RTHQV_ens_R[10,j]=time_10min_RTHQVD10_d10_R[j]

     time_10min_RTHQV_ens_mean_R[j]=time_10min_RTHQV_ens_R[2:11,j].mean()
     time_10min_RTHQV_ens_stdv_R[j]=np.std(time_10min_RTHQV_ens_R[2:11,j])
     time_10min_RTHQV_ens_meanP_R[j]=time_10min_RTHQV_ens_mean_R[j]+time_10min_RTHQV_ens_stdv_R[j]
     time_10min_RTHQV_ens_meanM_R[j]=time_10min_RTHQV_ens_mean_R[j]-time_10min_RTHQV_ens_stdv_R[j]

     precip_RTHQV_ens_R[0,j]=0.0
     precip_RTHQV_ens_R[1,j]=0.0
     precip_RTHQV_ens_R[2,j]=precip_RTHQVD2_d2_R[j]
     precip_RTHQV_ens_R[3,j]=precip_RTHQVD3_d3_R[j]
     precip_RTHQV_ens_R[4,j]=precip_RTHQVD4_d4_R[j]
     precip_RTHQV_ens_R[5,j]=precip_RTHQVD5_d5_R[j]
     precip_RTHQV_ens_R[6,j]=precip_RTHQVD6_d6_R[j]
     precip_RTHQV_ens_R[7,j]=precip_RTHQVD7_d7_R[j]
     precip_RTHQV_ens_R[8,j]=precip_RTHQVD8_d8_R[j]
     precip_RTHQV_ens_R[9,j]=precip_RTHQVD9_d9_R[j]
     precip_RTHQV_ens_R[10,j]=precip_RTHQVD10_d10_R[j]

     precip_RTHQV_ens_mean_R[j]=precip_RTHQV_ens_R[2:11,j].mean()
     precip_RTHQV_ens_stdv_R[j]=np.std(precip_RTHQV_ens_R[2:11,j])
     precip_RTHQV_ens_meanP_R[j]=precip_RTHQV_ens_mean_R[j]+precip_RTHQV_ens_stdv_R[j]
     precip_RTHQV_ens_meanM_R[j]=precip_RTHQV_ens_mean_R[j]-precip_RTHQV_ens_stdv_R[j]

     precip_RTHQV_ens_RC[:,j]=precip_RTHQV_ens_R[:,j]/0.18
     precip_RTHQV_ens_mean_RC[j]=precip_RTHQV_ens_RC[2:11,j].mean()
     precip_RTHQV_ens_stdv_RC[j]=np.std(precip_RTHQV_ens_RC[2:11,j])
     precip_RTHQV_ens_meanP_RC[j]=precip_RTHQV_ens_mean_RC[j]+precip_RTHQV_ens_stdv_RC[j]
     precip_RTHQV_ens_meanM_RC[j]=precip_RTHQV_ens_mean_RC[j]-precip_RTHQV_ens_stdv_RC[j]

     cld_MF_cb_RTHQV_ens_R[0,j]=0.0
     cld_MF_cb_RTHQV_ens_R[1,j]=0.0
     cld_MF_cb_RTHQV_ens_R[2,j]=cld_MF_cb_RTHQVD2_d2_R[j]
     cld_MF_cb_RTHQV_ens_R[3,j]=cld_MF_cb_RTHQVD3_d3_R[j]
     cld_MF_cb_RTHQV_ens_R[4,j]=cld_MF_cb_RTHQVD4_d4_R[j]
     cld_MF_cb_RTHQV_ens_R[5,j]=cld_MF_cb_RTHQVD5_d5_R[j]
     cld_MF_cb_RTHQV_ens_R[6,j]=cld_MF_cb_RTHQVD6_d6_R[j]
     cld_MF_cb_RTHQV_ens_R[7,j]=cld_MF_cb_RTHQVD7_d7_R[j]
     cld_MF_cb_RTHQV_ens_R[8,j]=cld_MF_cb_RTHQVD8_d8_R[j]
     cld_MF_cb_RTHQV_ens_R[9,j]=cld_MF_cb_RTHQVD9_d9_R[j]
     cld_MF_cb_RTHQV_ens_R[10,j]=cld_MF_cb_RTHQVD10_d10_R[j]

     cld_MF_cb_RTHQV_ens_mean_R[j]=cld_MF_cb_RTHQV_ens_R[2:11,j].mean()
     cld_MF_cb_RTHQV_ens_stdv_R[j]=np.std(cld_MF_cb_RTHQV_ens_R[2:11,j])
     cld_MF_cb_RTHQV_ens_meanP_R[j]=cld_MF_cb_RTHQV_ens_mean_R[j]+cld_MF_cb_RTHQV_ens_stdv_R[j]
     cld_MF_cb_RTHQV_ens_meanM_R[j]=cld_MF_cb_RTHQV_ens_mean_R[j]-cld_MF_cb_RTHQV_ens_stdv_R[j]

#np.save('time_10min_RTHQV_ens_mean',time_10min_RTHQV_ens_mean_R)
#np.save('time_10min_CTRL_f_ens_mean',time_10min_CTRL_f_ens_mean_R)
#np.save('time_10min_S65L200_ens_mean',time_10min_MHALF_f_ens_mean_R)
#np.save('time_10min_S195L600_ens_mean',time_10min_PHALF_f_ens_mean_R)

n_30min=48


#
#rthqv
#
#nr=18+n_30min 

nr=36+n_10min

time_10min_PHALF_f_ens_mean_R_shift0=N.zeros((100))
time_10min_MHALF_f_ens_mean_R_shift0=N.zeros((100))
time_10min_CTRL_f_ens_mean_R_shift0=N.zeros((100))
time_10min_PHALF_f_shift0=N.zeros((100))
time_10min_MHALF_f_shift0=N.zeros((100))
time_10min_CTRL_f_shift0=N.zeros((100))

time_10min_RTHQV_ens_mean_R_shift0=N.zeros((100))
time_10min_RTHQV_shift0=N.zeros((100))

time_10min_CTRL_f=time_10min_CTRL_f/24.0
time_10min_PHALF_f=time_10min_PHALF_f/24.0
time_10min_MHALF_f=time_10min_MHALF_f/24.0

time_10min_PHALF_f_ens_mean_R_shift0[45:100]=time_10min_PHALF_f_ens_mean_R[45:100]- 2.5
time_10min_MHALF_f_ens_mean_R_shift0[48:100]=time_10min_MHALF_f_ens_mean_R[48:100] -3.25
time_10min_CTRL_f_ens_mean_R_shift0[46:100]=time_10min_CTRL_f_ens_mean_R[46:100]- 2.75
time_10min_RTHQV_ens_mean_R_shift0[46:100]=time_10min_RTHQV_ens_mean_R[46:100]-2.75


time_10min_CTRL_f[90:182]=24.*time_10min_CTRL_f[90:182]-24.0
print time_10min_CTRL_f[90:182]
print hein
print np.max(ACuzt_MF_time_CTRL_f[1:85,90:182])

fig = plt.figure(figsize=(16,13))
cld_levs = (N.arange(-0.8, 1, 0.2))
w_levs   = (N.arange(-0.8, 1, 0.2))

plt.subplot(211)
x=time_10min_CTRL_f[90:182]
y=hein[1:85]
X, Y = np.meshgrid(x, y)
CS = plt.contourf(X, Y ,ACuzt_MF_time_CTRL_f[1:85,90:182], levels=[0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.08,0.10,0.14])
cbar = plt.colorbar()

CSl = plt.contour(X, Y ,BCuzt_MF_time_CTRL_f[1:85,90:182], colors='k', levels=[0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14])
plt.clabel(CSl, fontsize=9, inline=1, fmt='%3.3f')
plt.ylabel('Height (km)', fontsize='x-large')
plt.title('Mass Flux (kg/m$^2$/s)', fontsize='x-large')
plt.xlabel('Time (Hour) ', fontsize='x-large')
axes = plt.gca()
axes.set_xlim([0,24])
axes.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
axes.set_xticklabels([0,2,4,6,8,10,12,14,16,18,20,22,24])
axes = plt.gca()
axes.set_ylim([0, 14])
axes.set_yticks([0,2,4,6,8,10,12, 14])
axes.set_yticklabels([0,2,4,6,8,10,12, 14])
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(1, 13, r'a)', fontdict=font, color="black", fontsize='x-large')

plt.subplot(212)
x=time_10min_CTRL_f[90:182]
y=hein[1:85]
X, Y = np.meshgrid(x, y)
CS = plt.contourf(X, Y ,100.*ACuzt_frac_time_CTRL_f[1:85,90:182], levels=[0.5,1,3,5,7,9,11,16,21,26])
cbar = plt.colorbar()
CSl = plt.contour(X, Y ,100.*BCuzt_frac_time_CTRL_f[1:85,90:182], colors='k', levels=[0.5,1,3,5,7,9,11,13])
plt.clabel(CSl, fontsize=9, inline=1, fmt='%3.1f')
plt.title('$\sigma$ (%)', fontsize='x-large')
plt.xlabel('Time (Hour) ', fontsize='x-large')
plt.ylabel('Height (km)', fontsize='x-large')
axes = plt.gca()
axes.set_xlim([0,24])
axes.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
axes.set_xticklabels([0,2,4,6,8,10,12,14,16,18,20,22,24])

axes = plt.gca()
axes.set_ylim([0, 14])
axes.set_yticks([0,2,4,6,8,10,12, 14])
axes.set_yticklabels([0,2,4,6,8,10,12, 14])
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(1, 13, r'b)', fontdict=font, color="black", fontsize='x-large')
plt.savefig('Figure1.png')


fig = plt.figure(figsize=(16,20))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
plt.subplot(311)
plt.plot(time_10min_CTRL_f_ens_mean_R[0:100], precip_CTRL_f_ens_mean_R[0:100], lw=4, color="black", label='Control')
plt.fill_between(time_10min_CTRL_f_ens_mean_R[0:100], precip_CTRL_f_ens_meanM_R[0:100], precip_CTRL_f_ens_meanP_R[0:100], alpha=0.5, edgecolor='grey', facecolor='grey')

plt.plot(time_SF[0:100], (surf_flx_CTRL_f[0:100])/1000.0, lw=1,linestyle=':', color="black", label='(shf+lhf)/e3')


plt.ylabel('Precipitation (mm/h)', fontsize='x-large')
plt.xlabel('Time (hour) ', fontsize='x-large')
plt.xlim(0,15)
axes = plt.gca()
axes.set_xlim([0, 15])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])

plt.ylim(0,1.2)
axes = plt.gca()
axes.set_ylim([0, 1.2])
axes.set_yticks([0.,0.3,0.6,0.9,1.2])
axes.set_yticklabels([0.,0.3,0.6,0.9,1.2])
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 1.1, r'a)', fontdict=font, color="black", fontsize='x-large')


plt.subplot(312)
plt.plot(time_SF[0:100], (surf_flx_CTRL_f[0:100])/10000.0, lw=1,linestyle=':', color="black", label='(shf+lhf)/e4')
plt.plot(time_10min_CTRL_f_ens_mean_R[0:100], cld_MF_cb_CTRL_f_ens_mean_R[0:100], lw=4, color="black", label='Ens mean')
plt.fill_between(time_10min_CTRL_f_ens_mean_R[0:100], cld_MF_cb_CTRL_f_ens_meanM_R[0:100], cld_MF_cb_CTRL_f_ens_meanP_R[0:100], alpha=0.5, edgecolor='grey', facecolor='grey')

plt.plot(time_10min_CTRL_f_ens_mean_R[0:100], cldBCu_MF_cb_CTRL_f_ens_mean_R[0:100], lw=4, color="red", label='Ens mean')
plt.fill_between(time_10min_CTRL_f_ens_mean_R[0:100], cldBCu_MF_cb_CTRL_f_ens_meanM_R[0:100], cldBCu_MF_cb_CTRL_f_ens_meanP_R[0:100], alpha=0.5, edgecolor='lightcoral', facecolor='lightcoral')

plt.ylabel('Mass Flux (kg/m$^2$/s)', fontsize='x-large')
plt.xlabel('Time (hour) ', fontsize='x-large')
plt.xlim(0,15)
axes = plt.gca()
axes.set_xlim([0, 15])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])

plt.ylim(0,0.15)
axes = plt.gca()
axes.set_ylim([0, 0.15])
axes.set_yticks([0,0.05,0.1,0.15])
axes.set_yticklabels([0,0.05,0.1,0.15])
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.13, r'b)', fontdict=font, color="black", fontsize='x-large')

plt.subplot(313)
plt.plot(time_SF[0:100], (surf_flx_CTRL_f[0:100])/100.0,linestyle=':', lw=1, color="black", label='(shf+lhf)/e4')
plt.plot(time_10min_CTRL_f_ens_mean_R[0:100], 100.*cld_frac_cb_CTRL_f_ens_mean_R[0:100], lw=4, color="black", label='Ens mean')
plt.fill_between(time_10min_CTRL_f_ens_mean_R[0:100], 100.*cld_frac_cb_CTRL_f_ens_meanM_R[0:100], 100.*cld_frac_cb_CTRL_f_ens_meanP_R[0:100], alpha=0.5, edgecolor='grey', facecolor='grey')

plt.plot(time_10min_CTRL_f_ens_mean_R[0:100], 100.*cldBCu_frac_cb_CTRL_f_ens_mean_R[0:100], lw=4, color="red", label='Ens mean')
plt.fill_between(time_10min_CTRL_f_ens_mean_R[0:100], 100.*cldBCu_frac_cb_CTRL_f_ens_meanM_R[0:100], 100.*cldBCu_frac_cb_CTRL_f_ens_meanP_R[0:100], alpha=0.5, edgecolor='lightcoral', facecolor='lightcoral')

plt.ylabel('$\sigma$ (%)', fontsize='x-large')
plt.xlabel('Time (hour) ', fontsize='x-large')
plt.xlim(0,15)
axes = plt.gca()
axes.set_xlim([0, 15])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])


plt.ylim(0,10)
axes = plt.gca()
axes.set_ylim([0, 10])
axes.set_yticks([0,2,4,6,8,10])
axes.set_yticklabels([0,2,4,6,8,10])
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 9, r'c)', fontdict=font, color="black", fontsize='x-large')

plt.savefig('Figure2.png')



fig = plt.figure(figsize=(20,13))
ax = fig.add_axes([0.1,0.1,0.8,0.8])

plt.subplot(221)
plt.plot(time_10min_CTRL_f_ens_mean_R[0:100], precip_CTRL_f_ens_mean_RC[0:100], lw=4, color="black", label='C')
plt.fill_between(time_10min_CTRL_f_ens_mean_R[0:100], precip_CTRL_f_ens_meanM_RC[0:100], precip_CTRL_f_ens_meanP_RC[0:100], alpha=0.5, edgecolor='grey', facecolor='grey')

plt.plot(time_10min_PHALF_f_ens_mean_R[0:100], precip_PHALF_f_ens_mean_RC[0:100], lw=4, color="red", label='S')
plt.fill_between(time_10min_PHALF_f_ens_mean_R[0:100], precip_PHALF_f_ens_meanM_RC[0:100], precip_PHALF_f_ens_meanP_RC[0:100], alpha=0.5, edgecolor='lightcoral', facecolor='lightcoral')

plt.plot(time_10min_MHALF_f_ens_mean_R[0:100], precip_MHALF_f_ens_mean_RC[0:100], lw=4, color="blue", label='W')
plt.fill_between(time_10min_MHALF_f_ens_mean_R[0:100], precip_MHALF_f_ens_meanM_RC[0:100], precip_MHALF_f_ens_meanP_RC[0:100], alpha=0.5, edgecolor='cornflowerblue', facecolor='cornflowerblue')

plt.plot(time_SF[0:100], (surf_flx_CTRL_f[0:100])/200.0, lw=1,linestyle=':', color="black", label='(shf+lhf)/e3')
plt.plot(time_SF[0:100], (surf_flx_PHALF_f[0:100])/200.0, lw=1,linestyle=':', color="red", label='(shf+lhf)/e3')
plt.plot(time_SF[0:100], (surf_flx_MHALF_f[0:100])/200.0, lw=1,linestyle=':', color="blue", label='(shf+lhf)/e3')

plt.ylabel('Normalized Precipitation', fontsize='x-large')
plt.xlabel('Time (hour) ', fontsize='x-large')
plt.xlim(0,15)
axes = plt.gca()
axes.set_xlim([0, 15])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])

plt.ylim(0,6)
axes = plt.gca()
axes.set_ylim([0, 6])
axes.set_yticks([0,1,2,3,4,5,6])
axes.set_yticklabels([0,1,2,3,4,5,6])
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 5.6, r'a)', fontdict=font, color="black", fontsize='x-large')

plt.subplot(223)
plt.plot(time_SF[0:100], (surf_flx_CTRL_f[0:100])/10000.0, lw=1,linestyle=':', color="black", label='(shf+lhf)/e4')
plt.plot(time_10min_CTRL_f_ens_mean_R[0:100], cld_MF_cb_CTRL_f_ens_mean_R[0:100], lw=4, color="black", label='Ens mean')
plt.fill_between(time_10min_CTRL_f_ens_mean_R[0:100], cld_MF_cb_CTRL_f_ens_meanM_R[0:100], cld_MF_cb_CTRL_f_ens_meanP_R[0:100], alpha=0.5, edgecolor='grey', facecolor='grey')

plt.plot(time_SF[0:100], (surf_flx_PHALF_f[0:100])/10000.0, lw=1,linestyle=':', color="red", label='(shf+lhf)/e4')
plt.plot(time_10min_PHALF_f_ens_mean_R[0:100], cld_MF_cb_PHALF_f_ens_mean_R[0:100], lw=4, color="red", label='Ens mean')
plt.fill_between(time_10min_PHALF_f_ens_mean_R[0:100], cld_MF_cb_PHALF_f_ens_meanM_R[0:100], cld_MF_cb_PHALF_f_ens_meanP_R[0:100], alpha=0.5, edgecolor='lightcoral', facecolor='lightcoral')

plt.plot(time_SF[0:100], (surf_flx_MHALF_f[0:100])/10000.0, lw=1,linestyle=':', color="blue", label='(shf+lhf)/e4')
plt.plot(time_10min_MHALF_f_ens_mean_R[0:100], cld_MF_cb_MHALF_f_ens_mean_R[0:100], lw=4, color="blue", label='Ens mean')
plt.fill_between(time_10min_MHALF_f_ens_mean_R[0:100], cld_MF_cb_MHALF_f_ens_meanM_R[0:100], cld_MF_cb_MHALF_f_ens_meanP_R[0:100], alpha=0.5, edgecolor='cornflowerblue', facecolor='cornflowerblue')

plt.ylabel('Mass Flux (kg/m$^2$/s)', fontsize='x-large')
plt.xlabel('Time (hour) ', fontsize='x-large')
plt.xlim(0,15)
axes = plt.gca()
axes.set_xlim([0, 15])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])

plt.ylim(0,0.16)
axes = plt.gca()
axes.set_ylim([0, 0.16])
axes.set_yticks([0,0.04,0.08,0.12,0.16])
axes.set_yticklabels([0,0.04,0.08,0.12,0.16])
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.15, r'b)', fontdict=font, color="black", fontsize='x-large')


plt.subplot(222)
plt.plot(time_10min_CTRL_f_ens_mean_R_shift0[46:100], precip_CTRL_f_ens_mean_RC[46:100], lw=4, color="black", label='C')
plt.fill_between(time_10min_CTRL_f_ens_mean_R_shift0[46:100], precip_CTRL_f_ens_meanM_RC[46:100], precip_CTRL_f_ens_meanP_RC[46:100], alpha=0.5, edgecolor='grey', facecolor='grey')


plt.plot(time_10min_PHALF_f_ens_mean_R_shift0[45:100], precip_PHALF_f_ens_mean_RC[45:100], lw=4, color="red", label='S')
plt.fill_between(time_10min_PHALF_f_ens_mean_R_shift0[45:100], precip_PHALF_f_ens_meanM_RC[45:100], precip_PHALF_f_ens_meanP_RC[45:100], alpha=0.5, edgecolor='lightcoral', facecolor='lightcoral')


plt.plot(time_10min_MHALF_f_ens_mean_R_shift0[48:100], precip_MHALF_f_ens_mean_RC[48:100], lw=4, color="blue", label='W')
plt.fill_between(time_10min_MHALF_f_ens_mean_R_shift0[48:100], precip_MHALF_f_ens_meanM_RC[48:100], precip_MHALF_f_ens_meanP_RC[48:100], alpha=0.5, edgecolor='cornflowerblue', facecolor='cornflowerblue')
plt.legend()


plt.ylabel('Normalized Precipitation', fontsize='x-large')
plt.xlabel('Time after triggering (hour) ', fontsize='x-large')
plt.xlim(0,12)
axes = plt.gca()
axes.set_xlim([0,12])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10,11,12])

plt.ylim(0,6)
axes = plt.gca()
axes.set_ylim([0, 6])
axes.set_yticks([0,1,2,3,4,5,6])
axes.set_yticklabels([0,1,2,3,4,5,6])
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 5.6, r'c)', fontdict=font, color="black", fontsize='x-large')

plt.subplot(224)

plt.plot(time_10min_CTRL_f_ens_mean_R_shift0[46:100], cld_MF_cb_CTRL_f_ens_mean_R[46:100], lw=4, color="black", label='Ens mean')
plt.fill_between(time_10min_CTRL_f_ens_mean_R_shift0[46:100], cld_MF_cb_CTRL_f_ens_meanM_R[46:100], cld_MF_cb_CTRL_f_ens_meanP_R[46:100], alpha=0.5, edgecolor='grey', facecolor='grey')


plt.plot(time_10min_PHALF_f_ens_mean_R_shift0[45:100], cld_MF_cb_PHALF_f_ens_mean_R[45:100], lw=4, color="red", label='Ens mean')
plt.fill_between(time_10min_PHALF_f_ens_mean_R_shift0[45:100], cld_MF_cb_PHALF_f_ens_meanM_R[45:100], cld_MF_cb_PHALF_f_ens_meanP_R[45:100], alpha=0.5, edgecolor='lightcoral', facecolor='lightcoral')


plt.plot(time_10min_MHALF_f_ens_mean_R_shift0[48:100], cld_MF_cb_MHALF_f_ens_mean_R[48:100], lw=4, color="blue", label='Ens mean')
plt.fill_between(time_10min_MHALF_f_ens_mean_R_shift0[48:100], cld_MF_cb_MHALF_f_ens_meanM_R[48:100], cld_MF_cb_MHALF_f_ens_meanP_R[48:100], alpha=0.5, edgecolor='cornflowerblue', facecolor='cornflowerblue')

#plt.legend()



plt.ylabel('Mass Flux (kg/m$^2$/s)', fontsize='x-large')
plt.xlabel('Time after triggering (hour) ', fontsize='x-large')
plt.xlim(0,12)
axes = plt.gca()
axes.set_xlim([0,12])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10,11,12])




plt.ylim(0,0.16)
axes = plt.gca()
axes.set_ylim([0, 0.16])
axes.set_yticks([0,0.04,0.08,0.12,0.16])
axes.set_yticklabels([0,0.04,0.08,0.12,0.16])
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.15, r'd)', fontdict=font, color="black", fontsize='x-large')

plt.savefig('Figure3.png')


fig = plt.figure(figsize=(20,16))
cld_levs = (N.arange(-0.8, 1, 0.2))
w_levs   = (N.arange(-0.8, 1, 0.2))

plt.subplot(321)
y=hein[1:85]
x=length_y[:]/1000
X, Y = np.meshgrid(x, y)
CS = plt.contourf(X, Y ,thzy24h_snap_CTRL_f[1:85,:], levels=[-0.3,-0.2,-0.1,0,0.1,0.2,0.3], cmap="RdBu_r")
cbar = plt.colorbar()
plt.xlabel('X (km)')
plt.ylabel('Height (km)')
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(5, 12, r'a)', fontdict=font, color="black", fontsize='x-large')

plt.title(r'$\Delta\theta$ (K) at y=50 km ', fontsize='x-large')

plt.subplot(322)
y=hein[1:85]
x=length_y[:]/1000
X, Y = np.meshgrid(x, y)
CS = plt.contourf(X, Y ,qvzy24h_snap_CTRL_f[1:85,:], levels=[-1.,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0], cmap="RdBu_r")
cbar = plt.colorbar()

plt.xlabel('X (km)')
plt.ylabel('Height (km)')
plt.title(r'$\Delta$q$_v$ (g/kg) at y=50 km  ', fontsize='x-large')
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(5, 12, r'b)', fontdict=font, color="black", fontsize='x-large')



plt.subplot(323)
y=length_y[:]/1000
x=length_x[:]/1000
X, Y = np.meshgrid(x, y)
CS = plt.contourf(X, Y ,thxy24h_1km_snap_CTRL_f[:,:], levels=[-0.3,-0.2,-0.1,0,0.1,0.2,0.3], cmap="RdBu_r")
cbar = plt.colorbar()
plt.xlabel('X (km)')
plt.ylabel('Y (km)')

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(5, 90, r'c)', fontdict=font, color="black", fontsize='x-large')
plt.title(r'$\Delta\theta$ (K) at z=1 km ', fontsize='x-large')

plt.subplot(324)
y=length_y[:]/1000
x=length_x[:]/1000
X, Y = np.meshgrid(x, y)
CS = plt.contourf(X, Y ,1000.*qvxy24h_1km_snap_CTRL_f[:,:], levels=[-1.,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0], cmap="RdBu_r")
cbar = plt.colorbar()
plt.xlabel('X (km)')
plt.ylabel('Y (km)')

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(5, 90, r'd)', fontdict=font, color="black", fontsize='x-large')
plt.title(r'$\Delta$q$_v$ (g/kg) at z=1 km ', fontsize='x-large')


plt.subplot(325)
y=length_y[:]/1000
x=length_x[:]/1000
X, Y = np.meshgrid(x, y)
CS = plt.contourf(X, Y ,thxy24h_3km_snap_CTRL_f[:,:], levels=[-0.3,-0.2,-0.1,0,0.1,0.2,0.3], cmap="RdBu_r")
cbar = plt.colorbar()
plt.xlabel('X (km)')
plt.ylabel('Y (km)')

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(5, 90, r'e)', fontdict=font, color="black", fontsize='x-large')
plt.title(r'$\Delta\theta$ (K) at z=3 km ', fontsize='x-large')

plt.subplot(326)
y=length_y[:]/1000
x=length_x[:]/1000
X, Y = np.meshgrid(x, y)
CS = plt.contourf(X, Y ,1000.*qvxy24h_3km_snap_CTRL_f[:,:], levels=[-1.,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0], cmap="RdBu_r")
cbar = plt.colorbar()
plt.xlabel('X (km)')
plt.ylabel('Y (km)')

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(5, 90, r'f)', fontdict=font, color="black", fontsize='x-large')
plt.title(r'$\Delta$q$_v$ (g/kg) at z=3 km', fontsize='x-large')

plt.savefig('Figure7.png')









fig = plt.figure(figsize=(16,20))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
plt.subplot(311)
plt.plot(time_10min_CTRL_f_ens_mean_R[0:100], precip_CTRL_f_ens_mean_R[0:100], lw=4, color="black", label='C')
plt.fill_between(time_10min_CTRL_f_ens_mean_R[0:100], precip_CTRL_f_ens_meanM_R[0:100], precip_CTRL_f_ens_meanP_R[0:100], alpha=0.5, edgecolor='grey', facecolor='grey')

plt.plot(time_10min_RTHQV_ens_mean_R[0:100], precip_RTHQV_ens_mean_R[0:100], lw=4, color="red", label='H')
plt.fill_between(time_10min_RTHQV_ens_mean_R[0:100], precip_RTHQV_ens_meanM_R[0:100], precip_RTHQV_ens_meanP_R[0:100], alpha=0.5, edgecolor='lightcoral', facecolor='lightcoral')
plt.legend()
plt.plot(time_SF[0:100], (surf_flx_CTRL_f[0:100])/1000.0, lw=1,linestyle=':', color="black", label='(shf+lhf)/e3')


plt.ylabel('Precipitation (mm/h)', fontsize='x-large')
plt.xlabel('Time (hour) ', fontsize='x-large')
plt.xlim(0,15)
axes = plt.gca()
axes.set_xlim([0, 15])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])

plt.ylim(0,1.2)
axes = plt.gca()
axes.set_ylim([0, 1.2])
axes.set_yticks([0.,0.3,0.6,0.9,1.2])
axes.set_yticklabels([0.,0.3,0.6,0.9,1.2])
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 1.1, r'a)', fontdict=font, color="black", fontsize='x-large')


plt.subplot(312)
plt.plot(time_SF[0:100], (surf_flx_CTRL_f[0:100])/10000.0, lw=1,linestyle=':', color="black", label='(shf+lhf)/e4')
plt.plot(time_10min_CTRL_f_ens_mean_R[0:100], cld_MF_cb_CTRL_f_ens_mean_R[0:100], lw=4, color="black", label='Ens mean')
plt.fill_between(time_10min_CTRL_f_ens_mean_R[0:100], cld_MF_cb_CTRL_f_ens_meanM_R[0:100], cld_MF_cb_CTRL_f_ens_meanP_R[0:100], alpha=0.5, edgecolor='grey', facecolor='grey')


plt.plot(time_10min_RTHQV_ens_mean_R[0:100], cld_MF_cb_RTHQV_ens_mean_R[0:100], lw=4, color="red", label='Ens mean')
plt.fill_between(time_10min_RTHQV_ens_mean_R[0:100], cld_MF_cb_RTHQV_ens_meanM_R[0:100], cld_MF_cb_RTHQV_ens_meanP_R[0:100], alpha=0.5, edgecolor='lightcoral', facecolor='lightcoral')

plt.ylabel('Mass Flux (kg/m$^2$/s)', fontsize='x-large')
plt.xlabel('Time (hour) ', fontsize='x-large')
plt.xlim(0,15)
axes = plt.gca()
axes.set_xlim([0, 15])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])

plt.ylim(0,0.15)
axes = plt.gca()
axes.set_ylim([0, 0.15])
axes.set_yticks([0,0.05,0.1,0.15])
axes.set_yticklabels([0,0.05,0.1,0.15])
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 0.13, r'b)', fontdict=font, color="black", fontsize='x-large')

plt.savefig('Figure8.png')











