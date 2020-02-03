from mpl_toolkits.basemap import Basemap, cm
from scipy import interpolate
import numpy as N 
import matplotlib.pyplot as plt
import os 
import netCDF4
from netCDF4 import Dataset
import numpy
from scipy.stats import pearsonr
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab


hour_ctrl_mean_shit0=np.load('hour_ctrl_mean_shit0.npy')
hour_phalf_mean_shit0=np.load('hour_phalf_mean_shit0.npy')
hour_mhalf_mean_shit0=np.load('hour_mhalf_mean_shit0.npy')
hour_rthqv_mean_shit0=np.load('hour_rthqv_mean_shit0.npy')

area_freq_ctrlD2=np.load('area_freq_ctrlD2.npy')
area_freq_ctrlD3=np.load('area_freq_ctrlD3.npy')
area_freq_ctrlD4=np.load('area_freq_ctrlD4.npy')
area_freq_ctrlD5=np.load('area_freq_ctrlD5.npy')
area_freq_ctrlD6=np.load('area_freq_ctrlD6.npy')
area_freq_ctrlD7=np.load('area_freq_ctrlD7.npy')
area_freq_ctrlD8=np.load('area_freq_ctrlD8.npy')
area_freq_ctrlD9=np.load('area_freq_ctrlD9.npy')
area_freq_ctrlD10=np.load('area_freq_ctrlD10.npy')



nhr=33
nd=11

nhist=300

area_hist=N.zeros((nhist))
for i in N.arange(nhist-1):
     area_hist[i+1]=area_hist[i]+1

area_freq_ctrl_ens_gs1_mean=N.zeros((nhist))
area_freq_ctrl_ens_gs2_mean=N.zeros((nhist))
area_freq_ctrl_ens_gs3_mean=N.zeros((nhist))
area_freq_ctrl_ens_gs4_mean=N.zeros((nhist))
area_freq_ctrl_ens_ms_mean=N.zeros((nhist))
area_freq_ctrl_ens_ss_mean=N.zeros((nhist))
area_freq_ctrl_ens_gs1_meanM=N.zeros((nhist))
area_freq_ctrl_ens_gs2_meanM=N.zeros((nhist))
area_freq_ctrl_ens_gs3_meanM=N.zeros((nhist))
area_freq_ctrl_ens_gs4_meanM=N.zeros((nhist))
area_freq_ctrl_ens_ms_meanM=N.zeros((nhist))
area_freq_ctrl_ens_ss_meanM=N.zeros((nhist))
area_freq_ctrl_ens_gs1_meanP=N.zeros((nhist))
area_freq_ctrl_ens_gs2_meanP=N.zeros((nhist))
area_freq_ctrl_ens_gs3_meanP=N.zeros((nhist))
area_freq_ctrl_ens_gs4_meanP=N.zeros((nhist))
area_freq_ctrl_ens_ms_meanP=N.zeros((nhist))
area_freq_ctrl_ens_ss_meanP=N.zeros((nhist))

area_freq_ctrl_ens_gs1=N.zeros((nhist,11))
area_freq_ctrl_ens_gs2=N.zeros((nhist,11))
area_freq_ctrl_ens_gs3=N.zeros((nhist,11))
area_freq_ctrl_ens_gs4=N.zeros((nhist,11))
area_freq_ctrl_ens_ms=N.zeros((nhist,11))
area_freq_ctrl_ens_ss=N.zeros((nhist,11))

area_freq_ctrl_ens_gs1[1:nhist,2]=(area_freq_ctrlD2[3,1:nhist]+area_freq_ctrlD2[4,1:nhist]+area_freq_ctrlD2[5,1:nhist])/3
area_freq_ctrl_ens_gs2[1:nhist,2]=(area_freq_ctrlD2[6,1:nhist]+area_freq_ctrlD2[7,1:nhist])/2
area_freq_ctrl_ens_gs3[1:nhist,2]= (area_freq_ctrlD2[8,1:nhist]+area_freq_ctrlD2[9,1:nhist])/2
area_freq_ctrl_ens_gs4[1:nhist,2]=  (area_freq_ctrlD2[10,1:nhist]+area_freq_ctrlD2[11,1:nhist])/2
area_freq_ctrl_ens_ms[1:nhist,2]=(area_freq_ctrlD2[12,1:nhist]+area_freq_ctrlD2[13,1:nhist]+area_freq_ctrlD2[14,1:nhist]+area_freq_ctrlD2[15,1:nhist]+area_freq_ctrlD2[16,1:nhist]+area_freq_ctrlD2[17,1:nhist]+area_freq_ctrlD2[18,1:nhist]+area_freq_ctrlD2[19,1:nhist]+area_freq_ctrlD2[20,1:nhist]+area_freq_ctrlD2[21,1:nhist]+area_freq_ctrlD2[22,1:nhist]+area_freq_ctrlD2[23,1:nhist]+area_freq_ctrlD2[24,1:nhist])/13.0
area_freq_ctrl_ens_ss[1:nhist,2]=(area_freq_ctrlD2[25,1:nhist]+area_freq_ctrlD2[26,1:nhist]+area_freq_ctrlD2[27,1:nhist]+area_freq_ctrlD2[28,1:nhist])/4.0

area_freq_ctrl_ens_gs1[1:nhist,3]=(area_freq_ctrlD3[3,1:nhist]+area_freq_ctrlD3[4,1:nhist]+area_freq_ctrlD3[5,1:nhist])/3
area_freq_ctrl_ens_gs2[1:nhist,3]=(area_freq_ctrlD3[6,1:nhist]+area_freq_ctrlD3[7,1:nhist])/2
area_freq_ctrl_ens_gs3[1:nhist,3]= (area_freq_ctrlD3[8,1:nhist]+area_freq_ctrlD3[9,1:nhist])/2
area_freq_ctrl_ens_gs4[1:nhist,3]=  (area_freq_ctrlD3[10,1:nhist]+area_freq_ctrlD3[11,1:nhist])/2
area_freq_ctrl_ens_ms[1:nhist,3]=(area_freq_ctrlD3[12,1:nhist]+area_freq_ctrlD3[13,1:nhist]+area_freq_ctrlD3[14,1:nhist]+area_freq_ctrlD3[15,1:nhist]+area_freq_ctrlD3[16,1:nhist]+area_freq_ctrlD3[17,1:nhist]+area_freq_ctrlD3[18,1:nhist]+area_freq_ctrlD3[19,1:nhist]+area_freq_ctrlD3[20,1:nhist]+area_freq_ctrlD3[21,1:nhist]+area_freq_ctrlD3[22,1:nhist]+area_freq_ctrlD3[23,1:nhist]+area_freq_ctrlD3[24,1:nhist])/13.0
area_freq_ctrl_ens_ss[1:nhist,3]=(area_freq_ctrlD3[25,1:nhist]+area_freq_ctrlD3[26,1:nhist]+area_freq_ctrlD3[27,1:nhist]+area_freq_ctrlD3[28,1:nhist])/4.0

area_freq_ctrl_ens_gs1[1:nhist,4]=(area_freq_ctrlD4[3,1:nhist]+area_freq_ctrlD4[4,1:nhist]+area_freq_ctrlD4[5,1:nhist])/3
area_freq_ctrl_ens_gs2[1:nhist,4]=(area_freq_ctrlD4[6,1:nhist]+area_freq_ctrlD4[7,1:nhist])/2
area_freq_ctrl_ens_gs3[1:nhist,4]= (area_freq_ctrlD4[8,1:nhist]+area_freq_ctrlD4[9,1:nhist])/2
area_freq_ctrl_ens_gs4[1:nhist,4]=  (area_freq_ctrlD4[10,1:nhist]+area_freq_ctrlD4[11,1:nhist])/2
area_freq_ctrl_ens_ms[1:nhist,4]=(area_freq_ctrlD4[12,1:nhist]+area_freq_ctrlD4[13,1:nhist]+area_freq_ctrlD4[14,1:nhist]+area_freq_ctrlD4[15,1:nhist]+area_freq_ctrlD4[16,1:nhist]+area_freq_ctrlD4[17,1:nhist]+area_freq_ctrlD4[18,1:nhist]+area_freq_ctrlD4[19,1:nhist]+area_freq_ctrlD4[20,1:nhist]+area_freq_ctrlD4[21,1:nhist]+area_freq_ctrlD4[22,1:nhist]+area_freq_ctrlD4[23,1:nhist]+area_freq_ctrlD4[24,1:nhist])/13.0
area_freq_ctrl_ens_ss[1:nhist,4]=(area_freq_ctrlD4[25,1:nhist]+area_freq_ctrlD4[26,1:nhist]+area_freq_ctrlD4[27,1:nhist]+area_freq_ctrlD4[28,1:nhist])/4.0

area_freq_ctrl_ens_gs1[1:nhist,5]=(area_freq_ctrlD5[3,1:nhist]+area_freq_ctrlD5[4,1:nhist]+area_freq_ctrlD5[5,1:nhist])/3
area_freq_ctrl_ens_gs2[1:nhist,5]=(area_freq_ctrlD5[6,1:nhist]+area_freq_ctrlD5[7,1:nhist])/2
area_freq_ctrl_ens_gs3[1:nhist,5]= (area_freq_ctrlD5[8,1:nhist]+area_freq_ctrlD5[9,1:nhist])/2
area_freq_ctrl_ens_gs4[1:nhist,5]=  (area_freq_ctrlD5[10,1:nhist]+area_freq_ctrlD5[11,1:nhist])/2
area_freq_ctrl_ens_ms[1:nhist,5]=(area_freq_ctrlD5[12,1:nhist]+area_freq_ctrlD5[13,1:nhist]+area_freq_ctrlD5[14,1:nhist]+area_freq_ctrlD5[15,1:nhist]+area_freq_ctrlD5[16,1:nhist]+area_freq_ctrlD5[17,1:nhist]+area_freq_ctrlD5[18,1:nhist]+area_freq_ctrlD5[19,1:nhist]+area_freq_ctrlD5[20,1:nhist]+area_freq_ctrlD5[21,1:nhist]+area_freq_ctrlD5[22,1:nhist]+area_freq_ctrlD5[23,1:nhist]+area_freq_ctrlD5[24,1:nhist])/13.0
area_freq_ctrl_ens_ss[1:nhist,5]=(area_freq_ctrlD5[25,1:nhist]+area_freq_ctrlD5[26,1:nhist]+area_freq_ctrlD5[27,1:nhist]+area_freq_ctrlD5[28,1:nhist])/4.0

area_freq_ctrl_ens_gs1[1:nhist,6]=(area_freq_ctrlD6[3,1:nhist]+area_freq_ctrlD6[4,1:nhist]+area_freq_ctrlD6[5,1:nhist])/3
area_freq_ctrl_ens_gs2[1:nhist,6]=(area_freq_ctrlD6[6,1:nhist]+area_freq_ctrlD6[7,1:nhist])/2
area_freq_ctrl_ens_gs3[1:nhist,6]= (area_freq_ctrlD6[8,1:nhist]+area_freq_ctrlD6[9,1:nhist])/2
area_freq_ctrl_ens_gs4[1:nhist,6]=  (area_freq_ctrlD6[10,1:nhist]+area_freq_ctrlD6[11,1:nhist])/2
area_freq_ctrl_ens_ms[1:nhist,6]=(area_freq_ctrlD6[12,1:nhist]+area_freq_ctrlD6[13,1:nhist]+area_freq_ctrlD6[14,1:nhist]+area_freq_ctrlD6[15,1:nhist]+area_freq_ctrlD6[16,1:nhist]+area_freq_ctrlD6[17,1:nhist]+area_freq_ctrlD6[18,1:nhist]+area_freq_ctrlD6[19,1:nhist]+area_freq_ctrlD6[20,1:nhist]+area_freq_ctrlD6[21,1:nhist]+area_freq_ctrlD6[22,1:nhist]+area_freq_ctrlD6[23,1:nhist]+area_freq_ctrlD6[24,1:nhist])/13.0
area_freq_ctrl_ens_ss[1:nhist,6]=(area_freq_ctrlD6[25,1:nhist]+area_freq_ctrlD6[26,1:nhist]+area_freq_ctrlD6[27,1:nhist]+area_freq_ctrlD6[28,1:nhist])/4.0

area_freq_ctrl_ens_gs1[1:nhist,7]=(area_freq_ctrlD7[3,1:nhist]+area_freq_ctrlD7[4,1:nhist]+area_freq_ctrlD7[5,1:nhist])/3
area_freq_ctrl_ens_gs2[1:nhist,7]=(area_freq_ctrlD7[6,1:nhist]+area_freq_ctrlD7[7,1:nhist])/2
area_freq_ctrl_ens_gs3[1:nhist,7]= (area_freq_ctrlD7[8,1:nhist]+area_freq_ctrlD7[9,1:nhist])/2
area_freq_ctrl_ens_gs4[1:nhist,7]=  (area_freq_ctrlD7[10,1:nhist]+area_freq_ctrlD7[11,1:nhist])/2
area_freq_ctrl_ens_ms[1:nhist,7]=(area_freq_ctrlD7[12,1:nhist]+area_freq_ctrlD7[13,1:nhist]+area_freq_ctrlD7[14,1:nhist]+area_freq_ctrlD7[15,1:nhist]+area_freq_ctrlD7[16,1:nhist]+area_freq_ctrlD7[17,1:nhist]+area_freq_ctrlD7[18,1:nhist]+area_freq_ctrlD7[19,1:nhist]+area_freq_ctrlD7[20,1:nhist]+area_freq_ctrlD7[21,1:nhist]+area_freq_ctrlD7[22,1:nhist]+area_freq_ctrlD7[23,1:nhist]+area_freq_ctrlD7[24,1:nhist])/13.0
area_freq_ctrl_ens_ss[1:nhist,7]=(area_freq_ctrlD7[25,1:nhist]+area_freq_ctrlD7[26,1:nhist]+area_freq_ctrlD7[27,1:nhist]+area_freq_ctrlD7[28,1:nhist])/4.0

area_freq_ctrl_ens_gs1[1:nhist,8]=(area_freq_ctrlD8[3,1:nhist]+area_freq_ctrlD8[4,1:nhist]+area_freq_ctrlD8[5,1:nhist])/3
area_freq_ctrl_ens_gs2[1:nhist,8]=(area_freq_ctrlD8[6,1:nhist]+area_freq_ctrlD8[7,1:nhist])/2
area_freq_ctrl_ens_gs3[1:nhist,8]= (area_freq_ctrlD8[8,1:nhist]+area_freq_ctrlD8[9,1:nhist])/2
area_freq_ctrl_ens_gs4[1:nhist,8]=  (area_freq_ctrlD8[10,1:nhist]+area_freq_ctrlD8[11,1:nhist])/2
area_freq_ctrl_ens_ms[1:nhist,8]=(area_freq_ctrlD8[12,1:nhist]+area_freq_ctrlD8[13,1:nhist]+area_freq_ctrlD8[14,1:nhist]+area_freq_ctrlD8[15,1:nhist]+area_freq_ctrlD8[16,1:nhist]+area_freq_ctrlD8[17,1:nhist]+area_freq_ctrlD8[18,1:nhist]+area_freq_ctrlD8[19,1:nhist]+area_freq_ctrlD8[20,1:nhist]+area_freq_ctrlD8[21,1:nhist]+area_freq_ctrlD8[22,1:nhist]+area_freq_ctrlD8[23,1:nhist]+area_freq_ctrlD8[24,1:nhist])/13.0
area_freq_ctrl_ens_ss[1:nhist,8]=(area_freq_ctrlD8[25,1:nhist]+area_freq_ctrlD8[26,1:nhist]+area_freq_ctrlD8[27,1:nhist]+area_freq_ctrlD8[28,1:nhist])/4.0

area_freq_ctrl_ens_gs1[1:nhist,9]=(area_freq_ctrlD9[3,1:nhist]+area_freq_ctrlD9[4,1:nhist]+area_freq_ctrlD9[5,1:nhist])/3
area_freq_ctrl_ens_gs2[1:nhist,9]=(area_freq_ctrlD9[6,1:nhist]+area_freq_ctrlD9[7,1:nhist])/2
area_freq_ctrl_ens_gs3[1:nhist,9]= (area_freq_ctrlD9[8,1:nhist]+area_freq_ctrlD9[9,1:nhist])/2
area_freq_ctrl_ens_gs4[1:nhist,9]=  (area_freq_ctrlD9[10,1:nhist]+area_freq_ctrlD9[11,1:nhist])/2
area_freq_ctrl_ens_ms[1:nhist,9]=(area_freq_ctrlD9[12,1:nhist]+area_freq_ctrlD9[13,1:nhist]+area_freq_ctrlD9[14,1:nhist]+area_freq_ctrlD9[15,1:nhist]+area_freq_ctrlD9[16,1:nhist]+area_freq_ctrlD9[17,1:nhist]+area_freq_ctrlD9[18,1:nhist]+area_freq_ctrlD9[19,1:nhist]+area_freq_ctrlD9[20,1:nhist]+area_freq_ctrlD9[21,1:nhist]+area_freq_ctrlD9[22,1:nhist]+area_freq_ctrlD9[23,1:nhist]+area_freq_ctrlD9[24,1:nhist])/13.0
area_freq_ctrl_ens_ss[1:nhist,9]=(area_freq_ctrlD9[25,1:nhist]+area_freq_ctrlD9[26,1:nhist]+area_freq_ctrlD9[27,1:nhist]+area_freq_ctrlD9[28,1:nhist])/4.0

area_freq_ctrl_ens_gs1[1:nhist,10]=(area_freq_ctrlD10[3,1:nhist]+area_freq_ctrlD10[4,1:nhist]+area_freq_ctrlD10[5,1:nhist])/3
area_freq_ctrl_ens_gs2[1:nhist,10]=(area_freq_ctrlD10[6,1:nhist]+area_freq_ctrlD10[7,1:nhist])/2
area_freq_ctrl_ens_gs3[1:nhist,10]= (area_freq_ctrlD10[8,1:nhist]+area_freq_ctrlD10[9,1:nhist])/2
area_freq_ctrl_ens_gs4[1:nhist,10]=  (area_freq_ctrlD10[10,1:nhist]+area_freq_ctrlD10[11,1:nhist])/2
area_freq_ctrl_ens_ms[1:nhist,10]=(area_freq_ctrlD10[12,1:nhist]+area_freq_ctrlD10[13,1:nhist]+area_freq_ctrlD10[14,1:nhist]+area_freq_ctrlD10[15,1:nhist]+area_freq_ctrlD10[16,1:nhist]+area_freq_ctrlD10[17,1:nhist]+area_freq_ctrlD10[18,1:nhist]+area_freq_ctrlD10[19,1:nhist]+area_freq_ctrlD10[20,1:nhist]+area_freq_ctrlD10[21,1:nhist]+area_freq_ctrlD10[22,1:nhist]+area_freq_ctrlD10[23,1:nhist]+area_freq_ctrlD10[24,1:nhist])/13.0
area_freq_ctrl_ens_ss[1:nhist,10]=(area_freq_ctrlD10[25,1:nhist]+area_freq_ctrlD10[26,1:nhist]+area_freq_ctrlD10[27,1:nhist]+area_freq_ctrlD10[28,1:nhist])/4.0

for ih in N.arange(nhist):
     area_freq_ctrl_ens_gs1_mean[ih]=area_freq_ctrl_ens_gs1[ih,2:11].mean()
     area_freq_ctrl_ens_gs1_meanM[ih]=area_freq_ctrl_ens_gs1_mean[ih]-np.std(area_freq_ctrl_ens_gs1[ih,2:11])
     area_freq_ctrl_ens_gs1_meanP[ih]=area_freq_ctrl_ens_gs1_mean[ih]+np.std(area_freq_ctrl_ens_gs1[ih,2:11])
     area_freq_ctrl_ens_gs2_mean[ih]=area_freq_ctrl_ens_gs2[ih,2:11].mean()
     area_freq_ctrl_ens_gs2_meanM[ih]=area_freq_ctrl_ens_gs2_mean[ih]-np.std(area_freq_ctrl_ens_gs2[ih,2:11])
     area_freq_ctrl_ens_gs2_meanP[ih]=area_freq_ctrl_ens_gs2_mean[ih]+np.std(area_freq_ctrl_ens_gs2[ih,2:11])
     area_freq_ctrl_ens_gs3_mean[ih]=area_freq_ctrl_ens_gs3[ih,2:11].mean()
     area_freq_ctrl_ens_gs3_meanM[ih]=area_freq_ctrl_ens_gs3_mean[ih]-np.std(area_freq_ctrl_ens_gs3[ih,2:11])
     area_freq_ctrl_ens_gs3_meanP[ih]=area_freq_ctrl_ens_gs3_mean[ih]+np.std(area_freq_ctrl_ens_gs3[ih,2:11])
     area_freq_ctrl_ens_gs4_mean[ih]=area_freq_ctrl_ens_gs4[ih,2:11].mean()
     area_freq_ctrl_ens_gs4_meanM[ih]=area_freq_ctrl_ens_gs4_mean[ih]-np.std(area_freq_ctrl_ens_gs4[ih,2:11])
     area_freq_ctrl_ens_gs4_meanP[ih]=area_freq_ctrl_ens_gs4_mean[ih]+np.std(area_freq_ctrl_ens_gs4[ih,2:11])

     area_freq_ctrl_ens_ms_mean[ih]=area_freq_ctrl_ens_ms[ih,2:11].mean()
     area_freq_ctrl_ens_ms_meanM[ih]=area_freq_ctrl_ens_ms_mean[ih]-np.std(area_freq_ctrl_ens_ms[ih,2:11])
     area_freq_ctrl_ens_ms_meanP[ih]=area_freq_ctrl_ens_ms_mean[ih]+np.std(area_freq_ctrl_ens_ms[ih,2:11])

     area_freq_ctrl_ens_ss_mean[ih]=area_freq_ctrl_ens_ss[ih,2:11].mean()
     area_freq_ctrl_ens_ss_meanM[ih]=area_freq_ctrl_ens_ss_mean[ih]-np.std(area_freq_ctrl_ens_ss[ih,2:11])
     area_freq_ctrl_ens_ss_meanP[ih]=area_freq_ctrl_ens_ss_mean[ih]+np.std(area_freq_ctrl_ens_ss[ih,2:11])



cloud_stat_phalfD10=np.load('cloud_stat_phalfD10.npy')
cloud_stat_phalfD9=np.load('cloud_stat_phalfD9.npy')
cloud_stat_phalfD8=np.load('cloud_stat_phalfD8.npy')
cloud_stat_phalfD7=np.load('cloud_stat_phalfD7.npy')
cloud_stat_phalfD6=np.load('cloud_stat_phalfD6.npy')
cloud_stat_phalfD5=np.load('cloud_stat_phalfD5.npy')
cloud_stat_phalfD4=np.load('cloud_stat_phalfD4.npy')
cloud_stat_phalfD3=np.load('cloud_stat_phalfD3.npy')
cloud_stat_phalfD2=np.load('cloud_stat_phalfD2.npy')

meanradius_phalf_ens=N.zeros((nhr,11))
meanradius_phalf_ens[:,0]=0.0
meanradius_phalf_ens[:,1]=0.0
meanradius_phalf_ens[:,2]=cloud_stat_phalfD2[:,11]
meanradius_phalf_ens[:,3]=cloud_stat_phalfD3[:,11]
meanradius_phalf_ens[:,4]=cloud_stat_phalfD4[:,11]
meanradius_phalf_ens[:,5]=cloud_stat_phalfD5[:,11]
meanradius_phalf_ens[:,6]=cloud_stat_phalfD6[:,11]
meanradius_phalf_ens[:,7]=cloud_stat_phalfD7[:,11]
meanradius_phalf_ens[:,8]=cloud_stat_phalfD8[:,11]
meanradius_phalf_ens[:,9]=cloud_stat_phalfD9[:,11]
meanradius_phalf_ens[:,10]=cloud_stat_phalfD10[:,11]

meanarea_phalf_ens=N.zeros((nhr,11))
meanarea_phalf_ens[:,0]=0.0
meanarea_phalf_ens[:,1]=0.0
meanarea_phalf_ens[:,2]=cloud_stat_phalfD2[:,2]
meanarea_phalf_ens[:,3]=cloud_stat_phalfD3[:,2]
meanarea_phalf_ens[:,4]=cloud_stat_phalfD4[:,2]
meanarea_phalf_ens[:,5]=cloud_stat_phalfD5[:,2]
meanarea_phalf_ens[:,6]=cloud_stat_phalfD6[:,2]
meanarea_phalf_ens[:,7]=cloud_stat_phalfD7[:,2]
meanarea_phalf_ens[:,8]=cloud_stat_phalfD8[:,2]
meanarea_phalf_ens[:,9]=cloud_stat_phalfD9[:,2]
meanarea_phalf_ens[:,10]=cloud_stat_phalfD10[:,2]

medianarea_phalf_ens=N.zeros((nhr,11))
medianarea_phalf_ens[:,0]=0.0
medianarea_phalf_ens[:,1]=0.0
medianarea_phalf_ens[:,2]=cloud_stat_phalfD2[:,7]
medianarea_phalf_ens[:,3]=cloud_stat_phalfD3[:,7]
medianarea_phalf_ens[:,4]=cloud_stat_phalfD4[:,7]
medianarea_phalf_ens[:,5]=cloud_stat_phalfD5[:,7]
medianarea_phalf_ens[:,6]=cloud_stat_phalfD6[:,7]
medianarea_phalf_ens[:,7]=cloud_stat_phalfD7[:,7]
medianarea_phalf_ens[:,8]=cloud_stat_phalfD8[:,7]
medianarea_phalf_ens[:,9]=cloud_stat_phalfD9[:,7]
medianarea_phalf_ens[:,10]=cloud_stat_phalfD10[:,7]


num_phalf_ens=N.zeros((nhr,11))
num_phalf_ens[:,0]=0.0
num_phalf_ens[:,1]=0.0
num_phalf_ens[:,2]=cloud_stat_phalfD2[:,1]
num_phalf_ens[:,3]=cloud_stat_phalfD3[:,1]
num_phalf_ens[:,4]=cloud_stat_phalfD4[:,1]
num_phalf_ens[:,5]=cloud_stat_phalfD5[:,1]
num_phalf_ens[:,6]=cloud_stat_phalfD6[:,1]
num_phalf_ens[:,7]=cloud_stat_phalfD7[:,1]
num_phalf_ens[:,8]=cloud_stat_phalfD8[:,1]
num_phalf_ens[:,9]=cloud_stat_phalfD9[:,1]
num_phalf_ens[:,10]=cloud_stat_phalfD10[:,1]

stdvR_phalf_ens=N.zeros((nhr,11))
stdvR_phalf_ens[:,0]=0.0
stdvR_phalf_ens[:,1]=0.0
stdvR_phalf_ens[:,2]=cloud_stat_phalfD2[:,3]
stdvR_phalf_ens[:,3]=cloud_stat_phalfD3[:,3]
stdvR_phalf_ens[:,4]=cloud_stat_phalfD4[:,3]
stdvR_phalf_ens[:,5]=cloud_stat_phalfD5[:,3]
stdvR_phalf_ens[:,6]=cloud_stat_phalfD6[:,3]
stdvR_phalf_ens[:,7]=cloud_stat_phalfD7[:,3]
stdvR_phalf_ens[:,8]=cloud_stat_phalfD8[:,3]
stdvR_phalf_ens[:,9]=cloud_stat_phalfD9[:,3]
stdvR_phalf_ens[:,10]=cloud_stat_phalfD10[:,3]

num_phalf_ens_mean=N.zeros((nhr))
num_phalf_ens_meanP=N.zeros((nhr))
num_phalf_ens_meanM=N.zeros((nhr))
meanradius_phalf_ens_mean=N.zeros((nhr))
medianarea_phalf_ens_mean=N.zeros((nhr))
meanradius_phalf_ens_meanP=N.zeros((nhr))
meanradius_phalf_ens_meanM=N.zeros((nhr))
meanarea_phalf_ens_mean=N.zeros((nhr))
meanarea_phalf_ens_meanP=N.zeros((nhr))
meanarea_phalf_ens_meanM=N.zeros((nhr))
stdvR_phalf_ens_mean=N.zeros((nhr))
stdvR_phalf_ens_meanP=N.zeros((nhr))
stdvR_phalf_ens_meanM=N.zeros((nhr))
for j in N.arange(nhr):
     num_phalf_ens_mean[j]=num_phalf_ens[j,2:11].mean()
     num_phalf_ens_meanM[j]=num_phalf_ens_mean[j]-np.std(num_phalf_ens[j,2:11])
     num_phalf_ens_meanP[j]=num_phalf_ens_mean[j]+np.std(num_phalf_ens[j,2:11])
     meanradius_phalf_ens_mean[j]=meanradius_phalf_ens[j,2:11].mean()
     meanradius_phalf_ens_meanM[j]=meanradius_phalf_ens_mean[j]-np.std(meanradius_phalf_ens[j,2:11])
     meanradius_phalf_ens_meanP[j]=meanradius_phalf_ens_mean[j]+np.std(meanradius_phalf_ens[j,2:11])

     meanarea_phalf_ens_mean[j]=meanarea_phalf_ens[j,2:11].mean()
     meanarea_phalf_ens_meanM[j]=meanarea_phalf_ens_mean[j]-np.std(meanarea_phalf_ens[j,2:11])
     meanarea_phalf_ens_meanP[j]=meanarea_phalf_ens_mean[j]+np.std(meanarea_phalf_ens[j,2:11])
     medianarea_phalf_ens_mean[j]=medianarea_phalf_ens[j,2:11].mean()

     stdvR_phalf_ens_mean[j]=stdvR_phalf_ens[j,2:11].mean()
     stdvR_phalf_ens_meanM[j]=stdvR_phalf_ens_mean[j]-np.std(stdvR_phalf_ens[j,2:11])#/3.0
     stdvR_phalf_ens_meanP[j]=stdvR_phalf_ens_mean[j]+np.std(stdvR_phalf_ens[j,2:11])#/3.0


cloud_stat_mhalfD10=np.load('cloud_stat_mhalfD10.npy')
cloud_stat_mhalfD9=np.load('cloud_stat_mhalfD9.npy')
cloud_stat_mhalfD8=np.load('cloud_stat_mhalfD8.npy')
cloud_stat_mhalfD7=np.load('cloud_stat_mhalfD7.npy')
cloud_stat_mhalfD6=np.load('cloud_stat_mhalfD6.npy')
cloud_stat_mhalfD5=np.load('cloud_stat_mhalfD5.npy')
cloud_stat_mhalfD4=np.load('cloud_stat_mhalfD4.npy')
cloud_stat_mhalfD3=np.load('cloud_stat_mhalfD3.npy')
cloud_stat_mhalfD2=np.load('cloud_stat_mhalfD2.npy')

meanradius_mhalf_ens=N.zeros((nhr,11))
meanradius_mhalf_ens[:,0]=0.0
meanradius_mhalf_ens[:,1]=0.0
meanradius_mhalf_ens[:,2]=cloud_stat_mhalfD2[:,11]
meanradius_mhalf_ens[:,3]=cloud_stat_mhalfD3[:,11]
meanradius_mhalf_ens[:,4]=cloud_stat_mhalfD4[:,11]
meanradius_mhalf_ens[:,5]=cloud_stat_mhalfD5[:,11]
meanradius_mhalf_ens[:,6]=cloud_stat_mhalfD6[:,11]
meanradius_mhalf_ens[:,7]=cloud_stat_mhalfD7[:,11]
meanradius_mhalf_ens[:,8]=cloud_stat_mhalfD8[:,11]
meanradius_mhalf_ens[:,9]=cloud_stat_mhalfD9[:,11]
meanradius_mhalf_ens[:,10]=cloud_stat_mhalfD10[:,11]

meanarea_mhalf_ens=N.zeros((nhr,11))
meanarea_mhalf_ens[:,0]=0.0
meanarea_mhalf_ens[:,1]=0.0
meanarea_mhalf_ens[:,2]=cloud_stat_mhalfD2[:,2]
meanarea_mhalf_ens[:,3]=cloud_stat_mhalfD3[:,2]
meanarea_mhalf_ens[:,4]=cloud_stat_mhalfD4[:,2]
meanarea_mhalf_ens[:,5]=cloud_stat_mhalfD5[:,2]
meanarea_mhalf_ens[:,6]=cloud_stat_mhalfD6[:,2]
meanarea_mhalf_ens[:,7]=cloud_stat_mhalfD7[:,2]
meanarea_mhalf_ens[:,8]=cloud_stat_mhalfD8[:,2]
meanarea_mhalf_ens[:,9]=cloud_stat_mhalfD9[:,2]
meanarea_mhalf_ens[:,10]=cloud_stat_mhalfD10[:,2]

medianarea_mhalf_ens=N.zeros((nhr,11))
medianarea_mhalf_ens[:,0]=0.0
medianarea_mhalf_ens[:,1]=0.0
medianarea_mhalf_ens[:,2]=cloud_stat_mhalfD2[:,7]
medianarea_mhalf_ens[:,3]=cloud_stat_mhalfD3[:,7]
medianarea_mhalf_ens[:,4]=cloud_stat_mhalfD4[:,7]
medianarea_mhalf_ens[:,5]=cloud_stat_mhalfD5[:,7]
medianarea_mhalf_ens[:,6]=cloud_stat_mhalfD6[:,7]
medianarea_mhalf_ens[:,7]=cloud_stat_mhalfD7[:,7]
medianarea_mhalf_ens[:,8]=cloud_stat_mhalfD8[:,7]
medianarea_mhalf_ens[:,9]=cloud_stat_mhalfD9[:,7]
medianarea_mhalf_ens[:,10]=cloud_stat_mhalfD10[:,7]


num_mhalf_ens=N.zeros((nhr,11))
num_mhalf_ens[:,0]=0.0
num_mhalf_ens[:,1]=0.0
num_mhalf_ens[:,2]=cloud_stat_mhalfD2[:,1]
num_mhalf_ens[:,3]=cloud_stat_mhalfD3[:,1]
num_mhalf_ens[:,4]=cloud_stat_mhalfD4[:,1]
num_mhalf_ens[:,5]=cloud_stat_mhalfD5[:,1]
num_mhalf_ens[:,6]=cloud_stat_mhalfD6[:,1]
num_mhalf_ens[:,7]=cloud_stat_mhalfD7[:,1]
num_mhalf_ens[:,8]=cloud_stat_mhalfD8[:,1]
num_mhalf_ens[:,9]=cloud_stat_mhalfD9[:,1]
num_mhalf_ens[:,10]=cloud_stat_mhalfD10[:,1]

stdvR_mhalf_ens=N.zeros((nhr,11))
stdvR_mhalf_ens[:,0]=0.0
stdvR_mhalf_ens[:,1]=0.0
stdvR_mhalf_ens[:,2]=cloud_stat_mhalfD2[:,3]
stdvR_mhalf_ens[:,3]=cloud_stat_mhalfD3[:,3]
stdvR_mhalf_ens[:,4]=cloud_stat_mhalfD4[:,3]
stdvR_mhalf_ens[:,5]=cloud_stat_mhalfD5[:,3]
stdvR_mhalf_ens[:,6]=cloud_stat_mhalfD6[:,3]
stdvR_mhalf_ens[:,7]=cloud_stat_mhalfD7[:,3]
stdvR_mhalf_ens[:,8]=cloud_stat_mhalfD8[:,3]
stdvR_mhalf_ens[:,9]=cloud_stat_mhalfD9[:,3]
stdvR_mhalf_ens[:,10]=cloud_stat_mhalfD10[:,3]

num_mhalf_ens_mean=N.zeros((nhr))
num_mhalf_ens_meanP=N.zeros((nhr))
num_mhalf_ens_meanM=N.zeros((nhr))
meanradius_mhalf_ens_mean=N.zeros((nhr))
medianarea_mhalf_ens_mean=N.zeros((nhr))
meanradius_mhalf_ens_meanP=N.zeros((nhr))
meanradius_mhalf_ens_meanM=N.zeros((nhr))
meanarea_mhalf_ens_mean=N.zeros((nhr))
meanarea_mhalf_ens_meanP=N.zeros((nhr))
meanarea_mhalf_ens_meanM=N.zeros((nhr))
stdvR_mhalf_ens_mean=N.zeros((nhr))
stdvR_mhalf_ens_meanP=N.zeros((nhr))
stdvR_mhalf_ens_meanM=N.zeros((nhr))

for j in N.arange(nhr):
     num_mhalf_ens_mean[j]=num_mhalf_ens[j,2:11].mean()
     num_mhalf_ens_meanM[j]=num_mhalf_ens_mean[j]-np.std(num_mhalf_ens[j,2:11])
     num_mhalf_ens_meanP[j]=num_mhalf_ens_mean[j]+np.std(num_mhalf_ens[j,2:11])
     meanradius_mhalf_ens_mean[j]=meanradius_mhalf_ens[j,2:11].mean()
     meanradius_mhalf_ens_meanM[j]=meanradius_mhalf_ens_mean[j]-np.std(meanradius_mhalf_ens[j,2:11])
     meanradius_mhalf_ens_meanP[j]=meanradius_mhalf_ens_mean[j]+np.std(meanradius_mhalf_ens[j,2:11])

     meanarea_mhalf_ens_mean[j]=meanarea_mhalf_ens[j,2:11].mean()
     meanarea_mhalf_ens_meanM[j]=meanarea_mhalf_ens_mean[j]-np.std(meanarea_mhalf_ens[j,2:11])
     meanarea_mhalf_ens_meanP[j]=meanarea_mhalf_ens_mean[j]+np.std(meanarea_mhalf_ens[j,2:11])
     medianarea_mhalf_ens_mean[j]=medianarea_mhalf_ens[j,2:11].mean()

     stdvR_mhalf_ens_mean[j]=stdvR_mhalf_ens[j,2:11].mean()
     stdvR_mhalf_ens_meanM[j]=stdvR_mhalf_ens_mean[j]-np.std(stdvR_mhalf_ens[j,2:11])
     stdvR_mhalf_ens_meanP[j]=stdvR_mhalf_ens_mean[j]+np.std(stdvR_mhalf_ens[j,2:11])

cloud_stat_ctrlD10=np.load('cloud_stat_ctrlD10.npy')
cloud_stat_ctrlD9=np.load('cloud_stat_ctrlD9.npy')
cloud_stat_ctrlD8=np.load('cloud_stat_ctrlD8.npy')
cloud_stat_ctrlD7=np.load('cloud_stat_ctrlD7.npy')
cloud_stat_ctrlD6=np.load('cloud_stat_ctrlD6.npy')
cloud_stat_ctrlD5=np.load('cloud_stat_ctrlD5.npy')
cloud_stat_ctrlD4=np.load('cloud_stat_ctrlD4.npy')
cloud_stat_ctrlD3=np.load('cloud_stat_ctrlD3.npy')
cloud_stat_ctrlD2=np.load('cloud_stat_ctrlD2.npy')

meanradius_ctrl_ens=N.zeros((nhr,11))
meanradius_ctrl_ens[:,0]=0.0
meanradius_ctrl_ens[:,1]=0.0
meanradius_ctrl_ens[:,2]=cloud_stat_ctrlD2[:,11]
meanradius_ctrl_ens[:,3]=cloud_stat_ctrlD3[:,11]
meanradius_ctrl_ens[:,4]=cloud_stat_ctrlD4[:,11]
meanradius_ctrl_ens[:,5]=cloud_stat_ctrlD5[:,11]
meanradius_ctrl_ens[:,6]=cloud_stat_ctrlD6[:,11]
meanradius_ctrl_ens[:,7]=cloud_stat_ctrlD7[:,11]
meanradius_ctrl_ens[:,8]=cloud_stat_ctrlD8[:,11]
meanradius_ctrl_ens[:,9]=cloud_stat_ctrlD9[:,11]
meanradius_ctrl_ens[:,10]=cloud_stat_ctrlD10[:,11]

meanarea_ctrl_ens=N.zeros((nhr,11))
meanarea_ctrl_ens[:,0]=0.0
meanarea_ctrl_ens[:,1]=0.0
meanarea_ctrl_ens[:,2]=cloud_stat_ctrlD2[:,2]
meanarea_ctrl_ens[:,3]=cloud_stat_ctrlD3[:,2]
meanarea_ctrl_ens[:,4]=cloud_stat_ctrlD4[:,2]
meanarea_ctrl_ens[:,5]=cloud_stat_ctrlD5[:,2]
meanarea_ctrl_ens[:,6]=cloud_stat_ctrlD6[:,2]
meanarea_ctrl_ens[:,7]=cloud_stat_ctrlD7[:,2]
meanarea_ctrl_ens[:,8]=cloud_stat_ctrlD8[:,2]
meanarea_ctrl_ens[:,9]=cloud_stat_ctrlD9[:,2]
meanarea_ctrl_ens[:,10]=cloud_stat_ctrlD10[:,2]

medianarea_ctrl_ens=N.zeros((nhr,11))
medianarea_ctrl_ens[:,0]=0.0
medianarea_ctrl_ens[:,1]=0.0
medianarea_ctrl_ens[:,2]=cloud_stat_ctrlD2[:,7]
medianarea_ctrl_ens[:,3]=cloud_stat_ctrlD3[:,7]
medianarea_ctrl_ens[:,4]=cloud_stat_ctrlD4[:,7]
medianarea_ctrl_ens[:,5]=cloud_stat_ctrlD5[:,7]
medianarea_ctrl_ens[:,6]=cloud_stat_ctrlD6[:,7]
medianarea_ctrl_ens[:,7]=cloud_stat_ctrlD7[:,7]
medianarea_ctrl_ens[:,8]=cloud_stat_ctrlD8[:,7]
medianarea_ctrl_ens[:,9]=cloud_stat_ctrlD9[:,7]
medianarea_ctrl_ens[:,10]=cloud_stat_ctrlD10[:,7]

num_ctrl_ens=N.zeros((nhr,11))
num_ctrl_ens[:,0]=0.0
num_ctrl_ens[:,1]=0.0
num_ctrl_ens[:,2]=cloud_stat_ctrlD2[:,1]
num_ctrl_ens[:,3]=cloud_stat_ctrlD3[:,1]
num_ctrl_ens[:,4]=cloud_stat_ctrlD4[:,1]
num_ctrl_ens[:,5]=cloud_stat_ctrlD5[:,1]
num_ctrl_ens[:,6]=cloud_stat_ctrlD6[:,1]
num_ctrl_ens[:,7]=cloud_stat_ctrlD7[:,1]
num_ctrl_ens[:,8]=cloud_stat_ctrlD8[:,1]
num_ctrl_ens[:,9]=cloud_stat_ctrlD9[:,1]
num_ctrl_ens[:,10]=cloud_stat_ctrlD10[:,1]


stdvR_ctrl_ens=N.zeros((nhr,11))
stdvR_ctrl_ens[:,0]=0.0
stdvR_ctrl_ens[:,1]=0.0
stdvR_ctrl_ens[:,2]=cloud_stat_ctrlD2[:,3]
stdvR_ctrl_ens[:,3]=cloud_stat_ctrlD3[:,3]
stdvR_ctrl_ens[:,4]=cloud_stat_ctrlD4[:,3]
stdvR_ctrl_ens[:,5]=cloud_stat_ctrlD5[:,3]
stdvR_ctrl_ens[:,6]=cloud_stat_ctrlD6[:,3]
stdvR_ctrl_ens[:,7]=cloud_stat_ctrlD7[:,3]
stdvR_ctrl_ens[:,8]=cloud_stat_ctrlD8[:,3]
stdvR_ctrl_ens[:,9]=cloud_stat_ctrlD9[:,3]
stdvR_ctrl_ens[:,10]=cloud_stat_ctrlD10[:,3]

num_ctrl_ens_mean=N.zeros((nhr))
num_ctrl_ens_meanP=N.zeros((nhr))
num_ctrl_ens_meanM=N.zeros((nhr))
meanradius_ctrl_ens_mean=N.zeros((nhr))
medianarea_ctrl_ens_mean=N.zeros((nhr))
meanradius_ctrl_ens_meanP=N.zeros((nhr))
meanradius_ctrl_ens_meanM=N.zeros((nhr))
meanarea_ctrl_ens_mean=N.zeros((nhr))
meanarea_ctrl_ens_meanP=N.zeros((nhr))
meanarea_ctrl_ens_meanM=N.zeros((nhr))
stdvR_ctrl_ens_mean=N.zeros((nhr))
stdvR_ctrl_ens_meanP=N.zeros((nhr))
stdvR_ctrl_ens_meanM=N.zeros((nhr))


for j in N.arange(nhr):
     num_ctrl_ens_mean[j]=num_ctrl_ens[j,2:11].mean()
     num_ctrl_ens_meanM[j]=num_ctrl_ens_mean[j]-np.std(num_ctrl_ens[j,2:11])
     num_ctrl_ens_meanP[j]=num_ctrl_ens_mean[j]+np.std(num_ctrl_ens[j,2:11])
     meanradius_ctrl_ens_mean[j]=meanradius_ctrl_ens[j,2:11].mean()
     meanradius_ctrl_ens_meanM[j]=meanradius_ctrl_ens_mean[j]-np.std(meanradius_ctrl_ens[j,2:11])
     meanradius_ctrl_ens_meanP[j]=meanradius_ctrl_ens_mean[j]+np.std(meanradius_ctrl_ens[j,2:11])

     meanarea_ctrl_ens_mean[j]=meanarea_ctrl_ens[j,2:11].mean()
     meanarea_ctrl_ens_meanM[j]=meanarea_ctrl_ens_mean[j]-np.std(meanarea_ctrl_ens[j,2:11])
     meanarea_ctrl_ens_meanP[j]=meanarea_ctrl_ens_mean[j]+np.std(meanarea_ctrl_ens[j,2:11])
     medianarea_ctrl_ens_mean[j]=medianarea_ctrl_ens[j,2:11].mean()

     stdvR_ctrl_ens_mean[j]=stdvR_ctrl_ens[j,2:11].mean()
     stdvR_ctrl_ens_meanM[j]=stdvR_ctrl_ens_mean[j]-np.std(stdvR_ctrl_ens[j,2:11])
     stdvR_ctrl_ens_meanP[j]=stdvR_ctrl_ens_mean[j]+np.std(stdvR_ctrl_ens[j,2:11])


cloud_stat_rthqvD10=np.load('cloud_stat_rthqvD10.npy')
cloud_stat_rthqvD9=np.load('cloud_stat_rthqvD9.npy')
cloud_stat_rthqvD8=np.load('cloud_stat_rthqvD8.npy')
cloud_stat_rthqvD7=np.load('cloud_stat_rthqvD7.npy')
cloud_stat_rthqvD6=np.load('cloud_stat_rthqvD6.npy')
cloud_stat_rthqvD5=np.load('cloud_stat_rthqvD5.npy')
cloud_stat_rthqvD4=np.load('cloud_stat_rthqvD4.npy')
cloud_stat_rthqvD3=np.load('cloud_stat_rthqvD3.npy')
cloud_stat_rthqvD2=np.load('cloud_stat_rthqvD2.npy')

meanradius_rthqv_ens=N.zeros((nhr,11))
meanradius_rthqv_ens[:,0]=0.0
meanradius_rthqv_ens[:,1]=0.0
meanradius_rthqv_ens[:,2]=cloud_stat_rthqvD2[:,11]
meanradius_rthqv_ens[:,3]=cloud_stat_rthqvD3[:,11]
meanradius_rthqv_ens[:,4]=cloud_stat_rthqvD4[:,11]
meanradius_rthqv_ens[:,5]=cloud_stat_rthqvD5[:,11]
meanradius_rthqv_ens[:,6]=cloud_stat_rthqvD6[:,11]
meanradius_rthqv_ens[:,7]=cloud_stat_rthqvD7[:,11]
meanradius_rthqv_ens[:,8]=cloud_stat_rthqvD8[:,11]
meanradius_rthqv_ens[:,9]=cloud_stat_rthqvD9[:,11]
meanradius_rthqv_ens[:,10]=cloud_stat_rthqvD10[:,11]

meanarea_rthqv_ens=N.zeros((nhr,11))
meanarea_rthqv_ens[:,0]=0.0
meanarea_rthqv_ens[:,1]=0.0
meanarea_rthqv_ens[:,2]=cloud_stat_rthqvD2[:,2]
meanarea_rthqv_ens[:,3]=cloud_stat_rthqvD3[:,2]
meanarea_rthqv_ens[:,4]=cloud_stat_rthqvD4[:,2]
meanarea_rthqv_ens[:,5]=cloud_stat_rthqvD5[:,2]
meanarea_rthqv_ens[:,6]=cloud_stat_rthqvD6[:,2]
meanarea_rthqv_ens[:,7]=cloud_stat_rthqvD7[:,2]
meanarea_rthqv_ens[:,8]=cloud_stat_rthqvD8[:,2]
meanarea_rthqv_ens[:,9]=cloud_stat_rthqvD9[:,2]
meanarea_rthqv_ens[:,10]=cloud_stat_rthqvD10[:,2]

medianarea_rthqv_ens=N.zeros((nhr,11))
medianarea_rthqv_ens[:,0]=0.0
medianarea_rthqv_ens[:,1]=0.0
medianarea_rthqv_ens[:,2]=cloud_stat_rthqvD2[:,7]
medianarea_rthqv_ens[:,3]=cloud_stat_rthqvD3[:,7]
medianarea_rthqv_ens[:,4]=cloud_stat_rthqvD4[:,7]
medianarea_rthqv_ens[:,5]=cloud_stat_rthqvD5[:,7]
medianarea_rthqv_ens[:,6]=cloud_stat_rthqvD6[:,7]
medianarea_rthqv_ens[:,7]=cloud_stat_rthqvD7[:,7]
medianarea_rthqv_ens[:,8]=cloud_stat_rthqvD8[:,7]
medianarea_rthqv_ens[:,9]=cloud_stat_rthqvD9[:,7]
medianarea_rthqv_ens[:,10]=cloud_stat_rthqvD10[:,7]


num_rthqv_ens=N.zeros((nhr,11))
num_rthqv_ens[:,0]=0.0
num_rthqv_ens[:,1]=0.0
num_rthqv_ens[:,2]=cloud_stat_rthqvD2[:,1]
num_rthqv_ens[:,3]=cloud_stat_rthqvD3[:,1]
num_rthqv_ens[:,4]=cloud_stat_rthqvD4[:,1]
num_rthqv_ens[:,5]=cloud_stat_rthqvD5[:,1]
num_rthqv_ens[:,6]=cloud_stat_rthqvD6[:,1]
num_rthqv_ens[:,7]=cloud_stat_rthqvD7[:,1]
num_rthqv_ens[:,8]=cloud_stat_rthqvD8[:,1]
num_rthqv_ens[:,9]=cloud_stat_rthqvD9[:,1]
num_rthqv_ens[:,10]=cloud_stat_rthqvD10[:,1]

stdvR_rthqv_ens=N.zeros((nhr,11))
stdvR_rthqv_ens[:,0]=0.0
stdvR_rthqv_ens[:,1]=0.0
stdvR_rthqv_ens[:,2]=cloud_stat_rthqvD2[:,3]
stdvR_rthqv_ens[:,3]=cloud_stat_rthqvD3[:,3]
stdvR_rthqv_ens[:,4]=cloud_stat_rthqvD4[:,3]
stdvR_rthqv_ens[:,5]=cloud_stat_rthqvD5[:,3]
stdvR_rthqv_ens[:,6]=cloud_stat_rthqvD6[:,3]
stdvR_rthqv_ens[:,7]=cloud_stat_rthqvD7[:,3]
stdvR_rthqv_ens[:,8]=cloud_stat_rthqvD8[:,3]
stdvR_rthqv_ens[:,9]=cloud_stat_rthqvD9[:,3]
stdvR_rthqv_ens[:,10]=cloud_stat_rthqvD10[:,3]

num_rthqv_ens_mean=N.zeros((nhr))
num_rthqv_ens_meanP=N.zeros((nhr))
num_rthqv_ens_meanM=N.zeros((nhr))
meanradius_rthqv_ens_mean=N.zeros((nhr))
medianarea_rthqv_ens_mean=N.zeros((nhr))
meanradius_rthqv_ens_meanP=N.zeros((nhr))
meanradius_rthqv_ens_meanM=N.zeros((nhr))
meanarea_rthqv_ens_mean=N.zeros((nhr))
meanarea_rthqv_ens_meanP=N.zeros((nhr))
meanarea_rthqv_ens_meanM=N.zeros((nhr))
stdvR_rthqv_ens_mean=N.zeros((nhr))
stdvR_rthqv_ens_meanP=N.zeros((nhr))
stdvR_rthqv_ens_meanM=N.zeros((nhr))

for j in N.arange(nhr):
     num_rthqv_ens_mean[j]=num_rthqv_ens[j,2:11].mean()
     num_rthqv_ens_meanM[j]=num_rthqv_ens_mean[j]-np.std(num_rthqv_ens[j,2:11])
     num_rthqv_ens_meanP[j]=num_rthqv_ens_mean[j]+np.std(num_rthqv_ens[j,2:11])
     meanradius_rthqv_ens_mean[j]=meanradius_rthqv_ens[j,2:11].mean()
     meanradius_rthqv_ens_meanM[j]=meanradius_rthqv_ens_mean[j]-np.std(meanradius_rthqv_ens[j,2:11])
     meanradius_rthqv_ens_meanP[j]=meanradius_rthqv_ens_mean[j]+np.std(meanradius_rthqv_ens[j,2:11])

     meanarea_rthqv_ens_mean[j]=meanarea_rthqv_ens[j,2:11].mean()
     meanarea_rthqv_ens_meanM[j]=meanarea_rthqv_ens_mean[j]-np.std(meanarea_rthqv_ens[j,2:11])
     meanarea_rthqv_ens_meanP[j]=meanarea_rthqv_ens_mean[j]+np.std(meanarea_rthqv_ens[j,2:11])
     medianarea_rthqv_ens_mean[j]=medianarea_rthqv_ens[j,2:11].mean()

     stdvR_rthqv_ens_mean[j]=stdvR_rthqv_ens[j,2:11].mean()
     stdvR_rthqv_ens_meanM[j]=stdvR_rthqv_ens_mean[j]-np.std(stdvR_rthqv_ens[j,2:11])#/3.0
     stdvR_rthqv_ens_meanP[j]=stdvR_rthqv_ens_mean[j]+np.std(stdvR_rthqv_ens[j,2:11])#/3.0



GX=N.zeros((2))
GY=N.zeros((2))
GY[0]=0.0
GY[1]=900
GX[0]=2.5
GX[1]=2.5
MX=N.zeros((2))
MY=N.zeros((2))
MY[0]=0.0
MY[1]=900
MX[0]=6.5
MX[1]=6.5


fig = plt.figure(figsize=(16,20))
ax = fig.add_axes([0.1,0.1,0.8,0.8])

plt.subplot(411)
plt.plot(area_hist[1:nhist], area_freq_ctrl_ens_gs1_mean[1:nhist], linestyle=':',color="green",lw=2, label='[0.25-0.75] h')
plt.fill_between(area_hist[1:nhist], area_freq_ctrl_ens_gs1_meanM[1:nhist], area_freq_ctrl_ens_gs1_meanP[1:nhist], alpha=0.5, edgecolor='lightgreen', facecolor='lightgreen')
plt.plot(area_hist[1:nhist], area_freq_ctrl_ens_gs2_mean[1:nhist], linestyle='--', color="green",lw=2, label='[1-1.25] h')
plt.fill_between(area_hist[1:nhist], area_freq_ctrl_ens_gs2_meanM[1:nhist], area_freq_ctrl_ens_gs2_meanP[1:nhist], alpha=0.5, edgecolor='lightgreen', facecolor='lightgreen')
plt.plot(area_hist[1:nhist], area_freq_ctrl_ens_gs3_mean[1:nhist], linestyle='-.', color="green",lw=2, label='[1.5-1.75] h')
plt.fill_between(area_hist[1:nhist], area_freq_ctrl_ens_gs3_meanM[1:nhist], area_freq_ctrl_ens_gs3_meanP[1:nhist], alpha=0.5, edgecolor='lightgreen', facecolor='lightgreen')

plt.plot(area_hist[1:nhist], area_freq_ctrl_ens_gs4_mean[1:nhist], color="green",lw=2, label='[2-2.5] h')
plt.fill_between(area_hist[1:nhist], area_freq_ctrl_ens_gs4_meanM[1:nhist], area_freq_ctrl_ens_gs4_meanP[1:nhist], alpha=0.5, edgecolor='lightgreen', facecolor='lightgreen')

plt.plot(area_hist[1:nhist], area_freq_ctrl_ens_ms_mean[1:nhist] , color="black", lw=3,label='[2.75-6.5] h')
plt.fill_between(area_hist[1:nhist], area_freq_ctrl_ens_ms_meanM[1:nhist], area_freq_ctrl_ens_ms_meanP[1:nhist], alpha=0.5, edgecolor='grey', facecolor='grey')
plt.plot(area_hist[1:nhist], area_freq_ctrl_ens_ss_mean[1:nhist] , lw=3, color="blue", label='[6.75-9.5] h')
plt.fill_between(area_hist[1:nhist], area_freq_ctrl_ens_ss_meanM[1:nhist], area_freq_ctrl_ens_ss_meanP[1:nhist], alpha=0.5, edgecolor='cornflowerblue', facecolor='cornflowerblue')


plt.legend()

plt.xscale('log')

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(1.15, .55, r'a)', fontdict=font, color="black", fontsize='x-large')

plt.ylabel('PDF ', fontsize='x-large')
plt.xlabel('A$_i$ (km$^2$) ', fontsize='x-large')
plt.xlim(0,50)
plt.ylim(0,0.6)
axes = plt.gca()
axes.set_xlim([0, 50])
axes.set_xticks([1,4, 10,50])
axes.set_xticklabels([1,4,10,50])


plt.subplot(412)
plt.plot(hour_ctrl_mean_shit0[3:31], meanradius_ctrl_ens_mean[3:31],lw=2, color="black", label='C')
plt.fill_between(hour_ctrl_mean_shit0[3:31], meanradius_ctrl_ens_meanM[3:31], meanradius_ctrl_ens_meanP[3:31], alpha=0.5, edgecolor='grey', facecolor='grey')

plt.plot(hour_mhalf_mean_shit0[5:31], meanradius_mhalf_ens_mean[5:31],lw=2, color="blue", label='W')
plt.fill_between(hour_mhalf_mean_shit0[5:31], meanradius_mhalf_ens_meanM[5:31], meanradius_mhalf_ens_meanP[5:31], alpha=0.5, edgecolor='cornflowerblue', facecolor='cornflowerblue')

plt.plot(hour_phalf_mean_shit0[1:31], meanradius_phalf_ens_mean[1:31],lw=2, color="red", label='S')
plt.fill_between(hour_phalf_mean_shit0[1:31], meanradius_phalf_ens_meanM[1:31], meanradius_phalf_ens_meanP[1:31], alpha=0.5, edgecolor='lightcoral', facecolor='lightcoral')
plt.legend(bbox_to_anchor=(0.99, 0.1), loc=4, borderaxespad=0., fontsize='large')
#plt.fill_between(hour_ctrl_mean[3:12], meanradius_ctrl_ens_meanM[3:12], meanradius_ctrl_ens_meanP[3:12], alpha=0.5, edgecolor='lime', facecolor='lime')
#plt.fill_between(hour_ctrl_mean[25:31], meanradius_ctrl_ens_meanM[25:31], meanradius_ctrl_ens_meanP[25:31], alpha=0.5, edgecolor='cornflowerblue', facecolor='cornflowerblue')

#plt.plot(GX[:], GY[:],lw=2, color="black")
#plt.plot(MX[:], MY[:],lw=2, color="black")


plt.ylabel('R (km) ', fontsize='x-large')
plt.xlabel('Time after triggering (hour)', fontsize='x-large')
plt.xlim(0,10)
plt.xlim(0,10)
axes = plt.gca()
axes.set_xlim([0, 10])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10])

plt.ylim(0,2.5)
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 2.2, r'b)', fontdict=font, color="black", fontsize='x-large')

plt.subplot(413)

plt.plot(hour_ctrl_mean_shit0[3:31], stdvR_ctrl_ens_mean[3:31],lw=2, color="black", label='Control')
plt.fill_between(hour_ctrl_mean_shit0[3:31], stdvR_ctrl_ens_meanM[3:31], stdvR_ctrl_ens_meanP[3:31], alpha=0.5, edgecolor='grey', facecolor='grey')


plt.plot(hour_mhalf_mean_shit0[5:31], stdvR_mhalf_ens_mean[5:31],lw=2, color="blue", label='Control')
plt.fill_between(hour_mhalf_mean_shit0[5:31], stdvR_mhalf_ens_meanM[5:31], stdvR_mhalf_ens_meanP[5:31], alpha=0.5, edgecolor='cornflowerblue', facecolor='cornflowerblue')

plt.plot(hour_phalf_mean_shit0[1:31], stdvR_phalf_ens_mean[1:31],lw=2, color="red", label='Control')
plt.fill_between(hour_phalf_mean_shit0[1:31], stdvR_phalf_ens_meanM[1:31], stdvR_phalf_ens_meanP[1:31], alpha=0.5, edgecolor='lightcoral', facecolor='lightcoral')

plt.ylabel('$\sigma_R$ (km)', fontsize='x-large')
plt.xlabel('Time after triggering (hour)', fontsize='x-large')
plt.xlim(0,10)
plt.xlim(0,10)
axes = plt.gca()
axes.set_xlim([0, 10])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10])

plt.ylim(0,6)
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5,5, r'c)', fontdict=font, color="black", fontsize='x-large')






plt.subplot(414)
plt.plot(hour_ctrl_mean_shit0[3:31], num_ctrl_ens_mean[3:31],lw=2, color="black", label='C')
plt.fill_between(hour_ctrl_mean_shit0[3:31], num_ctrl_ens_meanM[3:31], num_ctrl_ens_meanP[3:31], alpha=0.5, edgecolor='grey', facecolor='grey')


plt.plot(hour_mhalf_mean_shit0[5:31], num_mhalf_ens_mean[5:31],lw=2, color="blue", label='S')
plt.fill_between(hour_mhalf_mean_shit0[5:31], num_mhalf_ens_meanM[5:31], num_mhalf_ens_meanP[5:31], alpha=0.5, edgecolor='cornflowerblue', facecolor='cornflowerblue')

plt.plot(hour_phalf_mean_shit0[1:31], num_phalf_ens_mean[1:31],lw=2, color="red", label='W')
plt.fill_between(hour_phalf_mean_shit0[1:31], num_phalf_ens_meanM[1:31], num_phalf_ens_meanP[1:31], alpha=0.5, edgecolor='lightcoral', facecolor='lightcoral')

plt.ylabel('N ', fontsize='x-large')
plt.xlabel('Time after triggering (hour)', fontsize='x-large')
plt.xlim(0,10)
plt.xlim(0,10)
axes = plt.gca()
axes.set_xlim([0, 10])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10])

plt.ylim(0,900)
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5,820, r'd)', fontdict=font, color="black", fontsize='x-large')

plt.savefig('Figure4.png')


fig = plt.figure(figsize=(16,20))
ax = fig.add_axes([0.1,0.1,0.8,0.8])

plt.subplot(411)
plt.plot(hour_ctrl_mean_shit0[3:31], meanradius_ctrl_ens_mean[3:31],lw=2, color="black", label='C')
plt.fill_between(hour_ctrl_mean_shit0[3:31], meanradius_ctrl_ens_meanM[3:31], meanradius_ctrl_ens_meanP[3:31], alpha=0.5, edgecolor='grey', facecolor='grey')

plt.plot(hour_rthqv_mean_shit0[4:31], meanradius_rthqv_ens_mean[4:31],lw=2, color="red", label='H')
plt.fill_between(hour_rthqv_mean_shit0[4:31], meanradius_rthqv_ens_meanM[4:31], meanradius_rthqv_ens_meanP[4:31], alpha=0.5, edgecolor='lightcoral', facecolor='lightcoral')
plt.legend()
plt.ylabel('R (km) ', fontsize='x-large')
plt.xlabel('Time after triggering (hour)', fontsize='x-large')
plt.xlim(0,10)
plt.xlim(0,10)
axes = plt.gca()
axes.set_xlim([0, 10])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10])

plt.ylim(0,2.5)
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5, 2.2, r'a)', fontdict=font, color="black", fontsize='x-large')


plt.subplot(412)

plt.plot(hour_ctrl_mean_shit0[3:31], stdvR_ctrl_ens_mean[3:31],lw=2, color="black", label='Control')
plt.fill_between(hour_ctrl_mean_shit0[3:31], stdvR_ctrl_ens_meanM[3:31], stdvR_ctrl_ens_meanP[3:31], alpha=0.5, edgecolor='grey', facecolor='grey')

plt.plot(hour_rthqv_mean_shit0[4:31], stdvR_rthqv_ens_mean[4:31],lw=2, color="red", label='Control')
plt.fill_between(hour_rthqv_mean_shit0[4:31], stdvR_rthqv_ens_meanM[4:31], stdvR_rthqv_ens_meanP[4:31], alpha=0.5, edgecolor='lightcoral', facecolor='lightcoral')

plt.ylabel('$\sigma_R$ (km)', fontsize='x-large')
plt.xlabel('Time after triggering (hour)', fontsize='x-large')
plt.xlim(0,10)
plt.xlim(0,10)
axes = plt.gca()
axes.set_xlim([0, 10])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10])

plt.ylim(0,4)
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5,3.5, r'b)', fontdict=font, color="black", fontsize='x-large')



plt.subplot(413)

plt.plot(hour_ctrl_mean_shit0[3:31], num_ctrl_ens_mean[3:31],lw=2, color="black", label='Control')
plt.fill_between(hour_ctrl_mean_shit0[3:31], num_ctrl_ens_meanM[3:31], num_ctrl_ens_meanP[3:31], alpha=0.5, edgecolor='grey', facecolor='grey')

plt.plot(hour_rthqv_mean_shit0[4:31], num_rthqv_ens_mean[4:31],lw=2, color="red", label='Control')
plt.fill_between(hour_rthqv_mean_shit0[4:31], num_rthqv_ens_meanM[4:31], num_rthqv_ens_meanP[4:31], alpha=0.5, edgecolor='lightcoral', facecolor='lightcoral')

plt.ylabel('N ', fontsize='x-large')
plt.xlabel('Time after triggering (hour)', fontsize='x-large')
plt.xlim(0,10)
plt.xlim(0,10)
axes = plt.gca()
axes.set_xlim([0, 10])
axes.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
axes.set_xticklabels([0,1,2,3,4,5,6,7,8,9,10])

plt.ylim(0,1200)
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
plt.text(0.5,1100, r'c)', fontdict=font, color="black", fontsize='x-large')

plt.savefig('Figure9.png')

