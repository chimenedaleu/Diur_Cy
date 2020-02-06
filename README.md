# Diur_Cy
Short descriptions for the scripts and data contain in github.com/chimenedaleu/Diur_Cy

The scripts names are: 
script_figure1_2_3_7_8.py for figures 1, 2, 3, 7, and 8 
script_figure5_6_10.py for figures 5, 6, and 10, and  
script_figure4_9.py for figures 4 and 9.

The probability of finding rain and the memory function are calculated using the script: memory_function.py

TheThe script analysis_netcdf_data.py is used to analysed the Netcdf files. The data are then saved to be used in the files script_figure1_2_3_7_8.py,  script_figure5_6_10.py, script_figure4_9.py, memory_function.py.

The Netcdf files are very heavy to be loaded here. There are available on directory: 

Data for the control simulation have the subscript "_ctrl" or "_CTRL"
Data for the weakly forced simulations have the subcript "_mhalf" of "_MHALF"
Data for the strongly forced simulations have the subcript "_phalf" of "_PHALF"
Data for the simulation following homogenization have the subcript "_rthqv of "_RTHQV"
    Data for the simulation following homogenization at night of the first forcing cycle have the subscript "_rthqvd2" or "_RTHQVD2" 
    to be compared with those of the control simulation on the second forcing cycle.
    
    
The following data are timeseries (depend on time only) output every 15 minutes:    
         modets_time_10min_CTRL.npy, surf_flx_CTRL.npy, surf_precip_time_CTRL.npy,
         cloud_MF_time_cb5_CTRL.npy, cloud_frac_time_cb5_CTRL.npy, BCu_MF_time_cb5_CTRL.npy,
         BCu_frac_time_cb5_CTRL.npy, precip_CTRL_f, ACuzt_frac_time_CTRL.npy, BCuzt_frac_time_CTRL.npy,
         ACuzt_MF_time_CTRL.npy, BCuzt_MF_time_CTRL.npy, time_10min_CTRL_f_ens_mean.npy,
         modets_time_10min_phalf.npy, surf_flx_phalf.npy, surf_precip_time_phalf.npy,cloud_MF_time_cb5_phalf.npy, 
         time_10min_phalf_ens_mean.npy,
         modets_time_10min_mhalf.npy, surf_flx_mhalf.npy, surf_precip_time_mhalf.npy,cloud_MF_time_cb5_mhalf.npy, 
         time_10min_mhalf_ens_mean.npy,
         modets_time_10min_RTHQVD2.npy, surf_precip_time_RTHQVD2.npy, cloud_MF_time_cb5_RTHQVD2.npy,
         modets_time_10min_RTHQVD3.npy, surf_precip_time_RTHQVD3.npy, cloud_MF_time_cb5_RTHQVD3.npy, 
         modets_time_10min_RTHQVD4.npy, surf_precip_time_RTHQVD4.npy, cloud_MF_time_cb5_RTHQVD4.npy, 
         modets_time_10min_RTHQVD5.npy, surf_precip_time_RTHQVD5.npy, cloud_MF_time_cb5_RTHQVD5.npy, 
         modets_time_10min_RTHQVD6.npy, surf_precip_time_RTHQVD6.npy, cloud_MF_time_cb5_RTHQVD6.npy,
         modets_time_10min_RTHQVD7.npy, surf_precip_time_RTHQVD7.npy, cloud_MF_time_cb5_RTHQVD7.npy, 
         modets_time_10min_RTHQVD8.npy, surf_precip_time_RTHQVD8.npy, cloud_MF_time_cb5_RTHQVD8.npy, 
         modets_time_10min_RTHQVD9.npy, surf_precip_time_RTHQVD9.npy, cloud_MF_time_cb5_RTHQVD9.npy, 
         modets_time_10min_RTHQVD10.npy, surf_precip_time_RTHQVD10.npy, cloud_MF_time_cb5_RTHQVD10.npy,
         time_10min_RTHQV_ens_mean.npy
         
  The following data depend on nx and ny (number of grid points in he x and y directions): horizontal slices at 24h
         thxy24h_1km_CTRL_snap.npy, thxy24h_3km_CTRL_snap.npy, qvxy24h_1km_CTRL_snap.npy, qvxy24h_3km_CTRL_snap.npy
   
  The following data depend on time, nx and ny (2D surface precipitation fields): surf_precip_xy_time_CTRL.npy,
         surf_precip_xy_time_mhalf.npy, surf_precip_xy_time_phalf.npy, surf_precip_xy_time_RTHQVD2.npy, 
         surf_precip_xy_time_RTHQVD3.npy, surf_precip_xy_time_RTHQVD4.npy, surf_precip_xy_time_RTHQVD5.npy, 
         surf_precip_xy_time_RTHQVD6.npy, surf_precip_xy_time_RTHQVD7.npy, surf_precip_xy_time_RTHQVD8.npy,
        surf_precip_xy_time_RTHQVD9.npy, surf_precip_xy_time_RTHQVD10.npy
        
  The following data depend on nz and ny (number of vertical levels and number of grid points in the y direction): 
  vertical slices at 24h    
         thzy24h_snap.npy, qvzy24h_snap.npy
        
   Values of P[R(A,t0)]   as a function of time
            proR_lag_ens_mean_PHALF.npy, proR_lag_ens_mean_MHALF.npy, proR_lag_ens_mean_RTHQV.npy,  proR_lag_ens_mean_CTRL.npy
        
        
 Values of M[R(A,t0,Dt)] as a function of time lag Dt and time after triggering, t0: (Dt, tO)    
            proRR_ProR2t_lag_dt_RTHQV.npy, proRR_ProR2t_lag_dt_MHALF.npy, proRR_ProR2t_lag_dt_PHALF.npy, 
            proRR_ProR2t_lag_dt_CTRL.npy, proRR_ProR2t_lag_dt_CTRL_xy.npy, proRR_ProR2t_lag_dt_CTRL_22.npy, 
            proRR_ProR2t_lag_dt_CTRL_1010.npy, proRR_ProR2t_lag_dt_CTRL_1515.npy, proRR_ProR2t_lag_dt_CTRL_2525.npy,
            proRR_ProR2t_lag_dt_CTRL_5050.npy
            
 Values of PDFs(nhist, time). The size of each bin=(Max(area)-Min(area))/nhist
              area_freq_ctrlD2.npy, area_freq_ctrlD5.npy, area_freq_ctrlD8.npy,
              area_freq_ctrlD3.npy, area_freq_ctrlD6.npy, area_freq_ctrlD9.npy
              area_freq_ctrlD4.npy, area_freq_ctrlD7.npy, area_freq_ctrlD10.npy
  
#These fields (e.g., cloud_stat_phalfD10) contain the mean statistiques of rainfall events as a function of time
#cloud_stat_phalfD10(time,11) contains the mean radius of rain cores as a function of time for day 10 
#cloud_stat_phalfD10(time,2) contains the mean area of rain cores as a function of time for the 10th diurnal cycle 
#cloud_stat_phalfD10(time,1) contains the number of rain cores as a function of time for the 10th diurnal cycle
#cloud_stat_phalfD10(time,3) contains the standard deviation of rain cores radii as a function of time for the 10th diurnal cycle
         cloud_stat_phalfD2.npy, cloud_stat_phalfD3.npy, cloud_stat_phalfD4.npy, 
          cloud_stat_phalfD5.npy, cloud_stat_phalfD7.npy, cloud_stat_phalfD9.npy, 
          cloud_stat_phalfD6.npy, cloud_stat_phalfD8.npy, cloud_stat_phalfD10.npy, 
          cloud_stat_mhalfD2.npy, cloud_stat_mhalfD3.npy, cloud_stat_mhalfD4.npy, 
          cloud_stat_mhalfD5.npy, cloud_stat_mhalfD7.npy, cloud_stat_mhalfD9.npy, 
          cloud_stat_mhalfD6.npy, cloud_stat_mhalfD8.npy, cloud_stat_mhalfD10.npy, 
          cloud_stat_ctrlD2.npy, cloud_stat_ctrlD5.npy, cloud_stat_ctrlfD8.npy, 
          cloud_stat_ctrlD3.npy, cloud_stat_ctrlD6.npy, cloud_stat_ctrlfD9.npy, 
          cloud_stat_ctrlD4.npy, cloud_stat_ctrlD7.npy, cloud_stat_ctrlfD10.npy, 
          cloud_stat_rthqvD2.npy, cloud_stat_rthqvD5.npy, cloud_stat_rthqvD8.npy, 
          cloud_stat_rthqvD3.npy, cloud_stat_rthqvD6.npy, cloud_stat_rthqvD9.npy, 
          cloud_stat_rthqvD4.npy, cloud_stat_rthqvD7.npy, cloud_stat_rthqvD10.npy, 
