# Climate risk profiles

Base scripts.
-*indices.R*: contains basic functions necessary to calculate agro-climatic indices
-*win_parallelization.R*: functions to parallelize in Windows SO

Order to follow.
-*00_obtain_soil_capacity.R*: calculate soil capacity as an input for posterior indices calculation
-*01_gcm_resampling_windows.R*: process GCM data, clipping to study area, and apply resampling
-*02_gcm_bias_correction_qmap.R*: develop quantile mapping as a technique for GCM bias correction
-*03_loading_obs_climate.R*: prepare historic climate observational data
-*04_loading_gcm_climate.R*: prepare GCM climate data
-*05_calc_indices.R*: calculate partial agro-climatic indices, which do not require solar radiation data
-*06_calc_indices_srad.R*: calculate agro-climatic indices, which requires of soil and solar radiation data
-*07_graphs_do_climatology.R*: do a climatology graph using observational climate data
-*08_graphs_do_time_series.R*: do time series per index/season using historic and future
-*09_graphs_do_maps.R*: do maps for historic, future, and differences