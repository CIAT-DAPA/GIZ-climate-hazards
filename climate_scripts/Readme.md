# Climate risk profiles

Spatial resolution: 5 km

### Base scripts

* **indices.R**: basic functions needed for calculating agro-climatic indices
* **win_parallelization.R**: functions to parallelize in Windows SO

### Running order

* **00_obtain_soil_capacity.R**: calculate soil capacity for all countries at pixel level. Just need to be run once
* **01_gcm_data_processing.R**: process GCM data (historic and future) at country level: clip to country extent and then resampling to 5 km. Just need to be run once per country. This script can be run for several countries at the time
* **02_loading_obs_climate.R**: prepare observational data (historical information of: tmax, tmin, prec, and srad) at county level in the required format to calculate agro-climatic indices. This script can be run for several counties at the time
* **03_loading_gcm_climate.R**: prepare GCM data at county level in the required format to calculate agro-climatic indices. This script can be run for several countries at the time
* **04_gcm_bias_correction_qmap.R**: develop bias correction by mean of quantile mapping at county level for GCM data. This script must be run for one county at the time
* **05_calc_indices.R**: calculate partial agro-climatic indices. Those indices do not require solar radiation and soil data
* **06_calc_indices_srad.R**: calculate agro-climatic indices related with solar radiation and soil data
* **07_graphs_do_climatology.R**: do a climatology graph using observational climate data
* **08_graphs_do_time_series.R**: do time series per index/seasons using historic and future information
* **09_graphs_do_maps.R**: do maps per index/seasons for historic, future, and changes