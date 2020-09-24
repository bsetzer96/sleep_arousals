# sleep_arousals

main script: arousal_average_mghthalamus.m
+ extracts timeseries of ROIs around arousals
+ loops through subjects and runs 
  - extractArous_mghthalamus.m extracts arousal times
  - extractArousalTS.m extracts window of ROIs around arousal
+ makes plots of subject averages and ROI total averages
+ Does a cross-correlation analysis of averaged timeseries
  - in crosscoranal.m
  - outputs heatmaps of lags, correlations, and a bar plot of average correlation
+ Does a cross-correlation analysis of individual arousal timeseries and then takes average lag
  - singleArousalLags.m
  - outputs heatmaps of lags & correlations, bar plots of average lags, and bar plot of individual lags
+ Does a bootstrap analysis where sample timeseries are averaged and then cross-correlated.
  - bootstrapanal.m
  - following code makes heatmaps and bar plots of lag and correlation results
