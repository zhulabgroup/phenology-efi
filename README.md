# EFI

v1.3 changes
* use same basis vectors
* change priors
* increase basis number

v1.2 changes
* one model for all time of the year with temporal correlation, instead of one model for each month
* forecast with changing basis vectors that are in +- 30 day of the year
* error propagation turned off
* remove gcc data with sd > 0.0001
* better workflow to update data by checking new data and appending to old data
* code available for out-of-sample test

v1.1 changes
* remove daily lags of weather data
* change size of moving window from 16 to 8 days for weather data
* change number of moving window from 4 to 8 for weather data