clc; clear all;
## Call function calerrsim(standardLen, iterations, idealDUTs11, idealDUTs21) 
## with the following parameters:
## - standardLen — length of the coplanar waveguide used as a standard
##    available values are 900, 1800, 3500, and 5250 micrometers. Check standard
##    lengths using standard substrate datasheet
## - iterations — the number of iterations of modeling. The more iterations, the 
##    more time consumed, the more data to analyze.
## - idealDUTs11 — reflection coefficient of the virtual device under test 
##    virtually measured by virtually calibrated virtual vector network analyzer.
## - idealDUTs21 — transmission coefficient of the virtual device under test 
##    virtually measured by virtually calibrated virtual vector network analyzer.
## NOTE: S21 and S11 are set in abs|LinMag|.
##    
## This code runs Monte-Carlo calculation of virtual measurement of the ideal DUT
## using ideal VNA and non-ideal calibration substrate. Calibration procedure is
## TRL using ideal Thru standard. The result of the calculation is the estimation
## of the errors in the measurement result.

## For example, choose one of the presets below ↓

##calerrsim(900, 50, 0.33333, 1)
##calerrsim(900, 50, 0.33333, 0.0003)
##calerrsim(900, 50, 0.005, 1)
##calerrsim(900, 50, 0.005, 0.0003)
##
##calerrsim(1800, 250, 0.33333, 1)
##calerrsim(1800, 50, 0.33333, 0.0003)
##calerrsim(1800, 50, 0.005, 1)
##calerrsim(1800, 50, 0.005, 0.0003)
##
##calerrsim(3500, 50, 0.33333, 1)
##calerrsim(3500, 50, 0.33333, 0.0003)
##calerrsim(3500, 50, 0.005, 1)
##calerrsim(3500, 50, 0.005, 0.0003)
##
##calerrsim(5250, 50, 0.33333, 1)
##calerrsim(5250, 50, 0.33333, 0.0003)
##calerrsim(5250, 50, 0.005, 1)

calerrsim(5250, 1000, 0.005, 1.0000)