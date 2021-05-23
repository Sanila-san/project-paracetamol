# project-paracetamol

## What is it
This code was originally created to calculate the measurement errors caused by the measurement error of the linear geometric dimensions of the measures on the substrate, and the residual errors of the vector network analyzer caused by this error. These calculations are based on the mathematical apparatus described in the following sources:

 1. Heinrich, W. (1993). Quasi-TEM description of MMIC coplanar lines including conductor-loss effects. _IEEE transactions on microwave theory and techniques_, _41_(1), 45-52.
 2. Recommendation MI 3411-2013. State system for ensuring the uniformity of measurements. Vector network analyzers. Method for determining the metrological characteristics

This is not a production script. This is scientific research tool used for solving special tasks in the context of research on radiotechnical metrology. The main keywords are: coplanar waveguide, devices on the dielectric substrate, vector network analyzer, TRL calibration.

## How it works

This code completely emulates measurement procedure using ideal VNA and ideal DUT using non-ideal calibration substrate. The imperfection of a calibration substrate lies in the fact that there are ideal standards on the substrate, the dimensions of which are measured with an error. It allows us to estimate the influence of coplanar waveguide characterization errors on the overall precision of S-parameters measurements.


The main function is **calerrsim.m**, which should be called from **batch.m**. If no batch calculation is needed, **cal_error_simulation.m** should be used instead of the previous two. Other .m-files are necessary for all ways to run the simulation.

Brief workflow:

 1. Initialization and setting necessary parameters
 2. Setting randomized parameters of calibration substrate
 3. Calculating error vectors of VNA
 4. Measuring ideal DUT with error
 5. Saving calculated data into TSV-file for further analysis

Materials used in the calibration substrate description are gold for conductors and alumina for the substrate. Measurement precision for geometric measurements is typical for Olympus LEXT OLS 5000 microscope and pretty close to plus-minus 1 um.

> Written with [StackEdit](https://stackedit.io/).
