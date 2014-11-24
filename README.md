IceChrono
=========

A statistical and physical model to optimize chronologies of deep polar ice cores.


Where can I get help on IceChrono?
----------------------------------

A mailing list has been set up on Google Groups:
https://groups.google.com/forum/?hl=en#!forum/icechrono


How to download IceChrono?
--------------------------

Go here:
https://github.com/parrenin/IceChrono/releases
and choose the release you want to download (usually the latest one).


What do I need to run IceChrono?
--------------------------------

IceChrono is a python2 software, using the following modules:
- sys
- os
- time
- math
- numpy
- matplotlib
- multiprocessing
- warnings
- scipy

If you installed a distribution of scipy:
http://scipy.org/install.html
all these modules should be included.

IceChrono has been tested on debian 7 (with an update of matplotlib>=0.11) and on ubuntu 13.10/14.04.


How to run IceChrono?
---------------------

Go into your IceChrono directory and run the following command:

`python IceChrono.py exp_directory`

where `exp_directory` is the name of the directory where you have set up your experiment (it contains the observations, the background scenarios and all the other parameters). For example, you can test IceChrono using the `AICC2012-VLR` experiment provided for your convenience:

`python IceChrono.py AICC2012-VLR`

It takes about 10 mn to run on a recent computer. It is an AICC2012-like experiment, albeit whith a Very Low Resolution.


What is the structure of an experiment directory?
-------------------------------------------------

You can have a look at the provided `AICC2012-LR` directory.

You have three general files:
- `parameters.py`                                           : contains general parameters for the experiment
- `parameters-CovariancePrior-AllDrillings-init.py`         : defines the covariance matrices of the background
- `parameters-AllDrillings.py`                              : defines drilling parameters that are the same for all drillings (there are overidded by drilling specific parameters).
- `parameters-CovarianceObservations-AllDrillings.py`       : defines the covariance of the observations that are the same for all drillings  (there are overidded by drilling specific parameters).
- `parameters-CovarianceObservations-AllDrillingCouples.py` : defines the covariance for the observations that are the same for all drilling couples  (there are overidded by drilling couple specific parameters).

Then you have one directory per drilling, which contains:
- `parameters.py`                           : all the drilling specific parameters
- `parameters-CovarianceObservations.py`    : this file allows to define the correlation of drilling specific observations
- `density-prior.txt`                       : depth / relative density
- `accu-prior.txt`                          : depth / background accu (in ice-equivalent)
- `LID-prior.txt`                           : depth / background Lock-in-Depth
- `thinning-prior.txt`                      : depth / background thinning function
- `ice_age.txt`                             : depth / age / sigma for ice age observations
- `gas_age.txt`                             : depth / age / sigma for gas age observations
- `Ddepth.txt`                              : depth / Delta-depth / sigma for Delta-depth observations

Then you have one directory per drilling couple, which contains:
- `parameters-CovarianceObservations.py`    : this file allows to define the correlation of drilling couple specific observations
- `ice_depth.txt`           : depth1 / depth2 / sigma on age for ice-ice stratigraphic links observations
- `gas_depth.txt`           : depth1 / depth2 / sigma on age for gas-gas stratigraphic links observations
- `icegas_depth.txt`        : depth1 / depth2 / sigma on age for ice-gas stratigraphic links observations
- `gasice_depth.txt`        : depth1 / depth2 / sigma on age for gas-ice stratigraphic links observations
- `ice_age_intervals.txt`   : depth\_top / depth\_bottom / duration / sigma for dated ice intervals observations
- `gas_age_intervals.txt`   : depth\_top / depth\_bottom / duration / sigma for dated gas intervals observations

You can use whatever units you want but they need to be consistent. For example, if you use meters for the depths and years for the dated horizons, you need to use meters per years for the accumulation rates.

 
What is the structure of the general `parameters.py` file?
--------------------------------------------------------

Have a look at the file `AICC2012-VLR/parameters.py`, it is commented.


What is the structure of a drilling `parameters.py` file?
---------------------------------------------------------

Have a look at the files `AICC2012-VLR/EDC/parameters.py`, it is commented.


How to set up correlation matrices for the observations?
--------------------------------------------------------

You need to know a little bit of python to do that.
Feel free to contact the author if you need assistance.

For drilling specific observations, you set up the correlation matrices in the file `parameters-CovarianceObservations.py` in the drilling directory.
- `self.icemarkers_correlation`     : for ice dated horizons
- `self.gasmarkers_correlation`     : for gas dated horizons
- `self.iceintervals_correlation`   : for ice dated intervals
- `self.gasintervals_correlation`   : for gas dated intervals
- `self.Ddepth_correlation`         : for Delta-depth observations

For drilling couple specific observations (stratigraphic links), you set up the correlation matrices in the file `parameters-CovarianceObservations.py` in the drilling couple directory.
- `self.iceicemarkers_correlation`  : for ice-ice stratigraphic links
- `self.gasgasmarkers_correlation`  : for gas-gas stratigraphic links
- `self.icegasmarkers_correlation`  : for ice-gas stratigraphic links
- `self.gasicemarkers_correlation`  : for gas-ice stratigraphic links


How to set up the `parameters-CovariancePrior-AllDrillings-init.py` file?
-------------------------------------------------------------------------

You need to know a little bit of python to do that.
Have a look at the `AICC2012-VLR` experiment, it is the easiest way to understand how it works.
Feel free to contact the author if you need assistance.

You need to define:
- `self.correlation_corr_a`     : the correlation matix for the accu correction function
- `self.correlation_corr_LID`   : the correlation matix for the LID correction function
- `self.correlation_corr_tau`   : the correlation matix for the thinning correction function
- `self.sigmap_corr_a`          : the standard deviation of the accu correction function
- `self.sigmap_corr_LID`        : the standard deviation of the LID correction function
- `self.sigmap_corr_tau`        : the standard deviation of the thinning correction function


