IceChrono
=========

A statistical and physical model to optimize chronologies of deep polar ice cores.


Where can I get help on IceChrono?
----------------------------------

A mailing list has been set up on Google Groups:  
https://groups.google.com/forum/?hl=en#!forum/icechrono  
You just need a google account to access this mailing list.

How to download IceChrono?
--------------------------

Go here:  
https://github.com/parrenin/IceChrono/releases  
and choose the release you want to download (usually the latest one).  
In the downloaded folder, you will find the following files:
- README.md		: is the current documentation of IceChrono.
- LICENCE		: is the IceChrono licence file.
- IceChrono.py		: is the main IceChrono program that you will run.
- IceChronoClasses	: is a set of instructions used by IceChrono.py
- Clean.py		: is a python script to clean all experiment sub-directories
- AICC2012-VLR		: is an example experiment directory for the AICC2012 dating experiment: it contains all the necessary numerical settings, prior information and observations for the different ice cores. It has a Very Low Resolution (VLR) and takes about 5 mn to run an a recent computer.


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
all these modules should be included by default.

IceChrono has been tested on debian 7 (with an update of matplotlib>=0.11), on debian 8, on ubuntu 13.10/14.04 and on Anaconda 2.1.0.


How to run IceChrono?
---------------------

There are two ways of runing IceChrono. Assuming `exp_directory` is the name of your experiment directory, you can: 

1) open an os shell, go into your IceChrono directory and run the following command:

```
python IceChrono.py exp_directory
```

2) You can open a python shell (for example if you are using a distribution of python like Anaconda), go to the IceChrono directory and run the following command:

```
import sys
sys.argv=['IceChrono.py','exp_directory']
execfile('IceChrono.py')
```

The `AICC2012-VLR` experiment directory is provided for you convenience. It is an AICC2012-like experiment, albeit whith a Very Low Resolution. It takes about 5 mn to run on a recent computer.


What are the outputs of a run:
------------------------------

If the run went correctly, it has created output files.

In the main directory, you have the following output file:
- `output.txt`		: only contains the program execution time for now.

In each drilling directory, you have the following output files:
- `output.txt`			: is the main output file. It gives you the posterior estimates and uncertainties of the three input variables (accumulation, LID and thinning) and of the output variables (ice age, gas age, Δdepth, etc.)
- `restart.txt`			: is a restart file, which can be used to restart the optimization experiment where it ended last time.
- `accumulation.pdf`	: is the accumulation figure (with prior and posterior estimates)
- `Ddepth.pdf`			: is the Δdepth figure (with prior estimates, observations and posterior estimates)
- `gas_age.pdf`			: is the gas age figure (with prior estimates, observations and posterior estimates)
- `gaslayerthick.pdf`	: is the gas layer thickness figure (with prior estimates, observations of dated intervals and posterior estimates)
- `ice_age.pdf`			: is the ice age figure (with prior estimates, observations and posterior estimates)
- `icelayerthick.pdf`	: is the ice layer thickness figure (with prior estimates, observations of dated intervals and posterior estimates)
- `LID.pdf`				: is the Lock-In Depth figure (with prior estimates, observations and posterior estimates)
- `thinning.pdf`		: is the thinning figure (with prior and posterior estimates)

In each drilling-couple directory, you have the following output files:
- `gas-gas.pdf`		: is the air-air stratigraphic links figure (with prior and posterior estimates)
- `gas-ice.pdf`		: is the air-ice stratigraphic links figure (with prior and posterior estimates)
- `ice-gas.pdf`		: is the ice-air stratigraphic links figure (with prior and posterior estimates)
- `ice-ice.pdf`		: is the ice-ice stratigraphic links figure (with prior and posterior estimates)

How to clean an experiment directory after a run?
-------------------------------------------------

If your run was successful, it has produced output files and figure files. If your experiment directory is placed under the IceChrono main directory, you can run the following command in an os shell:

```
python Clean.py
```

or in a python shell:

```
execfile('Clean.py')
```


What is the structure of an experiment directory?
-------------------------------------------------

You can have a look at the provided `AICC2012-LR` directory.

You have five general files:
- `parameters.py`                                           : contains general parameters for the experiment
- `parameters-AllDrillings.py`                              : defines drilling parameters that are the same for all drillings (there are overidded by drilling specific parameters).
- `parameters-CovariancePrior-AllDrillings-init.py`         : defines the covariance matrices of the background
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
- `ice_age_intervals.txt`		   			: depth\_top / depth\_bottom / duration / sigma for dated ice intervals observations
- `gas_age_intervals.txt`   		    	: depth\_top / depth\_bottom / duration / sigma for dated gas intervals observations

Then you have one directory per drilling couple, which contains:
- `parameters-CovarianceObservations.py`    : this file allows to define the correlation of drilling couple specific observations
- `ice_depth.txt`           : depth1 / depth2 / sigma on age for ice-ice stratigraphic links observations
- `gas_depth.txt`           : depth1 / depth2 / sigma on age for gas-gas stratigraphic links observations
- `icegas_depth.txt`        : depth1 / depth2 / sigma on age for ice-gas stratigraphic links observations
- `gasice_depth.txt`        : depth1 / depth2 / sigma on age for gas-ice stratigraphic links observations

A few things you need to know to use Icechrono:
1) You can use whatever units you want but they need to be consistent. For example, if you use meters for the depths and years for the dated horizons, you need to use meters per years for the accumulation rates. 
2) The drilling specific parameters override the general parameters for all drillings. In the very same way, the drilling-couple specific parameters override the general parameters for all drilling-couples.
3) The standard deviations defined in the parameters-Covariance*.py override the standard deviation defined in the observation or prior files.
4) Most of these files are optional. If there is no file for an certain type of observations, that means that there is no observation of this type. If a covariance matrix is not defined for a prior or an observation type, that means that the correlation matrix is supposed to be equal to identity and that the standard deviation is given in the prior or observation file.


 
What is the structure of the general `parameters.py` file?
--------------------------------------------------------

It contains the list of drillings, the optimization method to be used and some settings for the figures.
Have a look at the file `AICC2012-VLR/parameters.py`, it is commented.


What is the structure of a drilling `parameters.py` file?
---------------------------------------------------------

It defines age at the top of the core, the age equation grid and the correction functions grids. You can also define other parameters that are used to defined the covariance matrices of the priors.
Have a look at the files `AICC2012-VLR/EDC/parameters.py`, it is commented.


How to set up correlation matrices for the observations?
--------------------------------------------------------

You need to know a little bit of python to do that.
Feel free to send an email on the mailing list if you need assistance.

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
Feel free to send an email to the mailing list if you need assistance.

You need to define:
- `self.correlation_corr_a`     : the correlation matix for the accu correction function
- `self.correlation_corr_LID`   : the correlation matix for the LID correction function
- `self.correlation_corr_tau`   : the correlation matix for the thinning correction function
- `self.sigmap_corr_a`          : the standard deviation of the accu correction function
- `self.sigmap_corr_LID`        : the standard deviation of the LID correction function
- `self.sigmap_corr_tau`        : the standard deviation of the thinning correction function


