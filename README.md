IceChrono
=========

A statistical and physical model to optimize chronologies of deep polar ice cores.


What this manual is and is not?
-------------------------------

This manual is a documentation on how to use the IceChrono software.  
It is _not_ a description of the IceChrono principles. Please read to the scientific articles describing IceChrono for that purpose:  
Parrenin, F., Bazin, L., Capron, E., Landais, A., Lemieux-Dudon, B. and Masson-Delmotte, V.: IceChrono1: a probabilistic model to compute a common and optimal chronology for several ice cores, _Geosci. Model Dev._, 8(5), 1473–1492, doi:10.5194/gmd-8-1473-2015, 2015.  
It is _not_ an operating system or python documentation. Please use your operating system or python documentation instead.


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

IceChrono is a scientific python software, therefore you need a scipy distribution.  
IceChrono is developed and tested using the anaconda distribution, therefore we recommend it.  
Anaconda can be downloaded here (use the python2 version):  
http://continuum.io/downloads  

IceChrono probably works on other scipy distributions, provided they contain the following python modules:  
- sys
- os
- time
- math
- numpy
- matplotlib
- multiprocessing
- warnings
- scipy


How to run IceChrono?
---------------------

Assuming you use anaconda, you can go in the spyder shell and type the following commands:

```
cd path-to-IceChrono
import sys
sys.argv=['IceChrono.py','exp_directory']
execfile('IceChrono.py')
```

where `path-to-IceChrono` is the directory containing IceChrono and `exp_directory` is the name of your experiment directory. 
The `AICC2012-VLR` experiment directory is provided for you convenience. It is an AICC2012-like experiment, albeit whith a Very Low Resolution. It takes about 5 mn to run on a recent computer.

Alternatively, you can run IceChrono from your operating system shell. For example, on linux, enter the following commands in bash:

```
cd path-to-IceChrono
python IceChrono.py exp_directory
```


What are the outputs of a run:
------------------------------

If the run went correctly, it has created output files.

In the main directory, you have the following output file:
- `output.txt`		: only contains the program execution time for now.

In each drilling directory, you have the following output files:
- `output.txt`			: is the main output file. It gives you the posterior estimates and uncertainties of the three input variables (accumulation, LID and thinning) and of the output variables (ice age, air age, Δdepth, etc.). The header in the file tells you which column is which variable.
- `restart.txt`			: is a restart file, which can be used to start an optimization experiment from the result of a previous optimization experiment, for a faster convergence.
- `accumulation.pdf`	: is the accumulation figure (with prior and posterior estimates)
- `Ddepth.pdf`			: is the Δdepth figure (with prior estimates, observations and posterior estimates)
- `air_age.pdf`			: is the air age figure (with prior estimates, observations and posterior estimates)
- `airlayerthick.pdf`	: is the air layer thickness figure (with prior estimates, observations of dated intervals and posterior estimates)
- `ice_age.pdf`			: is the ice age figure (with prior estimates, observations and posterior estimates)
- `icelayerthick.pdf`	: is the ice layer thickness figure (with prior estimates, observations of dated intervals and posterior estimates)
- `LID.pdf`				: is the Lock-In Depth figure (with prior estimates, observations and posterior estimates)
- `thinning.pdf`		: is the thinning figure (with prior and posterior estimates)

In each drilling-pair directory, you have the following output files:
- `air-air.pdf`		: is the air-air stratigraphic links figure (with prior and posterior estimates)
- `air-ice.pdf`		: is the air-ice stratigraphic links figure (with prior and posterior estimates)
- `ice-air.pdf`		: is the ice-air stratigraphic links figure (with prior and posterior estimates)
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
You need to specify your prior scenarios for accumulation, LID and thinning (for each ice core) and your age observations.

You have five general files:
- `parameters.py`                                           : contains general parameters for the experiment
- `parameters-AllDrillings.py`                              : defines drilling parameters that are the same for all drillings (there are overidded by drilling specific parameters).
- `parameters-CovariancePrior-AllDrillings-init.py`         : defines the covariance matrices of the background
- `parameters-CovarianceObservations-AllDrillings.py`       : defines the covariance of the observations that are the same for all drillings  (there are overidded by drilling specific parameters).
- `parameters-CovarianceObservations-AllDrillingPairs.py` : defines the covariance for the observations that are the same for all drilling pairs  (there are overidded by drilling pair specific parameters).

Then you have one directory per drilling, which contains:
- `parameters.py`                           : all the drilling specific parameters
- `parameters-CovarianceObservations.py`    : this file allows to define the correlation of drilling specific observations
- `density-prior.txt`                       : depth / relative density
- `accu-prior.txt`                          : depth / background accu (in ice-equivalent) / standard deviation (opt, in %)
- `LID-prior.txt`                           : depth / background Lock-in-Depth / standard deviation (opt, in %)
- `thinning-prior.txt`                      : depth / background thinning function / standard deviation (opt, in %)
- `ice_age.txt`                             : depth / age / sigma for ice age observations
- `air_age.txt`                             : depth / age / sigma for air age observations
- `Ddepth.txt`                              : depth / Delta-depth / sigma for Delta-depth observations
- `ice_age_intervals.txt`		   			: depth\_top / depth\_bottom / duration / sigma for dated ice intervals observations
- `air_age_intervals.txt`   		    	: depth\_top / depth\_bottom / duration / sigma for dated air intervals observations

Then you have one directory per drilling pair, which contains:
- `parameters-CovarianceObservations.py`    : this file allows to define the correlation of drilling pair specific observations
- `ice_depth.txt`           : depth1 / depth2 / sigma on age for ice-ice stratigraphic links observations
- `air_depth.txt`           : depth1 / depth2 / sigma on age for air-air stratigraphic links observations
- `iceair_depth.txt`        : depth1 / depth2 / sigma on age for ice-air stratigraphic links observations
- `airice_depth.txt`        : depth1 / depth2 / sigma on age for air-ice stratigraphic links observations

A few things you need to know to use Icechrono:
1) You can use whatever units you want but they need to be consistent. For example, if you use meters for the depths and years for the dated horizons, you need to use meters per years for the accumulation rates. 
2) The drilling specific parameters override the general parameters for all drillings. In the very same way, the drilling-pair specific parameters override the general parameters for all drilling-pairs.
3) The standard deviations defined in the parameters-Covariance*.py override the standard deviation defined in the observation or prior files.
4) Most of these files are optional. If there is no file for an certain type of observations, that means that there is no observation of this type. If a covariance matrix is not defined for a prior or an observation type, that means that the correlation matrix is supposed to be equal to identity and that the standard deviation is given in the prior or observation file.


What is the structure of the general `parameters.py` file?
--------------------------------------------------------

It contains the list of drillings, the optimization method to be used and some settings for the figures.
It is where you define the names of your drillings.
Have a look at the file `AICC2012-VLR/parameters.py`, it is commented.


What is the structure of a drilling `parameters.py` file?
---------------------------------------------------------

It defines age at the top of the core, the unthinned depth at the top of the core, the age equation grid, the correction functions grids and the type of representation of the prior accu scenario (linear or staircase). You can also define other parameters that are used to defined the covariance matrices of the priors.
Have a look at the files `AICC2012-VLR/EDC/parameters.py`, it is commented.


How to set up the `parameters-CovarianceObservations.py` file?
--------------------------------------------------------------

You need to know a little bit of python to do that.
Feel free to send an email on the mailing list if you need assistance.

For drilling specific observations, you set up the correlation matrices in the file `parameters-CovarianceObservations.py` in the drilling directory.
- `self.icemarkers_correlation`     : for ice dated horizons
- `self.airmarkers_correlation`     : for air dated horizons
- `self.iceintervals_correlation`   : for ice dated intervals
- `self.airintervals_correlation`   : for air dated intervals
- `self.Ddepth_correlation`         : for Delta-depth observations

For drilling pair specific observations (stratigraphic links), you set up the correlation matrices in the file `parameters-CovarianceObservations.py` in the drilling pair directory.
- `self.iceicemarkers_correlation`  : for ice-ice stratigraphic links
- `self.airairmarkers_correlation`  : for air-air stratigraphic links
- `self.iceairmarkers_correlation`  : for ice-air stratigraphic links
- `self.airicemarkers_correlation`  : for air-ice stratigraphic links

Let us take a concrete example and assume we want a correlation matrix for ice dated horizons with ones in the diagonal and with a constant correlation factor k outside the diagonal, you can write:

```
self.icemarkers_correlation=k*np.ones((np.shape(self.icemarkers_correlation)))+(1-k)*np.diag(np.ones(np.shape(self.icemarkers_correlation)[0]))
```

Don't forget that if you find the use of python and the IceChrono internal variables too difficult, you can define your correlation matrices outside IceChrono and import them here by using for example the `np.loadtxt` function.


How to set up the `parameters-CovariancePrior-AllDrillings-init.py` file?
-------------------------------------------------------------------------

You need to know a little bit of python to do that and also to know a bit of the IceChrono internal variables.
Have a look at the `AICC2012-VLR` experiment, it is the easiest way to understand how it works.
Feel free to send an email to the mailing list if you need assistance.

You need to define:
- `self.correlation_corr_a`     : the correlation matix for the accu correction function
- `self.correlation_corr_LID`   : the correlation matix for the LID correction function
- `self.correlation_corr_tau`   : the correlation matix for the thinning correction function

Optionnally, if they have not been imported in the accu-prior.txt, LID-prior.txt and thinning-prior.txt files, you can also define:
- `self.sigmap_corr_a`          : the standard deviation of the accu correction function
- `self.sigmap_corr_LID`        : the standard deviation of the LID correction function
- `self.sigmap_corr_tau`        : the standard deviation of the thinning correction function


Let us take a concrete example and assume, as in the AICC2012 example, that the accumulation correlation linearly decreases to zero when the absolute value of the age difference of the accumulation corrections increases to lambda_a yr. We first define a matrix whose term (i,j) is equal to the age difference of the accumulation corrections `self.corr_a_age[i]-self.corr_a_age[j]`:

```
M=np.ones((np.size(self.corr_a_age),np.size(self.corr_a_age)))*self.corr_a_age-np.transpose(np.ones((np.size(self.corr_a_age),np.size(self.corr_a_age)))*self.corr_a_age)
```

And then we define the correlation matrix by interpolation a function that linearly decreases from 1 to 0 between 0 and self.lambda_a:

```
self.correlation_corr_a=np.interp(np.abs(M), np.array([0,self.lambda_a]),np.array([1, 0]))
```

Don't forget that if you find the use of python and the IceChrono internal variables too difficult, you can define your correlation matrices outside IceChrono and import them here by using for example the `np.loadtxt` function.


What to do if something goes wrong?
-----------------------------------

Please post an email on the mailing list with the error message appearing on the command line.
