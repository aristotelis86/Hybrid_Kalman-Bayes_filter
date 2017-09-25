# Hybrid Kalman-Bayes filter 

Hybrid filter for post-process improvement of forecasts provided by Numerical Wind/Wave Prediction (NWP) models.


## Overview

The main function of the program is HYBRID_KAL_BAYES_FILTER.M. All other functions are called from the main program. Running other scripts and/or functions individually is not recommended.
The program applies Hybrid Kalman-Bayes filter to a timeseries of forecasted values that require correction and produces output to the command line of MATLAB as well as to a txt file in the working folder.

### How to call it

In order for this program to operate, two arrays must be entered. One for observations and one for model values (original forecast).
So the basic usage would be:

	OUT = Hybrid_Kal_Bayes_filter(OBS,MODEL)
	
which creates the corrected time series OUT based on default parameter values. Specifically, the first 72 values of both OBS and MODEL are used as database in order to correct the following 24 values of MODEL. 
The database is changing using a 72 values moving window with a step of 24 values.
Practically, if one has a timeseries made by hourly values, history database consists of 3 days and the forecast concerns the 4th day.

### Further functionaliy

The function may receive further options to extend abilities or refine its functionality.
The format of the output is either OUT only or a matrix with three rows such as \[model_corr, model_orig, obs\]
where model_orig and obs are the initial MODEL and OBS, respectively, but without their first 72 values. A text file is also created in the working directory containing the exact same results.
This is done mainly for validating the model.

Options/Refinements
---------------------

- dim      ->   rank of the polynomial for Kalman, default value is {2} (works pretty good most of the time) 

- history    ->   change the length of the database, default value is {72}

- forecast    ->   change the length of forecast, default value is {24}

- kalmanIG    ->   If this option is given, KalmanIG is used instead of original Kalman filter. The value entered is passed as the multiplication factor when calculating variance in the training procedure. 

- Distribution  ->   Given this option, the program is forced to use the distribution provided instead of running tests to decide which one to use. Available distributions are: Weibull, Lognormal, Normal

- IntegSteps   ->    Change the accuracy of the integration during the bayesian steps. Default value is {100}. 
#!! Important Note !!# : Higher values increase runtime of the program dramatically!

- File     ->    With this option one can provide the exact path and name of the output text file.

- Limit   ->  This option must be followed by an array of two values indicating lower and higher values of indexes limiting the comparison figure. (For visualization of the results)


Example calls with some available options and their default values:

	[HYB,   MODEL,    OBS] = Hybrid_Kal_Bayes_filter(OBS,MODEL,’dim’, 2,’history’,72,’forecast’,24,’IntegSteps’,100);
	
	[HYB,   MODEL,    OBS] = Hybrid_Kal_Bayes_filter(OBS,MODEL,’Limit’,[36 100]);

Other examples:

	[HYB   MODEL    OBS] = Hybrid_Kal_Bayes_filter(OBS,MODEL,’kalmanIG’,4.0);

If kalmanIG is given without a value following it, the program crashes. A coefficient for the variation of the data must be applied. If kalmanIG is not given at all, classic Kalman filter is used.


	[HYB   MODEL    OBS] = Hybrid_Kal_Bayes_filter(OBS,MODEL,’Distribution’,’weibull’);

	[HYB   MODEL    OBS] = Hybrid_Kal_Bayes_filter(OBS,MODEL,’file’,’/PATH/TO/FILE.TXT’);



All of the options and the according arguments are NOT case sensitive, except from the path to the text file, which is OS-dependent.

