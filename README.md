############### Hybrid Kalman-Bayes filter ###########
Hybrid filter for post-process improvement of forecasts provided by NWP (Numerical Wind/Wave Prediction) models.
######################################################

Brief Manual
-------------

The main function of the program is HYBRID_KAL_BAYES_FILTER.M. All other functions are called from the main program and only one/two can operate as standalone but it is not recommended.
The program applies Hybrid Kalman-Bayes filter to a timeseries of forecasted values that require correction and produces output to the command line of MATLAB as well as to a txt file in the working folder.

--------------------------------------------------------------------------------------------------------------------

In order for this program to operate, two arrays must be entered. One for observations and one for model values (original forecast).
The function may receive further options to extend abilities or refine its functionality.

The general use of the program is the following:

The entire set of pairs OBS-MODEL should be entered to the program as:

	OUT = Hybrid_Kal_Bayes_filter(OBS,MODEL)
     
This is basically the automatic mode where default values of certain parameters are used and OUT is the resulting timeseries. Specifically, the first 72 values of both OBS (=OBS(1:72)) and MODEL (=MODEL(1:72)) are used as database in order to correct the following 24 values of MODEL (=MODEL(73:96)). There is a loop that passes through the entire length of the timeseries provided moving the window of 72 values with a step of 24. 
Practically, if one has a timeseries made by hourly values, history database consists of 3 days and the forecast concerns the 4th day.

Moreover, the output is either OUT only or a matrix with three rows such as [model_corr, model_orig, obs]
where model_orig and obs are the initial MODEL and OBS, respectively, but without their first 72 values. A text file is also created in the working directory containing the exact same results.

--------------------------------------------------------------------------------------------------------------------

Options/Refinements

- dim      ->   rank of the polynomial for Kalman, default value is {2} (works pretty good most of the time) 

- history    ->   change the length of the database, default value is {72}

- forecast    ->   change the length of forecast, default value is {24}

- kalmanIG    ->   If this option is given, KalmanIG is used instead of original Kalman filter. The value entered is passed as the multiplication factor when calculating variance in the training procedure. 

- Distribution  ->   Given this option, the program is forced to use the distribution provided instead of running tests to decide which one to use. Available distributions are: Weibull, Lognormal, Normal

- IntegSteps   ->    Change the accuracy of the integration during the bayesian steps. Default value is {100}. 
#!! Important Note !!# : Higher values increase runtime of the program dramatically!

- File     ->    With this option one can provide the exact path and name of the output text file.


Example with some available options and their default values:

[HYB   MODEL    OBS] = Hybrid_Kal_Bayes_filter(OBS,MODEL,’dim’, 2,’history’,72,’forecast’,24,’IntegSteps’,100);

- Limit   ->  This option must be followed by an array of two values indicating lower and higher values of indexes limiting the comparison figure.

[HYB   MODEL    OBS] = Hybrid_Kal_Bayes_filter(OBS,MODEL,’Limit’,[36 100]);

----------------------------------------------------------------------

Other examples:

[HYB   MODEL    OBS] = Hybrid_Kal_Bayes_filter(OBS,MODEL,’kalmanIG’,4.0);

If kalmanIG is given without a value following it, the program crashes. If kalmanIG is not given at all, classic Kalman filter is used.


[HYB   MODEL    OBS] = Hybrid_Kal_Bayes_filter(OBS,MODEL,’Distribution’,’weibull’);

[HYB   MODEL    OBS] = Hybrid_Kal_Bayes_filter(OBS,MODEL,’file’,’/PATH/TO/FILE.TXT’);



All of the options and the according arguments are NOT case sensitive, except from the path to the text file.

