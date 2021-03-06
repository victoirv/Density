This study was done to determine the effectiveness of linear impulse response modeling on mass density at the magnetic equator. The data to be predicted is from the "GOES Alfven Frequency and Mass Density Data 1980-1991" located here:
http://www.dartmouth.edu/~rdenton/

Also described in the paper here:
http://onlinelibrary.wiley.com/doi/10.1029/2009JA015243/abstract

It's cleaned up by changing fill values (usually 9999, but every column uses something different) to NaNs for uniformity's sake. The file I use is his on the website but with all of the headers removed. dlmread() supports this in the function, but I did it a while ago and just used that edited file here.

The solar wind measurements used as impulses are taken from a gap-filled dataset here:

http://onlinelibrary.wiley.com/doi/10.1002/2014GL059741/suppinfo

specifically, WGhourFS_72_13.txt. The data are then formatted from white space separated to comma separated in vim using %s/\s+/,/g and the leading comma removed with %s/^,//g and fixing the two status headers (8 status and 6 stat, both end up with a comma that they shouldn't have). Though MATLAB has means of reading in white space delimited files, none worked well and this ended up being easier and quicker.

By looking at GOES spacecraft 6, the years 1983-1992 can be covered for prediction and modeling.

All variables from the filled OMNI dataset were used as impulses with both 1 impulse coefficient, and 12 impulse coefficients (equivalent to 12 hours), to determine how much difference a complete linear model made in predictability.

The table of correlation coefficients is in table.txt, and the plots of all variables (in the format OMNI_{var}.png) as well as all predicted datasets, and their corresponding coefficients (format density_{Var}_{num persist coef}_{num impulse coef}.png), are in the figures/ folder.

This study so far seems to indicate that there is a value to using a linear impulse model for prediction. The statistical significance of any increase remains to be tested, but it certainly seems that the best linear predictor of solar wind density is Bz, and not just Bs. The next few things to look into are significance tests, nonlinear predictions, and a comparison of the prediction efficiency of large vs small Bz, and Bz with periods of precoditioning. 

Some interesting things to note from the following table (added here for ease of concatenation): 
-For 1 hour time steps, Bz outperforms B_south in predicting mass density. 
-Though F10.7 is still a significant predictor, it's not nearly as effective at 1 hour increments as it is on 27 day averages
-Looking at figures/avf107.png shows that the variables seem to have a strong dependence on hour. However, stripping the data to only predict mass density after 10AM UT only seems to affect models that include hour.

The table of correlations (prediction with 1 impulse coefficient vs with 12 coefficients) follows:

