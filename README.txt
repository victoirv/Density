This study was done to determine the effectiveness of linear impulse response modeling on mass density at the magnetic equator. The data to be predicted is from the "GOES Alfven Frequency and Mass Density Data 1980-1991" located here:
http://www.dartmouth.edu/~rdenton/

It's cleaned up by changing fill values (usually 9999, but every column uses something different) to NaNs for uniformity's sake. The file I use is his on the website but with all of the headers removed. dlmread() supports this in the function, but I did it a while ago and just used that edited file here.

The impulses are taken from OMNI spacecraft measurements, but since that data is very sparse in the years of interest, a reconstructed dataset is taken from this paper:

http://onlinelibrary.wiley.com/doi/10.1002/2014GL059741/suppinfo

specifically, WGhourFS_72_13.txt. The data are then formatted from white space separated to comma separated in vim using %s/\s+/,/g and the leading comma removed with %s/^,//g . Though MATLAB has means of reading in white space delimited files, none worked well and this ended up being easier and quicker.

By looking at GOES spacecraft 6, the years 1983-1992 can be covered for prediction and modeling.

All variables from the filled OMNI dataset were used as impulses with both 1 impulse coefficient, and 120 impulse coefficients, to determine how much difference a complete linear model made in predictability.

The table of correlation coefficients is in table.txt, and the plots of all variables (in the format OMNI_{var}.png) as well as all predicted datasets, and their corresponding coefficients (format density_{Var}_{num persist coef}_{num impulse coef}.png), are in the figures/ folder.

This study so far seems to indicate that there is a value to using a linear impulse model for prediction. The statistical significance of any increase remains to be tested.
