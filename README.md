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

A table of correlations (prediction with 1 impulse coefficient vs with 12 coefficients) follows:

<pre>
Input 	 CC1 	 CC12 	 PE(1) 	 PE(12)
W3	 0.00 	 0.07 	 0.00 	 0.01
Pdyn	 0.01 	 0.04 	 0.00 	 0.00
ByIMF	 0.01 	 0.03 	 0.00 	 0.00
Day	 0.02 	 0.20 	 0.00 	 0.04
dst	 0.02 	 0.05 	 0.00 	 0.00
G1	 0.02 	 0.07 	 0.00 	 0.00
DBS	 0.04 	 0.05 	 0.00 	 0.00
6stat	 0.04 	 0.04 	 0.00 	 0.00
8stat	 0.05 	 0.05 	 0.00 	 0.00
G3	 0.06 	 0.06 	 0.00 	 0.00
W6	 0.06 	 0.06 	 0.00 	 0.00
W4	 0.06 	 0.07 	 0.00 	 0.00
Den_P	 0.07 	 0.08 	 0.01 	 0.01
BS	 0.08 	 0.09 	 0.01 	 0.01
W1	 0.08 	 0.08 	 0.01 	 0.01
G2	 0.08 	 0.08 	 0.01 	 0.01
VBS	 0.08 	 0.09 	 0.01 	 0.01
W2	 0.08 	 0.09 	 0.01 	 0.01
Hr	 0.09 	 0.21 	 0.01 	 0.05
W5	 0.09 	 0.09 	 0.01 	 0.01
kp	 0.10 	 0.13 	 0.01 	 0.02
akp3	 0.11 	 0.15 	 0.01 	 0.02
DBz	 0.12 	 0.13 	 0.01 	 0.02
VBz	 0.12 	 0.16 	 0.02 	 0.03
BzIMF	 0.14 	 0.20 	 0.02 	 0.04
Bz1	 0.14 	 0.20 	 0.02 	 0.04
Bz2	 0.14 	 0.20 	 0.02 	 0.04
Bz3	 0.14 	 0.20 	 0.02 	 0.04
Bz4	 0.14 	 0.20 	 0.02 	 0.04
Bz5	 0.14 	 0.20 	 0.02 	 0.04
Bz6	 0.14 	 0.20 	 0.02 	 0.04
Hr+Bz	 0.16 	 0.28 	 0.03 	 0.08
V_SW	 0.19 	 0.20 	 0.04 	 0.04
Bz+V	 0.23 	 0.28 	 0.05 	 0.08
Year	 0.36 	 0.36 	 0.13 	 0.13
lnF107	 0.43 	 0.43 	 0.18 	 0.19
F107	 0.43 	 0.43 	 0.19 	 0.19
F107+V	 0.44 	 0.44 	 0.20 	 0.20
F107+Bz	 0.45 	 0.45 	 0.20 	 0.20

</pre>