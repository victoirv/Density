This study was done to determine the effectiveness of linear impulse response modeling on mass density at the magnetic equator. The data to be predicted is from the "GOES Alfven Frequency and Mass Density Data 1980-1991" located here:
http://www.dartmouth.edu/~rdenton/

Also described in the paper here:
http://onlinelibrary.wiley.com/doi/10.1029/2009JA015243/abstract

It's cleaned up by changing fill values (usually 9999, but every column uses something different) to NaNs for uniformity's sake. The file I use is his on the website but with all of the headers removed. dlmread() supports this in the function, but I did it a while ago and just used that edited file here.

The solar wind measurements are taken from a gap-filled dataset here:

http://onlinelibrary.wiley.com/doi/10.1002/2014GL059741/suppinfo

specifically, WGhourFS_72_13.txt. The data are then formatted from white space separated to comma separated in vim using %s/\s+/,/g and the leading comma removed with %s/^,//g . Though MATLAB has means of reading in white space delimited files, none worked well and this ended up being easier and quicker.

By looking at GOES spacecraft 6, the years 1983-1992 can be covered for prediction and modeling.

All variables from the filled OMNI dataset were used as impulses with both 1 impulse coefficient, and 120 impulse coefficients, to determine how much difference a complete linear model made in predictability.

The table of correlation coefficients is in table.txt, and the plots of all variables (in the format OMNI_{var}.png) as well as all predicted datasets, and their corresponding coefficients (format density_{Var}_{num persist coef}_{num impulse coef}.png), are in the figures/ folder.

This study so far seems to indicate that there is a value to using a linear impulse model for prediction. The statistical significance of any increase remains to be tested.

A table of correlations (prediction with 1 impulse coefficient vs with 120 coefficients) follows:

Variable 	 corr(1) 	 corr(120)
Year 	 0.40649 	 0.40873
Day 	 0.02062 	 0.04038
Hr 	 0.00546 	 0.04759
ByIMF 	 -0.01842 	 -0.02475
BzIMF 	 0.08010 	 0.18258
V_SW 	 0.09221 	 0.12063
Den_P 	 0.05621 	 0.10010
Pdyn 	 0.04620 	 0.07635
G1 	 0.05513 	 0.00313
G2 	 -0.00568 	 0.02843
G3 	 -0.01177 	 0.03263
8stat 	 0.02706 	 0.03505
kp 	 0.01068 	 0.08803
akp3 	 0.02380 	 0.08993
dst 	 -0.02562 	 0.00855
Bz1 	 0.08010 	 0.18258
Bz2 	 0.08010 	 0.18258
Bz3 	 0.08010 	 0.18258
Bz4 	 0.08010 	 0.18258
Bz5 	 0.08010 	 0.18258
Bz6 	 0.08010 	 0.18258
W1 	 0.00007 	 0.04389
W2 	 0.01350 	 0.06754
W3 	 -0.02017 	 0.00518
W4 	 -0.02006 	 0.01150
W5 	 0.00084 	 0.05064
W6 	 -0.01140 	 0.03591
6stat 	 0.01938 	 0.02466
VBS 	 -0.00430 	 0.03222
