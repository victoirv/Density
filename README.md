This study was done to determine the effectiveness of linear impulse response modeling on mass density at the magnetic equator. The data to be predicted is from the "GOES Alfven Frequency and Mass Density Data 1980-1991" located here:
http://www.dartmouth.edu/~rdenton/

Also described in the paper here:
http://onlinelibrary.wiley.com/doi/10.1029/2009JA015243/abstract

It's cleaned up by changing fill values (usually 9999, but every column uses something different) to NaNs for uniformity's sake. The file I use is his on the website but with all of the headers removed. dlmread() supports this in the function, but I did it a while ago and just used that edited file here.

The solar wind measurements are taken from a gap-filled dataset here:

http://onlinelibrary.wiley.com/doi/10.1002/2014GL059741/suppinfo

specifically, WGhourFS_72_13.txt. The data are then formatted from white space separated to comma separated in vim using %s/\s+/,/g and the leading comma removed with %s/^,//g and fixing the two status headers (8 status and 6 stat, both end up with a comma that they shouldn't have). Though MATLAB has means of reading in white space delimited files, none worked well and this ended up being easier and quicker.

By looking at GOES spacecraft 6, the years 1983-1992 can be covered for prediction and modeling.

All variables from the filled OMNI dataset were used as impulses with both 1 impulse coefficient, and 120 impulse coefficients, to determine how much difference a complete linear model made in predictability.

The table of correlation coefficients is in table.txt, and the plots of all variables (in the format OMNI_{var}.png) as well as all predicted datasets, and their corresponding coefficients (format density_{Var}_{num persist coef}_{num impulse coef}.png), are in the figures/ folder.

This study so far seems to indicate that there is a value to using a linear impulse model for prediction. The statistical significance of any increase remains to be tested.

A table of correlations (prediction with 1 impulse coefficient vs with 120 coefficients) follows:

Variable 	 corr(1) 	 corr(120) 	 eff(1) 	 eff(120)
Year 	 -0.00018 	 0.01011 	 -0.12737 	 -0.12499
Day 	 0.00154 	 0.03346 	 0.00000 	 0.00170
Hr 	 0.00601 	 0.05224 	 0.00004 	 0.00331
ByIMF 	 0.00370 	 0.00987 	 0.00001 	 0.00068
BzIMF 	 0.08715 	 0.19816 	 0.00749 	 0.03913
V_SW 	 0.07154 	 0.11575 	 0.00512 	 0.01397
Den_P 	 0.08113 	 0.11318 	 0.00658 	 0.01337
Pdyn 	 0.06386 	 0.07962 	 0.00407 	 0.00690
G1 	 0.03005 	 0.03768 	 0.00079 	 0.00169
G2 	 0.01095 	 0.06679 	 0.00010 	 0.00493
G3 	 0.00525 	 0.05628 	 0.00002 	 0.00371
8stat 	 0.03457 	 0.04296 	 0.00119 	 0.00242
kp 	 0.02326 	 0.11172 	 0.00054 	 0.01302
akp3 	 0.03872 	 0.11345 	 0.00150 	 0.01341
dst 	 0.04517 	 0.07900 	 0.00191 	 0.00668
Bz1 	 0.08715 	 0.19816 	 0.00749 	 0.03913
Bz2 	 0.08715 	 0.19816 	 0.00749 	 0.03913
Bz3 	 0.08715 	 0.19816 	 0.00749 	 0.03913
Bz4 	 0.08715 	 0.19816 	 0.00749 	 0.03913
Bz5 	 0.08715 	 0.19816 	 0.00749 	 0.03913
Bz6 	 0.08715 	 0.19816 	 0.00749 	 0.03913
W1 	 0.01283 	 0.07382 	 0.00015 	 0.00595
W2 	 0.02212 	 0.09056 	 0.00049 	 0.00873
W3 	 0.01593 	 0.03660 	 0.00018 	 0.00185
W4 	 0.00529 	 0.04793 	 -0.00001 	 0.00281
W5 	 0.01081 	 0.07604 	 0.00011 	 0.00630
W6 	 0.00700 	 0.05228 	 0.00005 	 0.00331
6stat 	 0.04376 	 0.05429 	 0.00191 	 0.00351
VBS 	 0.01070 	 0.06891 	 0.00010 	 0.00522
