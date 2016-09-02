* ~~Line 37: I suggest "larger/more positive 1.5 and 3-day averages".~~

* ~~Line 41: "and and" -> "and"~~

* ~~Line 70: Delete "in the".~~

* ~~Line 106: -> "after A period of"~~

* Line 129, "Dst events correspond to elevated rho_m before and after the event": I find this very confusing. Isn't your data showing that rho_m is elevated on the day of the event? **Fixed**

* Lines 133-134, "twice the standard deviation of the values used in computing the median divided by the square root of the number of values": I suggest that you give some motivation for this. I suspect that you are trying to characterize the uncertainty of the median value. **Fixed**

* Line 137, "as all available measurements were used instead of restricting to only values where a rho_m value also existed": Doesn't this introduce the danger that the rho_m values could correspond to different values of the indices? **Assuming he means "Could these indices not actually correspond to an onset", yes, but they'd be close, and if anything will trend towards less significance**

* Lines 130-131 and 144: You mentioned that you found an increase of 4.5 amu/cm^-3 where as Takahashi et al. found an increase of 10 amu/cm^-3. Could this be because of the different spacecraft used? You mentioned that GOES 2 found a much larger increase **GOES 2 doesn't cover 1989-1991, GOES 7 figure:**
![fig](paper/figures/PNGs/stormavs-dst-50-tak-GOES7.png)

* Line 151, "compute the median on each epoch day with replacement": To me, it's not exactly clear what you're doing here. Are you replacing the values by the measured values plus or minus a value based on the standard deviation? **Rephrase as "created by sampling (with replacement) the values used..."**  **Fixed**

* ~~Line 168, -> "which ARE shown in the top panel".~~

* ~~Line 180: -> "there IS statistically significant variation".~~

* Lines 185-187, "were separated into two parts ... including the hour of onset": The description of what you are doing here is extremely unclear and needs to be written. Is this what you mean: "depending on whether the 4 hour average for a particular event was less than or greater than the 4 hour average averaged over all events"? **"seperated into two parts based on the median value within a window. One window covered the four hours before onset and onset, and the other covered onset and the four hours after."**

* ~~Line 204, -> "were taken during"~~

* Lines 214-215, "whereas the 10 quiet Kp events considered by Denton et al. [2016] showed growth that lasted for at least 48 hours": Note that the events in my paper were extremely special meeting very strict criteria, a period of very low Kp following a period of not low Kp. Also, a better way to compare to my paper might be a superposed epoch analysis using minimum rather than maximum rho_m.

* Line 217, "The largest difference occurs when the separation is performed after onset": Seems like this would almost certainly be the case because the large rho_m starts at onset, so if you're comparing the range of large rho_m to something, this would give you the strongest effect. **Not sure what he means**

* ~~Line 219, "for in the time windows used": The sentence is not understandable. I'm guessing that what you mean is that if you look at the rate of increase in rho_m preceding the onset, you see similar rates of increase.~~

* ~~Line 222, "of if the": I suggest "of whether the".~~

* Line 225, "but similar peak values": The following paper notes that there seems to be an upper saturation limit for rho_m of about 300 amu/cm^3 regardless of whether the spacecraft is in the plasmasphere or plasmatrough:

Denton, R.E., K. Takahashi, M. F. Thomsen, J. E. Borovsky, H. J.
Singer, Y. Wang, J. Goldstein, P. C. Brandt, and B. W. Reinisch
(2014), Evolution of mass density and O+ concentration at
geostationary orbit during storm and quiet Events, J. Geophys. Res.,
119 (8), doi: 10.1002/2014JA019888.

* Line 227, "because the initial starting rho_m is higher for high F10.7, the associated growth rate is lower". Denton et al [2016] made this point also. See their Figure 9.

* Line 332, references: you might also consider adding: Denton, R. E., M. F. Thomsen, K. Takahashi, R. R. Anderson, and H. J. Singer (2011), Solar cycle dependence of bulk ion composition at geosynchronous orbit, J. Geophys. Res., 116, A03212, doi:10.1029/2010JA016027.

* Line 246, "statistically week**(sic)** enhancement": See the model in Denton et al. [2016].

* Line 257-258, "more positive values of B_z are associated with larger rho_m values": This was noted for the model of Denton et al. [2016]. We speculated that it might be due to the closed magnetosphere being more conducive for refilling.

* ~~Line 264: -> "may BE due".~~

* Line 265, "refilling mechanism described by Denton et al. [2016]": Do you mean the point I was making an item 22 **(two bullets up)** above? It's not clear here what you mean.

* Figure 1: A log scale would be better for the last panel showing rho_m. Densities vary by enormous factors. **Not particularly relevant**

* Figure 3: I see that there is a trend for F10.7 to be slightly decreasing over the time interval plotted. Could that be because these events are more likely during the declining phase of the solar cycle? 

* Figure 4: You could consider going to several hour averages rather than one hour, for better statistics.

* I noticed that the density was larger before the decrease in Dst. Maybe this is due to higher Dst during this interval or due to the "calm before the storm".

* The "calm before the storm'' in CIR/magnetosphere interactions: Occurrence statistics, solar wind statistics, and magnetospheric preconditioning By: Borovsky, JE; Steinberg, JT JOURNAL OF GEOPHYSICAL RESEARCH-SPACE PHYSICS   Volume: 111   Issue: A7     Article Number: A07S10   Published: JUN 16 2006

* How significant are the details of this pattern for rho_m at the bottom? Higher -> drop -> peak -> rebound -> second peak -> drop. What do you think it means?

* ~~All figures: The fonts are too small to easily read.~~ **Fixed? Made larger, at least**

* ~~Figure 7: No significant variation in other parameters. That probably means that there's a lot of apparently random variation.~~ 

* ~~Figure 8: I see a repeated peak about 24 hours after the initial peak, especially in the right plot, but also evident in the left plot. I think this is probably because high mass density structures corotate around Earth and you see them again 24 hours later.~~

