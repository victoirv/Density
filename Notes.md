## History
[May 18: Differences](#may-18-2016)
[May 16: significance tests and pressure behavior](#may-16-2016)


### May 18, 2016 ###
Trying to run the "differences" code. Compiles, but segfaults when running "make pcdiffvtk". Opens the data, interpolates, then crashes on writing output. Added mkdir -p output/Precondition/$(B)_minus_$(A) to the makefile which fixed it. That said, it still doesn't make vtk files (not sure if it's supposed to?)

"make images" fails with multiple errors of being unable to find non-specified results (e.g. Brian_Curtis_042213_1 when I only have _2 and _6), but might actually partially work if I can get the vtk files generated. Still investigating.

Also paraview file has disappeared from mag, so I'm replacing it (just kidding, don't have /var/www/tmp or sudo permissions), but might also change the link from mag to [the actual paraview website](http://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v4.2&type=binary&os=linux64&downloadFile=ParaView-4.2.0-Linux-64bit.tar.gz) in the makefile.

Created [bootstrap significance table for the daily binned Dst events](tables/DeltaBootstraps-case13.txt), showing a significant difference between median value on onset day vs the day after. See figure:
![Significant daily difference](figures/PNGs/DailyBootstrapDifferences-GOES6-case13.png)


### May 16, 2016 ###
The significance tests are now in the paper, and I've verified the t-tests as well as I can think to. I made two test datasets, each with 10 "events" of 100 time lags each and different means. Doing a t-test of 5 events of one mean vs 5 of the other returns 100% significant results, and doing 5 of one vs 5 more of the same mean returns 2-7% significant results, as expected.  The next thing I think I'd have to test is whether the variances are significantly different since the t-test assumes equal variances. 

The figure in question is this:
![Figure in question](paper/figures/PNGs/RhoBinnedBz-case24-t020-tf25-GOES6x.png)

Where that green dot indicates significance, and seems like it shouldn't be there. I made a histogram of the two distributions (events with larger Bz vs events with smaller Bz) at that one significant point, and they do look somewhat different (though hard to tell since I can't print figures with transparency):

![Example](paper/figures/PNGs/RhoBinnedBz-case24-t020-tf25-GOES6-histogram.png)

It was determined that the oddities were caused by significant differences in means and medians of the data: 

![Mean vs Median](paper/figures/RhoBinned/PNGs/MeanvsMedian.png)

This indicates that we can't use t-tests for differences in means, but must pursue bootstrapping (or some other median-based significance test such as Mann-Whitney (via ranksum function)) to determine actual levels of significance, since switching everything to means at this point would be arduous. Though these notes now probably look silly because the figures they link to now show the correct significance...

* * *

As for how Bz, Dst, and pressure behave during mass density events:

![storm with pressure](paper/figures/PNGs/stormavs-mass-GOES6-withPressure.png)

And again but with a lower mass density threshhold for events:

![storm with pressure smaller thresh](paper/figures/PNGs/stormavs-mass-gt20-GOES6-withPressure.png)

Can also look at the binned figures, but only a couple show any significant features:

* Bz binned by Dst
![Binned plots](figures/PNGs/HighLowDstBz-rhoeq20-GOES6-1983-1991.png)

* Bz binned by F10.7
![Binned plots](figures/PNGs/HighLowF107Bz-rhoeq20-GOES6-1983-1991.png)

* Dst binned by F10.7
![Binned plots](figures/PNGs/HighLowF107Dst-rhoeq20-GOES6-1983-1991.png)

* Mass Density binned by F10.7
![Binned plots](figures/PNGs/HighLowF107rhoeq-rhoeq20-GOES6-1983-1991.png)

Pressure doesn't seem to have much effect, though I also haven't looked at all possible binning combinations yet.


