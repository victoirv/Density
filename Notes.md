## History
[May 16: significance tests and pressure behavior](#may-16-2016)


### May 16, 2016 ###
The significance tests are now in the paper, and I've verified the t-tests as well as I can think to. I made two test datasets, each with 10 "events" of 100 time lags each and different means. Doing a t-test of 5 events of one mean vs 5 of the other returns 100% significant results, and doing 5 of one vs 5 more of the same mean returns 2-7% significant results, as expected.  The next thing I think I'd have to test is whether the variances are significantly different since the t-test assumes equal variances. 

The figure in question is this:
![Figure in question](paper/figures/PNGs/RhoBinnedBz-case24-t020-tf25-GOES6.png)

Where that green dot indicates significance, and seems like it shouldn't be there. I made a histogram of the two distributions (events with larger Bz vs events with smaller Bz) at that one significant point, and they do look somewhat different (though hard to tell since I can't print figures with transparency):

![Example](paper/figures/PNGs/RhoBinnedBz-case24-t020-tf25-GOES6-histogram.png)

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

