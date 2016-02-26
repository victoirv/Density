%Show ten events to make point about data availability
if(MakePaperPlots && stormcase==1)
    figure('Visible',visible);
    xevents=(-12:1:24);
    idx=arrayfun(@colon,starti(1:5)-12,starti(1:5)+24,'Uniform',false);
    plot(repmat(xevents,5,1)',reshape(MassDensitySpline([idx{:}]),37,5),'r')
    hold on;
    idx2=arrayfun(@colon,starti(end-4:end)-12,starti(end-4:end)+24,'Uniform',false);
    plot(repmat(xevents,5,1)',reshape(MassDensitySpline([idx2{:}]),37,5),'b')
    ylabel('\rho_{eq} (amu/cm^3)')
    xlabel('Hours from onset')
    title(sprintf('First 5 (red) and last 5 (blue) events from %d-%d',sy,ey));
    print -depsc2 -r200 paperfigures/TenEvents.eps
    print -dpng -r200 paperfigures/PNGs/TenEvents.png
end

%Compare the two densities
if(MakePlots)
    figure('Visible',visible); plot(FILLEDTime,MassDensitySpline);hold on; 
    plot(FILLEDTime,OMNIDensity,'r')
    legend('Denton','OMNI','Location','NorthEast')
    title('OMNI Density vs Denton Density')
    ylabel('Density')
    xlabel('Time')
    datetick
    print -dpng figures/densitycomp.png
    print -depsc2 figures/densitycomp.eps
end

%Showing 'detrending' by removing F10.7 influencex`
if(MakePaperPlots && removef107)
    figure('Visible',visible);
    plot(FILLEDTime,MassDensitySplineOriginal,'r.')
    hold on; plot(FILLEDTime,MassDensitySpline,'b.')
    legend('Original','F_{10.7} Removed');
    ylabel('\rho_{eq} (amu/cm^3)')
    xlabel('Year')
    datetick
    print -depsc2 -r200 paperfigures/f107removed.eps
end