
%Variables histogrammed against MLT
if(MakePaperPlots && stormcase==1)
    avrhos=zeros(1,24);
    for i=1:24
        avrhos(i)=nanmedian(MassDensity(round(MLT)==(i-1)));
    end
    figure('Visible',visible);
    plot(0:23,avrhos,'+-')
    ylabel('Median \rho_{eq} (amu/cm^3)')
    xlabel('MLT (hour)')
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[0:2:24])
    print -depsc2 -r200 paperfigures/rhoMLT.eps
    print -dpng -r200 paperfigures/PNGs/rhoMLT.png
    if(strcmp(visible,'off')),close(h);end;
    
    [~,~,~,DHr]=datevec(DentonTime);
    for i=1:24
        avrhos(i)=nanmedian(MassDensity(DHr==(i-1)));
    end
    h=figure('Visible',visible);
    plot(0:23,avrhos,'+-')
    ylabel('Median \rho_{eq} (amu/cm^3)')
    xlabel('Local Time (hour)')
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[0:4:23])
    print -depsc2 -r200 paperfigures/rhoLT.eps
    print -dpng -r200 paperfigures/PNGs/rhoLT.png
    if(strcmp(visible,'off')),close(h);end;
    
    for i=1:24
        avrhos(i)=nanmedian(AE(DHr==(i-1)));
    end
    figure('Visible',visible);
    plot(0:23,avrhos,'+-')
    ylabel('Median AE (nT)')
    xlabel('Local Time (hour)')
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[0:4:23])
    print -depsc2 -r200 paperfigures/AELT.eps
    print -dpng -r200 paperfigures/PNGs/AELT.png
    if(strcmp(visible,'off')),close(h);end;
    
    for i=1:24
        avrhos(i)=nanmedian(FILLED(FILLED(:,3)==(i-1),15));
    end
    figure('Visible',visible);
    plot(0:23,avrhos,'+-')
    ylabel('Median D_{st} (nT)')
    xlabel('Local Time (hour)')
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[0:4:23])
    print -depsc2 -r200 paperfigures/DstLT.eps
    print -dpng -r200 paperfigures/PNGs/DstLT.png
    if(strcmp(visible,'off')),close(h);end;
    

%Nans per hour
    h=figure('Visible',visible);
    hist(FILLED(isnan(MassDensitySpline),3),0:23)
    axis([-1 24 0 3000])
    set(gca,'XTick',0:2:24)
    xlabel('UT Hour of event start')
    ylabel('Frequency')
    print -depsc2 -r200 paperfigures/nansbyhour.eps
    print -dpng -r200 paperfigures/PNGs/nansbyhour.png
    if(strcmp(visible,'off')),close(h);end;

    h=figure('Visible',visible);
    subplot(2,1,1)
    plot(xa,AVs(:,3))
    axis tight;
    ylabel('Hour average')
    grid on
    subplot(2,1,2)
    plot(xa,AVnnans,'r')
    axis tight;
    set(gca,'XTick',0:2:24)
    ylabel('Data available')
    xlabel('Time from event start (hr)')
    grid on
    print -depsc2 -r200 paperfigures/nansbyhour_storm.eps
    print -dpng -r200 paperfigures/PNGs/nansbyhour_storm.png
    if(strcmp(visible,'off')),close(h);end;
end


if(MakePaperPlots && stormcase==16)
    h=figure('Visible',visible);
    hist(FILLED(starti,3),0:23)
    axis([-1 24 0 100])
    grid on
    xlabel('UT Hour')
    ylabel('Frequency')
    title(sprintf('%d events of AE > %d for %d-%d',length(starti),AECut,sy,ey));
    print -depsc2 -r200 paperfigures/AEbyhour.eps
    print -dpng -r200 paperfigures/PNGs/AEbyhour.png
    if(strcmp(visible,'off')),close(h);end;
end

if(MakePaperPlots && stormcase==1)
    h=figure('Visible',visible);
    hold on; 
    plot(FILLED(:,2),FILLED(:,15),'b.');
    plot(FILLED(starti,2),FILLED(starti,15),'r.','MarkerSize',10);
    legend({'All Data','Event Onset'})
    grid on
    xlabel('Day of Year')
    ylabel('D_{st} (nT)')
    print -depsc2 -r200 paperfigures/DoYDst.eps
    print -dpng -r200 paperfigures/PNGs/DoYDst.png
    if(strcmp(visible,'off')),close(h);end;
end