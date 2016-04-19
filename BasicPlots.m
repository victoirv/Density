
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


if(MakePaperPlots && stormcase==1)
    
    F107d=interptest(1:length(MassDensitySpline),FILLED(:,30),1:24*3:length(MassDensitySpline));
    MDd=interptest(1:length(MassDensitySpline),MassDensitySpline,1:24*3:length(MassDensitySpline));
    
    %Sliding window smoothing (instead of block averaging)
    %MDd2=ndnanfilter((1/73*ones(1,73)),1,MassDensitySpline);
    MDd2=ndnanfilter(MassDensitySpline,1,(1/73*ones(1,73)));
    F107d2=filter((1/73*ones(1,73)),1,FILLED(:,30));
    %figure; plot((0:(length(F107d)-1))*24*3,F107d)
    %hold on; plot(F107d2,'r')
    
    MDd27=ndnanfilter(MassDensitySpline,1,(1/649*ones(1,649)));
    F107d27=filter((1/649*ones(1,649)),1,FILLED(:,30));
    [cx, cf, ~,xnew,corr] = IR(MDd27,F107d27,0,12,0,0);
    
    
    %General (non-storm) trend
    [cx, cf, ~,xnew,corr] = IR(MDd2,F107d2,0,12,0,0);
    
    h=figure('Visible',visible);
    hold on; 
    plot(0:3:33,cf,'b');
    grid on
    xlabel('Time Lags (day)')
    ylabel('Impulse Response Coefficient')
    title('Predicting \rho_{eq} with F_{10.7}')
    print('-depsc2', '-r200', sprintf('paperfigures/F107IR-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('paperfigures/PNGs/F107IR-GOES%d.png',satnum));
    if(strcmp(visible,'off')),close(h);end;
    
    
    
    XOnsets=MDd;
    XOnsets(setdiff(1:end,floor(starti/24)))=NaN;
    [cx, cf, ~,xnew,corr] = IR(XOnsets,F107d,0,12,0,0);
    
        h=figure('Visible',visible);
    hold on; 
    plot(0:3:33,cf,'b');
    grid on
    xlabel('Time Lags (day)')
    ylabel('Impulse Response Coefficient')
    title('Predicting \rho_{eq} at onset of D_{st} event using F_{10.7} ')
    print('-depsc2', '-r200', sprintf('paperfigures/F107IR-onset-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('paperfigures/PNGs/F107IR-onset-GOES%d.png',satnum));
    if(strcmp(visible,'off')),close(h);end;
    
    
end