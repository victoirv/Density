
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


if(MakePaperPlots && stormcase==27)
    
    New1dTime=FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end);
    New3dTime=FILLEDTime(1):24*3*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end);
    New27dTime=FILLEDTime(1):24*27*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end);
    
    F107d=interptest(FILLEDTime,FILLED(:,30),New1dTime);
    MDd=interptest(FILLEDTime,MassDensitySpline,New1dTime);
    
    %Sliding window smoothing (instead of block averaging)
    %MDd2=ndnanfilter((1/73*ones(1,73)),1,MassDensitySpline);
    MDd2=ndnanfilter(MassDensitySpline,1,(1/73*ones(1,73)));
    F107d2=filter((1/73*ones(1,73)),1,FILLED(:,30));
    %figure; plot((0:(length(F107d)-1))*24*3,F107d)
    %hold on; plot(F107d2,'r')
    
    
    F1073d=interptest(FILLEDTime,FILLED(:,30),New3dTime);
    MD3d=interptest(FILLEDTime,MassDensitySpline,New3dTime);
    [cx, cf, cc,xnew,corr] = IR(log(MD3d),F1073d,0,12,0,0);
    
    
    
    figure; subplot(2,1,1); plot(New3dTime,log(MD3d)); hold on; plot(New3dTime,xnew,'r'); legend('Data','Model')
    datetick('keeplimits'); xlabel('Date'); ylabel('\rho_{eq} (amu/cm^3)');
    
    subplot(2,1,2);  plot(0:11,flipud(cf)); hold on; plot([0 11],[0 0],'k-.');
        ylabel('Impulse Response coefficient')
    xlabel('Lags (3 day)')
    title(sprintf('Coefficients for predicting \rho_{eq} with 12 3-day lags of F_{10.7} - CC: %2.2f',corr))
        print('-depsc2', '-r200', sprintf('figures/F1073dCoef-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('figures/F1073dCoef-GOES%d.png',satnum));
    
    
    
    F10727d=interptest(FILLEDTime,FILLED(:,30),New27dTime);
    MD27d1=interptest(FILLEDTime,MassDensitySpline,New27dTime);
    [cx, cf, ~,xnew,corr1] = IR(log(MD27d1),F10727d,0,12,0,0);
    
    figure; subplot(2,1,1); plot(New27dTime,log(MD27d1)); hold on; plot(New27dTime,xnew,'r'); legend('Data','Model')
    datetick('keeplimits'); xlabel('Date'); ylabel('\rho_{eq} (amu/cm^3)');
        subplot(2,1,2);  plot(0:11,flipud(cf)); hold on; plot([0 11],[0 0],'k-.');
        ylabel('Impulse Response coefficient')
    xlabel('Lags (27 day)')
    title(sprintf('Coefficients for predicting \rho_{eq} with 12 3-day lags of F_{10.7} - CC: %2.2f',corr1))
        print('-depsc2', '-r200', sprintf('figures/F10727dCoef-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('figures/F10727dCoef-GOES%d.png',satnum));
    
    %Verification
    for i=4:length(F10727d)
       MD27dv(i)=exp(1*F10727d(i-1)-0.5*F10727d(i-2)+0.2*F10727d(i-3));
    end
        [cx, cf, ~,xnew,corrv] = IR(log(MD27dv),F10727d,0,12,0,0);
    
    figure; subplot(2,1,1); plot(New27dTime,log(MD27dv)); hold on; plot(New27dTime,xnew,'r'); legend('Data','Model')
    datetick('keeplimits'); xlabel('Date'); ylabel('\rho_{eq} (amu/cm^3)');
        subplot(2,1,2);  plot(0:11,flipud(cf)); hold on; plot([0 11],[0 0],'k-.');
        ylabel('Impulse Response coefficient')
    xlabel('Lags (27 day)')
    title(sprintf('Coefficients for predicting \rho_{eq} with 12 3-day lags of F_{10.7} - CC: %2.2f',corrv))
        print('-depsc2', '-r200', sprintf('figures/F10727dVerif-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('figures/F10727dVerif-GOES%d.png',satnum));
    
    
    
    %648 = 24*27
    nterms=648;
    nbins=fix(length(MassDensitySpline)/nterms);
    MDd272=nanmean(reshape(MassDensitySpline(1:nbins*nterms),nterms,nbins));
    F107d27=nanmean(reshape(FILLED(1:nbins*nterms,30),nterms,nbins));
    [cx, cf, ~,xnew,corr2] = IR(log(MDd272),F107d27,0,12,0,0);

    MDd273=nanmedian(reshape(MassDensitySpline(1:nbins*nterms),nterms,nbins));
    F107d273=nanmedian(reshape(FILLED(1:nbins*nterms,30),nterms,nbins));
    [cx, cf, ~,xnew,corr3] = IR(log(MDd273),F107d273,0,12,0,0);
    
    MDd274=ndnanfilter(MassDensitySpline,@rectwin,649);
    F107d27=filter((1/649*ones(1,649)),1,FILLED(:,30));
    F107d274=ndnanfilter(FILLED(:,30),@rectwin,649);
    [cx, cf, ~,xnew,corr4] = IR(log(MDd274),F107d274,0,12,0,0);
    
    sprintf('Mine \treshape\t median\t filter\n %2.2f \t %2.2f \t %2.2f\t %2.2f \n',corr1, corr2,corr3,corr4)
    
    figure; plot(FILLEDTime, MDd274,'r','LineWidth',3); hold on; plot(New27dTime,MD27d1','LineWidth',3); plot(New27dTime(1:end-1),MDd272 ,'k','LineWidth',3);  plot(New27dTime(1:end-1),MDd273 ,'g','LineWidth',3);
    datetick;
    ylabel('\rho_{eq} (amu/cm^3)')
    xlabel('Date')
    legend('Filter','Mine','Mean Reshape','Median Reshape')
        print('-depsc2', '-r200', sprintf('figures/InterpStyles-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('figures/InterpStyles-GOES%d.png',satnum));
    
    
    %General (non-storm) trend
    [cx, cf, ~,xnew,corr] = IR(MDd2,F107d2,0,12,0,0);
    
    h=figure('Visible',visible);
    hold on; 
    plot(0:3:33,flipud(cf),'b');
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
    plot(0:3:33,flipud(cf),'b');
    grid on
    xlabel('Time Lags (day)')
    ylabel('Impulse Response Coefficient')
    title('Predicting \rho_{eq} at onset of D_{st} event using F_{10.7} ')
    print('-depsc2', '-r200', sprintf('paperfigures/F107IR-onset-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('paperfigures/PNGs/F107IR-onset-GOES%d.png',satnum));
    if(strcmp(visible,'off')),close(h);end;
    
    
end