
%Variables histogrammed against MLT
if(MakePaperPlots && stormcase==1)
    avrhos=zeros(1,24);
    for i=1:24
        avrhos(i)=nanmedian(MassDensity(round(MLT)==(i-1)));
    end
    figure('Visible',visible);
    plot(0:23,avrhos,'+-')
    axis([-1 24 5 35])
    ylabel('Median \rho_{eq} (amu/cm^3)')
    xlabel('MLT (hour)')
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[0:2:24])
    print -depsc2 -r200 figures/rhoMLT.eps
    print -dpng -r200 figures/PNGs/rhoMLT.png
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
    print -depsc2 -r200 figures/rhoLT.eps
    print -dpng -r200 figures/PNGs/rhoLT.png
    if(strcmp(visible,'off')),close(h);end;
    
    for i=1:24
        avrhos(i)=nanmedian(AE(DHr==(i-1)));
    end
    figure('Visible',visible);
    plot(0:23,avrhos,'+-')
    ylabel('Median AE (nT)')
    xlabel('Local Time (hour)')
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[0:4:23])
    print -depsc2 -r200 figures/AELT.eps
    print -dpng -r200 figures/PNGs/AELT.png
    if(strcmp(visible,'off')),close(h);end;
    
    for i=1:24
        avrhos(i)=nanmedian(FILLED(FILLED(:,3)==(i-1),15));
    end
    figure('Visible',visible);
    plot(0:23,avrhos,'+-')
    ylabel('Median D_{st} (nT)')
    xlabel('Local Time (hour)')
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[0:4:23])
    print -depsc2 -r200 figures/DstLT.eps
    print -dpng -r200 figures/PNGs/DstLT.png
    if(strcmp(visible,'off')),close(h);end;
    
    
    %Nans per hour
    h=figure('Visible',visible);
    hist(FILLED(isnan(MassDensitySpline),3),0:23)
    axis([-1 24 0 3000])
    set(gca,'XTick',0:2:24)
    xlabel('UT Hour')
    ylabel('Frequency')
    print -depsc2 -r200 figures/nansbyhour.eps
    print -dpng -r200 figures/PNGs/nansbyhour.png
    if(strcmp(visible,'off')),close(h);end;
    
    %data per MLT
    h=figure('Visible',visible);
    %MLTinterp=FILLED(:,29);
    %MLTinterp(isnan(MLTinterp)) = interp1(find(~isnan(MLTinterp)), MLTinterp(~isnan(MLTinterp)), find(isnan(MLTinterp)), 'linear'); 
    %hist(MLTinterp,0:23)
    hist(MLT,0:23)
    axis([-1 24 0 8000])
    set(gca,'XTick',0:2:24)
    xlabel('MLT')
    ylabel('Count')
    print -depsc2 -r200 figures/databyMLT.eps
    print -dpng -r200 figures/PNGs/databyMLT.png
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
    print -depsc2 -r200 figures/nansbyhour_storm.eps
    print -dpng -r200 figures/PNGs/nansbyhour_storm.png
    if(strcmp(visible,'off')),close(h);end;
    
    
    
    
    h=figure('Visible',visible);
    hold on;
    for i=1:10
    plot(xa,i/3+normc(AVMat(i,:,5)'));
    end
    plot(xa,normc(nanmedian(AVMat(1:10,:,5))'),'r','LineWidth',2);
    %plot(xa,normc(AVs(:,5)),'r','LineWidth',2);
    ylabel('Normalized D_{st}')
    xlabel('Hours since onset')
    grid on;
    print -depsc2 -r200 figures/epochexample.eps
    print -dpng -r200 figures/PNGs/epochexample.png
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
    print -depsc2 -r200 figures/AEbyhour.eps
    print -dpng -r200 figures/PNGs/AEbyhour.png
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
    print -depsc2 -r200 figures/DoYDst.eps
    print -dpng -r200 figures/PNGs/DoYDst.png
    if(strcmp(visible,'off')),close(h);end;
end


if(MakePaperPlots && stormcase==27)
    
    New1dTime=FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end);
    New3dTime=FILLEDTime(1):24*3*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end);
    New27dTime=FILLEDTime(1):24*27*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end);
    
    F107d=interptest(FILLEDTime,FILLED(:,30),New1dTime);
    Dstd=interptest(FILLEDTime,FILLED(:,15),New1dTime);
    [MDd, NUsed]=interptest(FILLEDTime,MassDensitySpline,New1dTime);
    
    
    
    save('data/1dData','New1dTime','F107d','Dstd','MDd','NUsed')
    
    
    figure; plot(NUsed,MDd,'+')
    coef=[NUsed(~isnan(MDd)) ones(length(NUsed(~isnan(MDd))),1)]\MDd(~isnan(MDd));
    cc=corrcoef(MDd,NUsed,'rows','pairwise');
    hold on; plot([0 24],[coef(2) coef(2)+coef(1)*24],'r-.','LineWidth',3)
    ylabel('\rho_{eq}')
    xlabel('Valid hourly points in daily median')
    title(sprintf('Linear correlation: %2.2f',cc(1,2)));
    print('-depsc2', '-r200', sprintf('figures/MDvsValid-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('figures/MDvsValid-GOES%d.png',satnum));
    
    
    %Sliding window smoothing (instead of block averaging)
    %MDd2=ndnanfilter((1/73*ones(1,73)),1,MassDensitySpline);
    MDd2=ndnanfilter(MassDensitySpline,1,(1/73*ones(1,73)));
    F107d2=filter((1/73*ones(1,73)),1,FILLED(:,30));
    %figure; plot((0:(length(F107d)-1))*24*3,F107d)
    %hold on; plot(F107d2,'r')
    
    
    F1073d=interptest(FILLEDTime,FILLED(:,30),New3dTime);
    MD3d=interptest(FILLEDTime,MassDensitySpline,New3dTime);
    [cx, cf, cc,xnew,corr, ~, cxsd, cfsd, ccsd] = IR(log(MD3d),F1073d,0,12,0,0,100);
    
    figure; subplot(2,1,1); plot(New3dTime,log(MD3d)); hold on; plot(New3dTime,xnew,'r'); legend('Data','Model')
    datetick('keeplimits'); xlabel('Date'); ylabel('\rho_{eq} (amu/cm^3)');
    
    subplot(2,1,2);  plot(0:11,flipud(cf)); hold on;
    plot(0:11,[flipud(cf)+flipud(cfsd) flipud(cf)-flipud(cfsd)],'r-.')
    plot([0 11],[0 0],'k-.');
    xlim([0 11])
    ylabel('Impulse Response coefficient')
    xlabel('Lags (3 day)')
    title(sprintf('Average 100-bootstrap-sample coefficients for predicting \\rho_{eq} with 12 3-day lags of F_{10.7} - CC: %2.2f',corr))
    print('-depsc2', '-r200', sprintf('figures/F1073dCoef-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('figures/F1073dCoef-GOES%d.png',satnum));
    
    %Redo but with same time span as 27 day
    F1073d=interptest(FILLEDTime,FILLED(:,30),New3dTime);
    MD3d=interptest(FILLEDTime,MassDensitySpline,New3dTime);
    [cx, cf, cc,xnew,corr, ~, cxsd, cfsd, ccsd] = IR(log(MD3d),F1073d,0,100,0,0,100);
    
    F10727d=interptest(FILLEDTime,FILLED(:,30),New27dTime);
    MD27d1=interptest(FILLEDTime,MassDensitySpline,New27dTime);
    [cx, cf1, ~,xnew1,corr1, ~, cxsd, cfsd1, ccsd] = IR(log(MD27d1),F10727d,0,12,0,0,100);
    
    figure; subplot(2,1,1); plot(New27dTime,log(MD27d1)); hold on; plot(New27dTime,xnew1,'r'); plot(New3dTime,xnew,'k');
    legend('Data',sprintf('27-day model - CC: %2.2f',corr1),sprintf('3-day model - CC: %2.2f',corr))
    datetick('keeplimits'); xlabel('Date'); ylabel('\rho_{eq} (amu/cm^3)');
    subplot(2,1,2);  hold on; plot([0 27*11],[0 0],'m-.','LineWidth',3);
    plot(0:27:27*11,flipud(cf1),'r','LineWidth',3); 
    plot(0:27:27*11,[flipud(cf1)+flipud(cfsd1) flipud(cf1)-flipud(cfsd1)],'r-.')
    plot(0:3:27*11,flipud(cf),'k','LineWidth',3); 
    plot(0:3:27*11,[flipud(cf)+flipud(cfsd) flipud(cf)-flipud(cfsd)],'k-.')
    xlim([0 27*11])
    ylabel('Impulse Response coefficient')
    xlabel('Lags (days)')
    title('Average 100-bootstrap-sample coefficients for predicting \\rho_{eq} with F_{10.7}')
    print('-depsc2', '-r200', sprintf('figures/F10727dCoef-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('figures/F10727dCoef-GOES%d.png',satnum));
    
    
    [cx, cf2, ~,xnew2,corr2, ~, cxsd, cfsd2, ccsd] = IR(log(interptest(FILLEDTime,MassDensitySpline,FILLEDTime(1):24*13.5*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)))...
        ,interptest(FILLEDTime,FILLED(:,30),FILLEDTime(1):24*13.5*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end))...
        ,0,23,0,0,100);
    [cx, cf3, ~,xnew3,corr3, ~, cxsd, cfsd3, ccsd] = IR(log(interptest(FILLEDTime,MassDensitySpline,FILLEDTime(1):24*9*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)))...
        ,interptest(FILLEDTime,FILLED(:,30),FILLEDTime(1):24*9*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end))...
        ,0,34,0,0,100);
    
    figure; subplot(2,1,1); plot(New27dTime,log(MD27d1)); 
    hold on; 
    plot(New27dTime,xnew1,'r'); 
    plot(FILLEDTime(1):24*13.5*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end),xnew2,'k');
    plot(FILLEDTime(1):24*9*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end),xnew3,'g');
    
    legend('Data','27-day Model','13.5-day Model','9-day Model')
    datetick('keeplimits'); xlabel('Date'); ylabel('\rho_{eq} (amu/cm^3)');
    subplot(2,1,2);  hold on;  
    plot([0 27*11],[0 0],'m-.','LineWidth',3);
    h1=plot(0:27:27*11,flipud(cf1),'r-','LineWidth',3); plot(0:27:27*11,[flipud(cf1)+flipud(cfsd1) flipud(cf1)-flipud(cfsd1)],'r-.')
    h2=plot(0:13.5:27*11,flipud(cf2),'k-','LineWidth',3); plot(0:13.5:27*11,[flipud(cf2)+flipud(cfsd2) flipud(cf2)-flipud(cfsd2)],'k-.')
    h3=plot(0:9:27*11,flipud(cf3),'g-','LineWidth',3); plot(0:9:27*11,[flipud(cf3)+flipud(cfsd3) flipud(cf3)-flipud(cfsd3)],'g-.')
    legend([h1 h2 h3],{'27-day Model','13.5-day Model','9-day Model'})

    
    xlim([0 27*11])
    ylabel('Impulse Response coefficient')
    xlabel('Lags (days)')
    title('Average 100-bootstrap-sample coefficients for predicting \rho_{eq} with F_{10.7}')
    print('-depsc2', '-r200', sprintf('figures/F107MultiDayCoef-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('figures/F107MultiDayCoef-GOES%d.png',satnum));
    
    
    
    
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
    title(sprintf('Coefficients for verifying IR prediction - CC: %2.2f',corrv))
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
    print('-depsc2', '-r200', sprintf('figures/F107IR-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('figures/PNGs/F107IR-GOES%d.png',satnum));
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
    print('-depsc2', '-r200', sprintf('figures/F107IR-onset-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('figures/PNGs/F107IR-onset-GOES%d.png',satnum));
    if(strcmp(visible,'off')),close(h);end;
    
    
end


if(MakePaperPlots && stormcase==28 && sy==1989 && ey==1989) %For a short time period plot, make quick IR model plot
    
    VBS=(FILLED(:,6).*(abs(FILLED(:,5))-FILLED(:,5))./2)./100;
    [~,~, ~,xnew,corr] = IR(FILLED(:,15),VBS,0,12,0,0);
    
    h=figure('Visible',visible);
    orient tall;
    hold on;
    h(1)=subplot('position',subplotstack(2,1)); 
    
    plot(FILLEDTime,FILLED(:,15),'k','LineWidth',1.5); 
    hold on;
    plot(FILLEDTime,xnew,'b','LineWidth',1.5);
    
    [~,~, ~,xnew,corr2] = IR(FILLED(:,15),VBS,1,0,0,0); %persistence 
    
    plot(FILLEDTime,xnew,'r','LineWidth',1.5);
    text(0.01,0.1,'D_{st} (nT)','Units','normalized','FontSize',14);
    h_legend=legend('Measured',sprintf('IR Model - CC %2.2f',corr),sprintf('Persistence Model - %2.2f',corr2),'Location','SouthEast');
    set(h_legend,'FontSize',14)
    h(2)=subplot('position',subplotstack(2,2)); 
    
    plot(FILLEDTime,VBS,'k','LineWidth',1.5);
    text(0.01,0.92,'vB_S (V/Km)','Units','normalized','FontSize',14);
    set(findobj('type','axes'),'xticklabel',{[]});
set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on')
axis tight;
datetick('x','keeplimits')
set(findobj('type','axes'),'xtick',get(h(end),'xtick'))
linkaxes(h,'x')
axis tight;
    xlabel('Date of 1989')
    
    print('-depsc2', '-r200', sprintf('figures/BasicModelExample-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('figures/PNGs/BasicModelExample-GOES%d.png',satnum));
    %if(strcmp(visible,'off')),close(h);end;
    
    
    
    
   
    
    
end

if(MakePaperPlots && stormcase==28 && sy~=ey)
     %Persistence correlation plot
    lags=200;
    for i=0:lags
         [~,~, ~,~,corr(i+1)] = IR(FILLED(:,15),FILLED(:,9),1,0,i,0); %persistence
    end
    
    h=figure('Visible',visible);
    %figure
    plot(0:lags,corr,'LineWidth',2)
    ylabel('Correlation')
    xlabel('Lags in hours')
    print('-depsc2', '-r200', 'figures/PersistentCorrelation.eps');
    print('-dpng', '-r200', 'figures/PNGs/PersistentCorrelation.png');
    if(strcmp(visible,'off')),close(h);end;
    
end