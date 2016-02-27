

if(MakePaperPlots && stormcase==1)
    stackplot(FILLEDTime, [FILLED(:,[5,6,15,30]) MassDensitySpline'],{'B_z (nT)','V_{SW} (km/s)','D_{st} (nT)','F_{10.7} (s.f.u.)','\rho_{eq} (amu/cm^3)'},satnum,[3 DSTCut])
    
    if(satnum==6 && sy<=1989 && ey>=1989) %Make March 1989 Geomagnetic storm plot (if necessary time range still exists)
        stackplot(FILLEDTime, [FILLED(:,[5,6,15,30]) MassDensitySpline'],{'B_z (nT)','V_{SW} (km/s)','D_{st} (nT)','F_{10.7} (s.f.u.)','\rho_{eq} (amu/cm^3)'},satnum,[],[],[datenum('Mar-10-1989') datenum('Mar-18-1989')])
    end

    h=figure('Visible',visible);
    plot(-timewidth:1:timewidth*2,AVMDMat,'.')
    hold on;
    xa=(-timewidth:LongTimeScale:timewidth*2);
    plot(xa,AVMDs(1,:),'r','LineWidth',3);
    hold on; plot(xa,AVMDMatBars(:,:),'r-.','LineWidth',2);
    print -depsc2 -r200 paperfigures/allstorms.eps
    print -dpng -r200 paperfigures/PNGs/allstorms.png    
    if(strcmp(visible,'off')),close(h);end;
    
    
    h=figure('Visible',visible);    
    for i=1:length(AVMDMat(1,:))
        top(i)=nanmean(AVMDMat(AVMat(:,i,5)>0,i));
        topbar(i)=nanstd(AVMDMat(AVMat(:,i,5)>0,i));
        bottom(i)=nanmean(AVMDMat(AVMat(:,i,5)<0,i));
        bottombar(i)=nanstd(AVMDMat(AVMat(:,i,5)<0,i));
    end
    
    plot(xa,top)
    hold on; 
    plot(xa,bottom,'r')
    plot(xa,[top+topbar; top-topbar],'b-.')
    plot(xa,[bottom+bottombar; bottom-bottombar],'r-.')
    grid on;
    legend('B_z > 0','B_z < 0')
    ylabel('\rho_{eq} (amu/cm^3)')
    xlabel('Time since onset (hours)')
    print('-depsc2',sprintf('paperfigures/RhoBinnedBz-%d.eps',satnum));
    print('-dpng', '-r200', sprintf('paperfigures/PNGs/RhoBinnedBz-%d.png',satnum));  
    
    
    
    h=figure('Visible',visible);
    hist(nanmean(AVMDMat(:,20:24),2),0:5:100)
    h1 = findobj(gca,'Type','patch');
    %set(h1,'FaceColor','r','EdgeColor','w','facealpha',0.75)
    hold on;
    hist(nanmean(AVMDMat(:,25:28),2),0:5:100)
    h2 = findobj(gca,'Type','patch');
    %set(h2,'facealpha',0.75);
    legend('Before and at Onset','After Onset')
    xlim([-5,105])
    xlabel('Average \rho_{eq} (amu/cm^3)')
    ylabel('Count')
    title('Average \rho_{eq} Four Hours Before and After Storm Onset')
    %Printing segfaults for some reason...
    print('-depsc2',sprintf('paperfigures/rhobeforeafter-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('paperfigures/PNGs/rhobeforeafter-GOES%d.png',satnum));    
    if(strcmp(visible,'off')),close(h);end;
    
    
    
    h=figure('Visible',visible);
    for i=1:100
        before(i)=nanmean(randsample(reshape(AVMDMat(:,21:25),[],1),4));
        after(i)=nanmean(randsample(reshape(AVMDMat(:,26:29),[],1),4));
    end
    hist(before,0:5:100)
    hold on;
    hist(after,0:5:100)
    legend('Before and at Onset','After Onset')
    xlim([-5,105])
    xlabel('Average \rho_{eq} (amu/cm^3)')
    ylabel('Count')
    title('Average \rho_{eq} Four Hours Before and After Storm Onset')
    %Printing segfaults for some reason...
    print('-depsc2',sprintf('paperfigures/rhobeforeafter-boot-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('paperfigures/PNGs/rhobeforeafter-boot-GOES%d.png',satnum));    
    if(strcmp(visible,'off')),close(h);end;
    
    
    
end

if(MakePaperPlots && MakeBinPlots)
    if(DSTCut<0)
        plotthresh=sprintf('< %d',DSTCut);
        stormtype='D_{st}'; stormunits='nT';
    elseif(AECut>0)
        plotthresh=sprintf('> %d',AECut);
        stormtype='AE'; stormunits='nT';
    elseif(stormcase==17)
        plotthresh='';
        stormtype='random'; stormunits='';
    else
        plotthresh=sprintf('> %d',MDCut);
        stormtype='\rho_{eq}'; stormunits='amu/cm^3';
    end
    
    %plot dst, sort by f10.7
    binplot(AVMat(:,:,15),FILLED(:,30),starti,timewidth,LongTimeScale,plotthresh,{'D_{st}';'F_{10.7}';stormtype},{'nT';'s.f.u';stormunits},[sy; ey],satnum,visible);
    %plot rho, sort f10.7
    binplot(AVMDMat(:,:),FILLED(:,30),starti,timewidth,LongTimeScale,plotthresh,{'\rho_{eq}';'F_{10.7}';stormtype},{'amu/cm^3';'s.f.u';stormunits},[sy; ey],satnum,visible);
    %plot Bz, sort f10.7
    binplot(AVMat(:,:,5),FILLED(:,30),starti,timewidth,LongTimeScale,plotthresh,{'B_z';'F_{10.7}';stormtype},{'nT';'s.f.u';stormunits},[sy; ey],satnum,visible);
    %plot Bz, sort rho
    binplot(AVMat(:,:,5),MassDensitySpline,starti,timewidth,LongTimeScale,plotthresh,{'B_z';'\rho_{eq}';stormtype},{'nT';'amu/cm^3';stormunits},[sy; ey],satnum,visible);
    %plot rho, sort dst
    binplot(AVMDMat(:,:),FILLED(:,15),starti,timewidth,LongTimeScale,plotthresh,{'\rho_{eq}';'D_{st}';stormtype},{'amu/cm^3';'nT';stormunits},[sy; ey],satnum,visible);
    %plot Bz, sort dst
    binplot(AVMat(:,:,5),FILLED(:,15),starti,timewidth,LongTimeScale,plotthresh,{'B_z';'D_{st}';stormtype},{'nT';'nT';stormunits},[sy; ey],satnum,visible);
end

if(MakePaperPlots && MakeDstThreshPlot)
    stormtype='D_{st}'; stormunits='nT';
    DCs=-30:-5:-90;
    s=zeros(1,length(DCs));
    st=s;
    numevents=s;
    DCi=1;
    for DC=DCs
        storms=diff([0 (FILLED(:,15)<DC)' 0]);
        starti=find(storms>0);
        endi=find(storms<0)-1;
        stormi=1;
        AVMDMat=[];
        while(starti(1)-maxwidth<1)
            starti(1)=[]; endi(1)=[];
        end
        while(endi(end)+maxwidth>length(MassDensitySpline))
            starti(end)=[]; endi(end)=[];
        end
        for i=1:length(starti)
            AVMDMat(stormi,:)=MassDensitySpline((starti(i)-timewidth):starti(i)+timewidth*2);
            stormi=stormi+1;
        end
        plotthresh=sprintf('< %d',DC);
        [s(DCi),st(DCi)]=twobinplot(AVMDMat(:,:),FILLED(:,30),starti,timewidth,LongTimeScale,plotthresh,{'\rho_{eq}';'F_{10.7}';stormtype},{'amu/cm^3';'s.f.u';stormunits},[sy; ey],satnum,visible); 
        numevents(DCi)=length(starti);
        DCi=DCi+1;
        
    end
    h=figure('visible',visible);
    plot(DCs,s,'r')
    hold on; plot(DCs,s+st,'r-.'); plot(DCs,s-st,'r-.');
    %plot(DCs,log10(numevents),'b-.')
    plot(DCs,zeros(1,length(DCs)),'k-','LineWidth',3)
    %legend('\rho_{eq}','log(# events)');
    xlabel('D_{st} threshold (nT)');
    ylabel('\rho_{eq} (amu/cm^3)')
    title(sprintf('High F_{10.7} events - Low F_{10.7} events, both with baselines removed, %d-%d',sy,ey))
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on')
    print('-depsc2','-r200', sprintf('paperfigures/DstRhoThresh-%d-%d.eps',sy,ey));
    print('-dpng','-r200', sprintf('paperfigures/PNGs/DstRhoThresh-%d-%d.png',sy,ey));
    if(strcmp(visible,'off')),close(h);end;
end


if(MakePaperPlots && MakeRandThreshPlot)    
    stormtype='Random'; stormunits='units';
    DCs=800:-100:300;
    s=zeros(1,length(DCs));
    st=s;
    numevents=s;
    DCi=1;
    for DC=DCs
        starti = randsample(1:length(MassDensitySpline),DC,false);
        storms=zeros(1,length(MassDensitySpline));
        storms(starti)=1; storms(starti+1)=-1;
        while(sum(storms)~=0)
            starti = randsample(1:length(MassDensitySpline),DC,false);
            storms=zeros(1,length(MassDensitySpline));
            storms(starti)=1; storms(starti+1)=-1;
        end
        starti=find(storms>0);
        endi=find(storms<0)-1;
        stormi=1;
        AVMDMat=[];
        while(starti(1)-maxwidth<1)
            starti(1)=[]; endi(1)=[];
        end
        while(endi(end)+maxwidth>length(MassDensitySpline))
            starti(end)=[]; endi(end)=[];
        end
        for i=1:length(starti)
            AVMDMat(stormi,:)=MassDensitySpline((starti(i)-timewidth):starti(i)+timewidth*2);
            stormi=stormi+1;
        end
        plotthresh=sprintf('< %d',DC);
        [s(DCi),st(DCi)]=twobinplot(AVMDMat(:,:),FILLED(:,30),starti,timewidth,LongTimeScale,plotthresh,{'\rho_{eq}';'F_{10.7}';stormtype},{'amu/cm^3';'s.f.u';stormunits},[sy; ey],satnum,visible); 
        numevents(DCi)=length(starti);
        DCi=DCi+1;
        
    end
    h=figure('visible',visible);
    plot(DCs,s,'r')
    hold on; plot(DCs,s+st,'r-.'); plot(DCs,s-st,'r-.');
    %plot(DCs,log10(numevents),'b-.')
    plot(DCs,zeros(1,length(DCs)),'k-','LineWidth',3)
    %legend('\rho_{eq}','log(# events)');
    xlabel('Random Sample Size');
    ylabel('\rho_{eq} (amu/cm^3)')
    title(sprintf('High F_{10.7} events - Low F_{10.7} events, both with baselines removed, %d-%d',sy,ey))
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on')
    print('-depsc2','-r200', sprintf('paperfigures/RandRhoThresh-%d-%d.eps',sy,ey));
    print('-dpng','-r200', sprintf('paperfigures/PNGs/RandRhoThresh-%d-%d.png',sy,ey));
    if(strcmp(visible,'off')),close(h);end;
end

%Make main stack plots
if(MakePaperPlots)
    xa=(-timewidth:LongTimeScale:timewidth*2)./LongTimeScale;
    h=figure('Visible',visible);
    orient tall;
    
    h1=subplot('position',subplotstack(5,1));plot(xa,AVs(:,5),'+-','LineWidth',2); text(0.01,0.9,'B_z (nT)','Units','normalized','FontSize',14); %ylabel('B_z (nT)'); %Bz
    hold on; plot(xa,AVMatBars(:,:,5),'r-.'); ylim(yranges(yr,1,:))
    h2=subplot('position',subplotstack(5,2));plot(xa,AVs(:,6),'+-','LineWidth',2); text(0.01,0.9,'V_{SW} (km/s)','Units','normalized','FontSize',14); %ylabel('V_{SW} (km/s)');%V_sw
    hold on; plot(xa,AVMatBars(:,:,6),'r-.'); ylim(yranges(yr,2,:))
    h3=subplot('position',subplotstack(5,3));plot(xa,AVs(:,15),'+-','LineWidth',2); text(0.01,0.9,'D_{st} (nT)','Units','normalized','FontSize',14); %ylabel('D_{st} (nT)'); %dst
    hold on; plot(xa,AVMatBars(:,:,15),'r-.'); ylim(yranges(yr,3,:))
    if(DSTCut<0), plot([xa(1) xa(end)],[DSTCut DSTCut],'k-.','LineWidth',2); hold off; end
    h4=subplot('position',subplotstack(5,4));plot(xa,AVs(:,30),'+-','LineWidth',2); text(0.01,0.9,'F_{10.7} (s.f.u)','Units','normalized','FontSize',14); %ylabel('F10.7 (s.f.u.)');%f107
    hold on; plot(xa,AVMatBars(:,:,30),'r-.'); ylim(yranges(yr,4,:))
    set(findobj('type','axes'),'xticklabel',{[]})
    xv=[xa(1) xa(end)];
    subplot('position',subplotstack(5,5)); [AX,H5,H6]=plotyy(xa,AVMDs(1,:),xa,AVnnans,'plot','bar');
    hold on; plot(xa,AVMDMatBars(:,:),'r-.'); %ylim(AX(1),[0 90])
    if(MDCut>0), plot([xa(1) xa(end)],[MDCut MDCut],'k-.','LineWidth',2); end
    set(H5,'LineWidth',2);set(AX(2),'Xlim',xv); set(AX(1),'Xlim',xv);  set(H5,'marker','+','color','red'); set(AX(1),'YColor','r'); set(AX(2),'YColor',[0 0.5 0.5]); set(get(H6,'child'),'FaceColor',[0 0.5 0.5]); uistack(AX(1)); set(AX(1),'Color','none'); set(AX(2),'Color','w');
    text(0.01,0.85,'\rho_{eq} (amu/cm^3)','Units','normalized','FontSize',14); %ylabel(AX(1),'\rho_{eq} (amu/cm^3)'); 
    ylabel(AX(2),'# valid hourly values');
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[-timewidth:timewidth/2:timewidth*2]./LongTimeScale)
    linkaxes([AX h1 h2 h3 h4],'x')
    if(DSTCut<0), title(h1,sprintf('%d D_{st} < %d nT events; %d-%d',length(starti),DSTCut,sy,ey)); end
    if(AECut<0), title(h1,sprintf('%d AE > %d events; %d-%d',length(starti),AECut,sy,ey)); end
    if(MDCut>0), title(h1,sprintf('%d \\rho_{eq} > %d amu/cm^3 events; %d-%d',length(starti),MDCut,sy,ey)); end
    if(LongTimeScale>1),xlabel('Time from start of event (day)');
    else xlabel('Time from start of event (hour)'); end
    fprintf('Average of %d storms %s (%d to %d)\n',length(duration),durationcaveat, sy,ey);
    print('-depsc2','-r200',figurename);
    if(strcmp(visible,'off')),close(h);end;
end



%DST vs rho_eq for 1 hour and 1 day
if(MakePaperPlots && stormcase==1) 
    F107Day=interptest(FILLEDTime,FILLED(:,30),FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end));
    MDDay=interptest(FILLEDTime,MassDensitySpline',FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end));
    F10727Day=interptest(FILLEDTime,FILLED(:,30),FILLEDTime(1):24*27*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end));
    MD27Day=interptest(FILLEDTime,MassDensitySpline',FILLEDTime(1):24*27*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end));
    
    h=figure('Visible',visible);
    
    plot(FILLED(~isnan(MassDensitySpline),30),NaN.*log10(MassDensitySpline(~isnan(MassDensitySpline))),'.','MarkerSize',15);
    hold on;
    plot(F107Day,NaN.*log10(MDDay),'r.','MarkerSize',20);
    plot(F10727Day,NaN.*log10(MD27Day),'go','LineWidth',1.5,'MarkerSize',12) ;
    
    plot(FILLED(~isnan(MassDensitySpline),30),log10(MassDensitySpline(~isnan(MassDensitySpline))),'.','MarkerSize',3);
    hold on;
    plot(F107Day,log10(MDDay),'r.','MarkerSize',8);    
    plot(F10727Day,log10(MD27Day),'go','LineWidth',1.5,'MarkerSize',12) ;
    xlabel('F_{10.7} (s.f.u.)','FontSize',16)
    
    ylabel('log_{10}[\rho_{eq} (amu/cm^3)]','FontSize',16)
    
    cc1=corrcoef(FILLED(~isnan(MassDensitySpline),30),log10(MassDensitySpline(~isnan(MassDensitySpline))),'rows','pairwise');
    cc1=cc1(1,2);
    cc2=corrcoef(F107Day,log10(MDDay),'rows','pairwise');
    cc2=cc2(1,2);
    cc3=corrcoef(F10727Day,log10(MD27Day),'rows','pairwise');
    cc3=cc3(1,2);
    
    legend(sprintf('1-Hour Median; cc = %2.2f',cc1),sprintf('1-Day Medians; cc = %2.2f',cc2),sprintf('27-Day Medians; cc = %2.2f',cc3),'Location','SouthEast');

    legend boxoff;
    axis tight;
    print('-depsc2',sprintf('paperfigures/ccplot-GOES%d.eps',satnum)); 
    print('-dpng','-r200',sprintf('paperfigures/PNGs/ccplot-GOES%d.png',satnum)); 
    if(strcmp(visible,'off')),close(h);end;
    
 
    
    
    
    
end