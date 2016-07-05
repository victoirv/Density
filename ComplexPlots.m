

if(MakePaperPlots && stormcase==28)
    stackplot(FILLEDTime, [FILLED(:,[5,6,15,30]) MassDensitySpline'],{'B_z (nT)','V_{SW} (km/s)','D_{st} (nT)','F_{10.7} (s.f.u.)','\rho_{eq} (amu/cm^3)'},satnum,[],[],cuttime)
end
    

if(MakePaperPlots && stormcase==1)
    stackplot(FILLEDTime, [FILLED(:,[5,6,15,30]) MassDensitySpline'],{'B_z (nT)','V_{SW} (km/s)','D_{st} (nT)','F_{10.7} (s.f.u.)','\rho_{eq} (amu/cm^3)'},satnum,[3 DSTCut; 5 20],yranges(yr,:,:))
    
    if(satnum==6 && sy<=1989 && ey>=1989) %Make March 1989 Geomagnetic storm plot (if necessary time range still exists)
        stackplot(FILLEDTime, [FILLED(:,[5,6,15,30]) MassDensitySpline'],{'B_z (nT)','V_{SW} (km/s)','D_{st} (nT)','F_{10.7} (s.f.u.)','\rho_{eq} (amu/cm^3)'},satnum,[],[],[datenum('Mar-10-1989') datenum('Mar-18-1989')])
    end

    h=figure('Visible',visible);
    plot(-timewidth:1:timewidth*2,AVMDMat,'.')
    hold on;
    xa=(-timewidth:LongTimeScale:timewidth*2);
    plot(xa,AVMDs(1,:),'r','LineWidth',3);
    hold on; plot(xa,AVMDMatBars(:,:),'r-.','LineWidth',2);
    print -depsc2 -r200 figures/allstorms.eps
    print -dpng -r200 figures/PNGs/allstorms.png    
    if(strcmp(visible,'off')),close(h);end;
    
end

if(MakePaperPlots && stormcase==26) %Kp Stack plot
    stackplot(FILLEDTime, [FILLED(:,[5,6,13,30]) MassDensitySpline'],{'B_z (nT)','V_{SW} (km/s)','Kp','F_{10.7} (s.f.u.)','\rho_{eq} (amu/cm^3)'},satnum,[3 6],yranges(yr,:,:))
end


if(MakePaperPlots && (stormcase==2 || stormcase==24 || stormcase==1)) 
    tws=[20:25; 25:30];
    for varnum=[5 6 8 13 15 30]
    %varnum=30; %5 is Bz, 15 is dst, 8 is p, 6 is Vsw, 13 is kp, 30 is f10.7
    varname=headers{varnum};
    varunit=units{varnum};
    for i=1:2
        tw=tws(i,:);

        topcut=nanmedian(nanmedian(AVMat(:,tw,varnum),2));
        bottomcut=nanmedian(nanmedian(AVMat(:,tw,varnum),2));
        h=figure('Visible',visible);   
        top=nanmedian(AVMDMat(nanmedian(AVMat(:,tw,varnum),2)>=topcut,:));
        bottom=nanmedian(AVMDMat(nanmedian(AVMat(:,tw,varnum),2)<bottomcut,:));
        topbar=nanstd(AVMDMat(nanmedian(AVMat(:,tw,varnum),2)>=topcut,:))./sqrt(sum(~isnan(AVMDMat(nanmedian(AVMat(:,tw,varnum),2)>=topcut,:))));
        bottombar=nanstd(AVMDMat(nanmedian(AVMat(:,tw,varnum),2)<bottomcut,:))./sqrt(sum(~isnan(AVMDMat(nanmedian(AVMat(:,tw,varnum),2)<bottomcut,:))));
        tvals=ones(1,length(top));
        for j=1:length(top)
            [p,tvals(j)]=ranksum(AVMDMat(nanmedian(AVMat(:,tw,varnum),2)>=topcut,j),AVMDMat(nanmedian(AVMat(:,tw,varnum),2)<topcut,j));
        end
        %tvals=ttest2(AVMDMat(nanmedian(AVMat(:,tw,varnum),2)>=topcut,:),AVMDMat(nanmedian(AVMat(:,tw,varnum),2)<bottomcut,:));
        tvals(tvals==0)=NaN;
        tvalid=sum(~isnan(AVMDMat(nanmedian(AVMat(:,tw,varnum),2)>=topcut,:)));
        tvalid2=sum(~isnan(AVMDMat(nanmedian(AVMat(:,tw,varnum),2)<bottomcut,:)));
        %{
        for i=1:length(AVMDMat(1,:))
            top(i)=nanmean(AVMDMat(AVMat(:,i,varnum)>0,i));
            topbar(i)=nanstd(AVMDMat(AVMat(:,i,varnum)>0,i));
            bottom(i)=nanmean(AVMDMat(AVMat(:,i,varnum)<0,i));
            bottombar(i)=nanstd(AVMDMat(AVMat(:,i,varnum)<0,i));
        end
        %}
        
        if(nansum(tvals)==1) %That one plot with only 1 significant value
           h=figure('Visible',visible); [hhist, nhist]=hist(AVMDMat(nanmedian(AVMat(:,tw,varnum),2)>=topcut,tvals==1),1:2.5:120); hhhist=bar(nhist,hhist); set(hhhist,'FaceColor','none','EdgeColor','r'); 
           hold on;[hhist2, nhist2]=hist(AVMDMat(nanmedian(AVMat(:,tw,varnum),2)<topcut,tvals==1),1:2.5:120); hhhist2=bar(nhist2,hhist2); set(hhhist2,'FaceColor','none','EdgeColor','b'); 
           plot(nhist,hhist,'r','LineWidth',2); plot(nhist2,hhist2,'b','LineWidth',2); 
           legend(sprintf('Events with higher median %s',varname),sprintf('Events with lower median %s',varname))
           xlabel('\rho_{eq} (amu/cm^3)')
           ylabel('Number of events');
           print('-depsc2',sprintf('paper/figures/RhoBinned/RhoBinned%s-case%d-t0%d-tf%d-GOES%d-histogram.eps',varname,stormcase,tw(1),tw(end),satnum));
           print('-dpng',sprintf('paper/figures/RhoBinned/PNGs/RhoBinned%s-case%d-t0%d-tf%d-GOES%d-histogram.png',varname,stormcase,tw(1),tw(end),satnum));
        end
        
        
        h=figure('Visible',visible);
        hold on;
        plot(xa,tvalid);
        plot(xa,tvalid2,'r');
        ylims=get(gca,'YLim');
        ylabel('Valid events')
        xlabel('Time since onset (hours)')
        plot(xa,tvals+ylims(1)-1,'.','MarkerSize',15,'Color',[0.3 0.8 0.3])
        print('-depsc2',sprintf('paper/figures/RhoBinned/RhoBinned%s-case%d-t0%d-tf%d-GOES%d-valid.eps',varname,stormcase,tw(1),tw(end),satnum));
        print('-dpng','-r200',sprintf('paper/figures/RhoBinned/PNGs/RhoBinned%s-case%d-t0%d-tf%d-GOES%d-valid.png',varname,stormcase,tw(1),tw(end),satnum));
        
        
        t1=normrnd(ones(20,100),0.1);
        t1(1:10,:)=normrnd(ones(10,100)+4,0.1);
        fprintf('T-test verification, should be insignificant, percent sig: %2.2f %% \n',sum(ttest2(t1(1:5,:),t1(6:10,:))));
        fprintf('T-test verification, should be significant, percent sig: %2.2f %% \n',sum(ttest2(t1(1:5,:),t1(11:15,:))));
        
        
        h=figure('Visible',visible);
        plot(xa,nanmedian(AVMDMat),'r','LineWidth',3)
        hold on;

        plot(xa,top,'b','LineWidth',2)    
        plot(xa,bottom,'k','LineWidth',2)%'Color',[0.3 0.8 0.3])
        plot(xa,[top+topbar; top-topbar],'b-.')
        plot(xa,[bottom+bottombar; bottom-bottombar],'k-.')%,'Color',[0.3 0.8 0.3])
        ylims=get(gca,'YLim');
        plot(xa,tvals+ylims(1)-1,'.','MarkerSize',15,'Color',[0.3 0.8 0.3])
        rect=rectangle('Position',[tw(1)-timewidth-1 ylims(1) tw(end)-tw(1) ylims(2)-ylims(1)]) ;
        set(rect,'FaceColor',[0.9 0.9 0.9])
        uistack(rect,'bottom')
        grid on;
        if(tw(1)==25)
            text(0.01,0.97,'(b)','Units','normalized','FontSize',14)
        else
            text(0.01,0.97,'(a)','Units','normalized','FontSize',14)
        end
        set(gca,'xtick',[-timewidth:timewidth/2:timewidth*2]./LongTimeScale)
        xlim([-timewidth timewidth*2]./LongTimeScale)
        lh=legend('All \rho_{eq} events ',sprintf('%s  \\geq %2.2f %s; %d events ',varname,topcut,varunit,sum(nanmedian(AVMat(:,tw,varnum),2)>=topcut)),sprintf('%s < %2.2f %s; %d events ',varname,bottomcut,varunit,sum(nanmedian(AVMat(:,tw,varnum),2)<bottomcut)));
        set(lh,'box','off');
        title(sprintf('\\rho_{eq} events; GOES %d; %d-%d',satnum,sy,ey));
        ylabel('\rho_{eq} (amu/cm^3)')
        xlabel('Time since onset (hours)')
        print('-depsc2',sprintf('paper/figures/RhoBinned/RhoBinned%s-case%d-t0%d-tf%d-GOES%d.eps',varname,stormcase,tw(1),tw(end),satnum));
        print('-dpng','-r200',sprintf('paper/figures/RhoBinned/PNGs/RhoBinned%s-case%d-t0%d-tf%d-GOES%d.png',varname,stormcase,tw(1),tw(end),satnum));
    end
    end
    
    
    
    h=figure('Visible',visible);
    hist(nanmean(AVMDMat(:,20:24),2),0:5:100)
    h1 = findobj(gca,'Type','patch');
    set(h1,'FaceColor','r','EdgeColor','w');%,'facealpha',0.75
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
    print('-depsc2',sprintf('figures/rhobeforeafter-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('figures/PNGs/rhobeforeafter-GOES%d.png',satnum));    
    if(strcmp(visible,'off')),close(h);end;
    
    
    
    h=figure('Visible',visible);
    for i=1:100
        before(i)=nanmean(randsample(reshape(AVMDMat(:,21:25),[],1),4));
        after(i)=nanmean(randsample(reshape(AVMDMat(:,26:29),[],1),4));
    end
    [n,x]=hist(before,0:2:50);
    stairs(x,n,'LineWidth',2)
    hold on;
    h1 = findobj(gca,'Type','patch');
    set(h1,'FaceColor','r','EdgeColor','w');
    [n,x]=hist(after,0:2:50);
    stairs(x,n,'Color','r');
    legend('Before and at Onset','After Onset')
    %xlim([-5,105])
    axis tight;
    xlabel('Average \rho_{eq} (amu/cm^3)')
    ylabel('Count')
    title('Average \rho_{eq} Four Hours Before and After Storm Onset')
    %Printing segfaults for some reason...
    print('-depsc2',sprintf('figures/rhobeforeafter-boot-GOES%d.eps',satnum));
    print('-dpng', '-r200', sprintf('figures/PNGs/rhobeforeafter-boot-GOES%d.png',satnum));    
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
    %plot Bz, sort Pdyn
    binplot(AVMat(:,:,5),FILLED(:,8),starti,timewidth,LongTimeScale,plotthresh,{'B_z';'P_{dyn}';stormtype},{'nT';'nPa';stormunits},[sy; ey],satnum,visible);
    %plot rhoeq, sort Pdyn
    binplot(AVMDMat(:,:),FILLED(:,8),starti,timewidth,LongTimeScale,plotthresh,{'\rho_{eq}';'P_{dyn}';stormtype},{'amu/cm^3';'nPa';stormunits},[sy; ey],satnum,visible);
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
    print('-depsc2','-r200', sprintf('figures/DstRhoThresh-%d-%d.eps',sy,ey));
    print('-dpng','-r200', sprintf('figures/PNGs/DstRhoThresh-%d-%d.png',sy,ey));
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
    print('-depsc2','-r200', sprintf('figures/RandRhoThresh-%d-%d.eps',sy,ey));
    print('-dpng','-r200', sprintf('figures/PNGs/RandRhoThresh-%d-%d.png',sy,ey));
    if(strcmp(visible,'off')),close(h);end;
end

%Make main stack plots
if(MakePaperPlots)
    xa=(-timewidth:LongTimeScale:timewidth*2)./LongTimeScale;
    h=figure('Visible',visible);
    orient tall;
    
    h1=subplot('position',subplotstack(5,1));plot(xa,AVs(:,5),'+-','LineWidth',2); text(0.01,0.85,'B_z (nT)','Units','normalized','FontSize',14); %ylabel('B_z (nT)'); %Bz
    hold on; plot(xa,AVMatBars(:,:,5),'r-.'); ylim(yranges(yr,1,:))
    h2=subplot('position',subplotstack(5,2));plot(xa,AVs(:,6),'+-','LineWidth',2); text(0.01,0.85,'V_{SW} (km/s)','Units','normalized','FontSize',14); %ylabel('V_{SW} (km/s)');%V_sw
    hold on; plot(xa,AVMatBars(:,:,6),'r-.'); ylim(yranges(yr,2,:))
    h3=subplot('position',subplotstack(5,3));plot(xa,AVs(:,15),'+-','LineWidth',2); text(0.01,0.85,'D_{st} (nT)','Units','normalized','FontSize',14); %ylabel('D_{st} (nT)'); %dst
    hold on; plot(xa,AVMatBars(:,:,15),'r-.'); ylim(yranges(yr,3,:))
    if(DSTCut<0), plot([xa(1) xa(end)],[DSTCut DSTCut],'k-.','LineWidth',2); hold off; end
    h4=subplot('position',subplotstack(5,4));plot(xa,AVs(:,30),'+-','LineWidth',2); text(0.01,0.85,'F_{10.7} (s.f.u)','Units','normalized','FontSize',14); %ylabel('F10.7 (s.f.u.)');%f107
    hold on; plot(xa,AVMatBars(:,:,30),'r-.'); ylim(yranges(yr,4,:))
    set(findobj('type','axes'),'xticklabel',{[]})
    xv=[xa(1) xa(end)];
    subplot('position',subplotstack(5,5)); 
    [AX,H5,H6]=plotyy(xa,AVMDs(1,:),xa,AVnnans,'plot','bar');
    hold on; plot(xa,AVMDMatBars(:,:),'r-.'); %ylim(AX(1),yranges(yr,5,:))
    if(MDCut>0), plot([xa(1) xa(end)],[MDCut MDCut],'k-.','LineWidth',2); end
    set(AX(1),'Xlim',xv); set(AX(1),'YColor','r'); set(AX(1),'Color','none'); set(AX(1),'Ylim',yranges(yr,5,:),'YTick',round(linspace(yranges(yr,5,1),yranges(yr,5,2),length(get(AX(2),'YTick'))))); 
    set(AX(2),'Xlim',xv); set(AX(2),'YColor',[0 0.5 0.5]); set(AX(2),'Color','w');
    set(H5,'LineWidth',2);   set(H5,'marker','+','color','red'); set(get(H6,'child'),'FaceColor',[0 0.5 0.5]); uistack(AX(1));  
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
    set(findall(gcf,'-property','FontSize'),'FontSize',14) %Make all fonts same size 
    print('-depsc2','-r200',figurename);
    print('-dpng','-r200',strrep(strrep(figurename,'eps','png'),'figures','figures/PNGs'));
    if(strcmp(visible,'off')),close(h);end;
end

if(MakePaperPlots && MDCut>0) %Add a stack plot for Pressure for mass events
    xa=(-timewidth:LongTimeScale:timewidth*2)./LongTimeScale;
    h=figure('Visible',visible);
    orient tall;
    
    h1=subplot('position',subplotstack(5,1));plot(xa,AVs(:,5),'+-','LineWidth',2); text(0.01,0.85,'B_z (nT)','Units','normalized','FontSize',14); %ylabel('B_z (nT)'); %Bz
    hold on; plot(xa,AVMatBars(:,:,5),'r-.'); ylim(yranges(yr,1,:))
    h2=subplot('position',subplotstack(5,2));plot(xa,AVs(:,6),'+-','LineWidth',2); text(0.01,0.85,'V_{SW} (km/s)','Units','normalized','FontSize',14); %ylabel('V_{SW} (km/s)');%V_sw
    hold on; plot(xa,AVMatBars(:,:,6),'r-.'); ylim(yranges(yr,2,:))
    h3=subplot('position',subplotstack(5,3));plot(xa,AVs(:,15),'+-','LineWidth',2); text(0.01,0.85,'D_{st} (nT)','Units','normalized','FontSize',14); %ylabel('D_{st} (nT)'); %dst
    hold on; plot(xa,AVMatBars(:,:,15),'r-.'); ylim(yranges(yr,3,:))
    if(DSTCut<0), plot([xa(1) xa(end)],[DSTCut DSTCut],'k-.','LineWidth',2); hold off; end
    h4=subplot('position',subplotstack(5,4));plot(xa,AVs(:,8),'+-','LineWidth',2); text(0.01,0.85,'P_{dyn} (nPa)','Units','normalized','FontSize',14); %ylabel('F10.7 (s.f.u.)');%f107
    hold on; plot(xa,AVMatBars(:,:,8),'r-.'); %ylim(yranges(yr,4,:))
    set(findobj('type','axes'),'xticklabel',{[]})
    xv=[xa(1) xa(end)];
    subplot('position',subplotstack(5,5)); 
    [AX,H5,H6]=plotyy(xa,AVMDs(1,:),xa,AVnnans,'plot','bar');
    hold on; plot(xa,AVMDMatBars(:,:),'r-.'); %ylim(AX(1),yranges(yr,5,:))
    if(MDCut>0), plot([xa(1) xa(end)],[MDCut MDCut],'k-.','LineWidth',2); end
    set(AX(1),'Xlim',xv); set(AX(1),'YColor','r'); set(AX(1),'Color','none'); set(AX(1),'Ylim',yranges(yr,5,:),'YTick',round(linspace(yranges(yr,5,1),yranges(yr,5,2),length(get(AX(2),'YTick'))))); 
    set(AX(2),'Xlim',xv); set(AX(2),'YColor',[0 0.5 0.5]); set(AX(2),'Color','w');
    set(H5,'LineWidth',2);   set(H5,'marker','+','color','red'); set(get(H6,'child'),'FaceColor',[0 0.5 0.5]); uistack(AX(1));  
    text(0.01,0.85,'\rho_{eq} (amu/cm^3)','Units','normalized','FontSize',14); %ylabel(AX(1),'\rho_{eq} (amu/cm^3)'); 
    ylabel(AX(2),'# valid hourly values');
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[-timewidth:timewidth/2:timewidth*2]./LongTimeScale)
    linkaxes([AX h1 h2 h3 h4],'x')
    if(MDCut>0), title(h1,sprintf('%d \\rho_{eq} > %d amu/cm^3 events; %d-%d',length(starti),MDCut,sy,ey)); end
    if(LongTimeScale>1),xlabel('Time from start of event (day)');
    else xlabel('Time from start of event (hour)'); end
    fprintf('Average of %d storms %s (%d to %d)\n',length(duration),durationcaveat, sy,ey);
    set(findall(gcf,'-property','FontSize'),'FontSize',14) %Make all fonts same size
    print('-depsc2','-r200',strrep(figurename,'.eps','-withPressure.eps'));
    print('-dpng','-r200',strrep(strrep(figurename,'.eps','-withPressure.png'),'figures','figures/PNGs'));
    if(strcmp(visible,'off')),close(h);end;
    
    
end

if(stormcase==29) %NNBinaryOnset for DST
    %{
   F107Day=interptest(FILLEDTime,FILLED(:,30),FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)); 
   KpDay=interptest(FILLEDTime,FILLED(:,13),FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)); 
   VSWDay=interptest(FILLEDTime,FILLED(:,6),FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)); 
   DstDay=interptest(FILLEDTime,FILLED(:,15),FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)); 
   MDDay=interptest(FILLEDTime,MassDensitySpline',FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end));
    %}
    
    
    %Do it with hourly data first since that's already well defined for
    %events
    
    stormstarts=storms(1:end-1);
    stormstarts(stormstarts==-1)=0;
    %predict=NNBinaryOnset(FILLED(:,[6 13 15 30]),stormstarts');
    NNBinaryOnset([FILLED(:,[6 13 15 30]) circshift(FILLED(:,[6 13 15 30]),1) circshift(FILLED(:,[6 13 15 30]),2) circshift(FILLED(:,[6 13 15 30]),3)],stormstarts','dst-hourly');
    
   
    %Do daily medians
    stormstartsDay=stormstarts;
    while(mod(length(stormstartsDay),24)~=0)
        stormstartsDay(end+1)=0;
    end
    stormstartsDay=sum(reshape(stormstartsDay,24,[]));
    stormstartsDay(stormstartsDay>1)=1; %Any day with multiple storms is just marked as one event. Might be worth leaving off?
   
   F107Day=interptest(FILLEDTime,FILLED(:,30),FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)); 
   KpDay=interptest(FILLEDTime,FILLED(:,13),FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)); 
   VSWDay=interptest(FILLEDTime,FILLED(:,6),FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)); 
   DstDay=interptest(FILLEDTime,FILLED(:,15),FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)); 
   
   NNBinaryOnset([[F107Day KpDay VSWDay] circshift([F107Day KpDay VSWDay],1) circshift([F107Day KpDay VSWDay],2) circshift([F107Day KpDay VSWDay],3)],stormstartsDay','dst-daily');
   
    
   %Test random sort
   NNBinaryOnset([[F107Day KpDay VSWDay] circshift([F107Day KpDay VSWDay],1) circshift([F107Day KpDay VSWDay],2) circshift([F107Day KpDay VSWDay],3)],randsample(stormstartsDay,length(stormstartsDay))','dst-randomdaily');

end

if(stormcase==30) %NNBinary Onset for rhoeq, and only storm time
    stormstarts=storms(1:end-1);
    stormstarts(stormstarts==-1)=0;
    
    xi=find(stormstarts);
    NNBinaryOnset([FILLED([xi xi-1 xi-2 xi-3],[6 13 15 30]) MassDensitySpline([xi xi-1 xi-2 xi-3])'] ,stormstarts([xi xi-1 xi-2 xi-3])','hourly-withreq',{'F_{10.7}','Kp','D_{st}','V_{sw}','\rho_{eq}'});
    NNBinaryOnset(FILLED([xi xi-1 xi-2 xi-3],[6 13 15 30]) ,stormstarts([xi xi-1 xi-2 xi-3])','hourly',{'F_{10.7}','Kp','D_{st}','V_{sw}'});
  
    NNBinaryOnset([FILLED(:,[6 13 15 30]) MassDensitySpline' circshift([FILLED(:,[6 13 15 30]) MassDensitySpline'],1) circshift([FILLED(:,[6 13 15 30]) MassDensitySpline'],2) circshift([FILLED(:,[6 13 15 30]) MassDensitySpline'],3)],stormstarts','full-hourly-withreq');
    NNBinaryOnset([FILLED(:,[6 13 15 30]) circshift(FILLED(:,[6 13 15 30]),1) circshift(FILLED(:,[6 13 15 30]),2) circshift(FILLED(:,[6 13 15 30]),3)],stormstarts','full-hourly');
    
    stormstartsDay=stormstarts;
    while(mod(length(stormstartsDay),24)~=0)
        stormstartsDay(end+1)=0;
    end
    stormstartsDay=sum(reshape(stormstartsDay,24,[]));
    stormstartsDay(stormstartsDay>1)=1; %Any day with multiple storms is just marked as one event. Might be worth leaving off?
    
    xi=find(stormstartsDay);
   
   F107Day=interptest(FILLEDTime,FILLED(:,30),FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)); 
   KpDay=interptest(FILLEDTime,FILLED(:,13),FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)); 
   VSWDay=interptest(FILLEDTime,FILLED(:,6),FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)); 
   DstDay=interptest(FILLEDTime,FILLED(:,15),FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)); 
   MDDay=interptest(FILLEDTime,MassDensitySpline,FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end)); 
   
   NNBinaryOnset([F107Day([xi xi-1 xi-2 xi-3]) KpDay([xi xi-1 xi-2 xi-3]) DstDay([xi xi-1 xi-2 xi-3]) VSWDay([xi xi-1 xi-2 xi-3])] ,stormstartsDay([xi xi-1 xi-2 xi-3])','daily',{'F_{10.7}','Kp','D_{st}','V_{sw}'});
   NNBinaryOnset([F107Day([xi xi-1 xi-2 xi-3]) KpDay([xi xi-1 xi-2 xi-3]) VSWDay([xi xi-1 xi-2 xi-3]) MDDay([xi xi-1 xi-2 xi-3])] ,stormstartsDay([xi xi-1 xi-2 xi-3])','daily-withreq',{'F_{10.7}','Kp','V_{sw}','\rho_{eq}'});
   
   NNBinaryOnset([[F107Day KpDay VSWDay MDDay] circshift([F107Day KpDay VSWDay MDDay],1) circshift([F107Day KpDay VSWDay MDDay],2) circshift([F107Day KpDay VSWDay MDDay],3)],stormstartsDay','full-daily-withreq');
   NNBinaryOnset([[F107Day KpDay VSWDay] circshift([F107Day KpDay VSWDay],1) circshift([F107Day KpDay VSWDay],2) circshift([F107Day KpDay VSWDay],3)],stormstartsDay','full-daily');
   
   
   %Shuffled target
   shufflexi=randsample(xi,length(xi));
   NNBinaryOnset([F107Day([xi xi-1 xi-2 xi-3]) KpDay([xi xi-1 xi-2 xi-3]) VSWDay([xi xi-1 xi-2 xi-3])] ,stormstartsDay([shufflexi shufflexi-1 shufflexi-2 shufflexi-3])','randomdaily',{'F_{10.7}','Kp','V_{sw}'});
   
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
    print('-depsc2',sprintf('paper/figures/ccplot-GOES%d.eps',satnum)); 
    print('-dpng','-r200',sprintf('paper/figures/PNGs/ccplot-GOES%d.png',satnum)); 
    if(strcmp(visible,'off')),close(h);end;
    
end

if(MakePaperPlots && (stormcase==10 || stormcase==13)) %1-day takahashi. Do for any day plot
    %Bootstrap test for same daily means
    xa1d=(-timewidth:LongTimeScale:timewidth*2)./LongTimeScale;
    AV1d=interptest(-timewidth:1:timewidth*2,AVMDMat',-timewidth:LongTimeScale:timewidth*2)';
    AVm1=AV1d(:,4);
    AV0=AV1d(:,5);
    AV1=AV1d(:,6);
    
    dloops=1000;
    delta=zeros(4,dloops);
    for i=1:dloops
        rs1=randsample(AVm1,sum(~isnan(AVm1)),'true');
        rs2=randsample(AV0,sum(~isnan(AV0)),'true');
        
        delta(1,i)=(nanmedian(rs1)-nanmedian(rs2));
       
        rs1=randsample(AVm1,sum(~isnan(AVm1)),'true');
        rs2=randsample(AV1,sum(~isnan(AV1)),'true');
        delta(2,i)=(nanmedian(rs1)-nanmedian(rs2));
        
        rs1=randsample(AV1,sum(~isnan(AV1)),'true');
        rs2=randsample(AV0,sum(~isnan(AV0)),'true');
        delta(3,i)=(nanmedian(rs1)-nanmedian(rs2));
        
        rs1=randsample(AV0,sum(~isnan(AV0)),'true');
        rs2=randsample(AV0,sum(~isnan(AV0)),'true');
        delta(4,i)=(nanmedian(rs1)-nanmedian(rs2));
    end
    
    
    table=fopen(sprintf('tables/DeltaBootstraps-case%d.txt',stormcase),'w');
    Dm10=sum(delta(1,:)>=0)/length(delta(1,:));
    Dm11=sum(delta(2,:)<=0)/length(delta(2,:));
    D10=sum(delta(3,:)>=0)/length(delta(3,:));
    %[~,p1]=ttest2(AVm1,AV0);
    %[~,p2]=ttest2(AVm1,AV1);
    %[~,p3]=ttest2(AV1,AV0);
    
    fprintf(table,'\\begin{tabular}{|c|cc|}\n \\hline \n');
    fprintf(table,'Days & Difference (amu/cm$^3$) & \\%%  \\\\ \\hline\n');
    
    fprintf(table,'-1 0 & %2.2f & $%2.2f\\%%\\geq 0$  \\\\ \n',nanmedian(AVm1)-nanmedian(AV0),Dm10*100);
    fprintf(table,'-1 1 & %2.2f & $%2.2f\\%%\\leq 0$  \\\\ \n',nanmedian(AVm1)-nanmedian(AV1),Dm11*100);
    fprintf(table,' 1 0 & %2.2f & $%2.2f\\%%\\geq 0$  \\\\ \n',nanmedian(AV1)-nanmedian(AV0),D10*100);
    
    fprintf(table,'\\hline\n\\end{tabular}\n');
    
    fclose(table);
    
    
    g=figure('Visible',visible);
    h(1)=subplot('position',subplotstack(3,1));
    hist(delta(1,:),[-7.5:0.5:5.5])
    hold on;
    plot(nanmedian(AVm1)-nanmedian(AV0),0,'.','MarkerSize',20)
    plot([0 0],get(gca,'ylim'),'k-','LineWidth',1.5)
    text(0.01,0.9,'1 day pre-onset and onset','Units','normalized','FontSize',14)
    text(1,max(get(gca,'ylim'))*3/4,sprintf('%2.2f%% \\geq 0',Dm10*100))
    ylabel('Count')
    title('Bootstrapped differences of daily \rho_{eq} medians around event onset')
    
    h(2)=subplot('position',subplotstack(3,2));
    hist(delta(2,:),[-7.5:0.5:5.5])
    hold on;
    plot(nanmedian(AVm1)-nanmedian(AV1),0,'.','MarkerSize',20)
    plot([0 0],get(gca,'ylim'),'k-','LineWidth',1.5)
    text(0.01,0.9,'1 day pre-onset and 1 day post-onset','Units','normalized','FontSize',14)
    text(-2.5,max(get(gca,'ylim'))*3/4,sprintf('%2.2f%% \\leq 0',Dm11*100))
    ylabel('Count')
    
   h(3)=subplot('position',subplotstack(3,3));
    hist(delta(3,:),[-7:0.5:5])
    hold on;
    plot(nanmedian(AV1)-nanmedian(AV0),0,'.','MarkerSize',20)
    plot([0 0],get(gca,'ylim'),'k-','LineWidth',1.5)
    text(0.01,0.9,'1 day post-onset and onset','Units','normalized','FontSize',14)
    text(1,max(get(gca,'ylim'))*3/4,sprintf('%2.2f%% \\geq 0',D10*100))
    ylabel('Count')
    
    set(findobj('type','axes'),'xticklabel',{[]});
set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on')
axis tight;
%set(findobj('type','axes'),'xtick',get(h(end),'xtick'))
set(findobj('type','axes'),'xtick',[-7:1:5])
linkaxes(h,'x')
axis tight;
set(h(end),'xticklabel',[-7:1:5])
xlabel('Difference (amu/cm^3)','FontSize',14)
    
    print('-depsc2', '-r200', sprintf('figures/DailyBootstrapDifferences-GOES%d-case%d.eps',satnum,stormcase))
    print('-dpng', '-r200', sprintf('figures/PNGs/DailyBootstrapDifferences-GOES%d-case%d.png',satnum,stormcase))
    if(strcmp(visible,'off')),close(g);end;
    
    

    
end



if(MakePaperPlots && stormcase==25)
    fprintf('Starting 27 day plot of all satellites\n')
    
    load('data/OMNIdata.mat')

    F107(OMNITime>datenum('01-Jan-1992'))=[];
    OMNITime(OMNITime>datenum('01-Jan-1992'))=[];
    NewTime=OMNITime(1):24*27*(OMNITime(2)-OMNITime(1)):OMNITime(end);
    
    F10727Day=interptest(OMNITime,F107,NewTime);
    
    
    M2=load('data/DentonDensityAndTime_2.mat');
    M5=load('data/DentonDensityAndTime_5.mat');
    M6=load('data/DentonDensityAndTime_6.mat');
    M7=load('data/DentonDensityAndTime_7.mat');
    
    %Times 6 since they're 10 minute differences, times 24 and 27 to go to
    %hours and then 27-days
    NewMTime2=M2.DentonTime(1):24*27*6*(M2.DentonTime(2)-M2.DentonTime(1)):M2.DentonTime(end);
        NewMTime5=M5.DentonTime(1):24*27*6*(M5.DentonTime(2)-M5.DentonTime(1)):M5.DentonTime(end); 
        NewMTime6=M6.DentonTime(1):24*27*6*(M6.DentonTime(2)-M6.DentonTime(1)):M6.DentonTime(end);
        NewMTime7=M7.DentonTime(1):24*27*6*(M7.DentonTime(2)-M7.DentonTime(1)):M7.DentonTime(end);
    
    MD2=interptest(M2.DentonTime,M2.MassDensity,NewMTime2);
    MD5=interptest(M5.DentonTime,M5.MassDensity,NewMTime5);
    MD6=interptest(M6.DentonTime,M6.MassDensity,NewMTime6);
    MD7=interptest(M7.DentonTime,M7.MassDensity,NewMTime7);
    
    
    h=figure('Visible',visible,'Position',[0 0 1000 500]);
    set(gcf,'PaperPositionMode','auto')
    x=F10727Day;
    y=log10(MD6);
    [AX,H1,H2]=plotyy(NewTime,x,NewMTime6,y,'plot','plot');
    
    Colors=[[0 0.8 0.8]; [0.8 0 0.8]; [0 0 0]; [0.5 0.9 0.2]];
     
    set(H1,'marker','.','MarkerSize',15,'color','red'); set(AX(1),'YColor','r'); set(AX(2),'XTick',[]);
    set(H2,'LineStyle','-','color',Colors(3,:)); set(AX(2),'YColor','k');
    hold(AX(2));
    H3=plot(AX(2),NewMTime2,log10(MD2),'-','Color',Colors(1,:),'LineWidth',1.25);
    H4=plot(AX(2),NewMTime5,log10(MD5),'-','Color',Colors(2,:));
    H5=plot(AX(2),NewMTime7,log10(MD7),'-','Color',Colors(4,:));
    plot(AX(2),[NewMTime2(1) NewMTime2(end)],[1.5 1.5],'Color',Colors(1,:),'LineWidth',3)
    plot(AX(2),[NewMTime5(1) NewMTime5(end)],[1.52 1.52],'Color',Colors(2,:),'LineWidth',3)
    plot(AX(2),[NewMTime6(1) NewMTime6(end)],[1.54 1.54],'Color',Colors(3,:),'LineWidth',3)
    plot(AX(2),[NewMTime7(1) NewMTime7(end)],[1.56 1.56],'Color',Colors(4,:),'LineWidth',3)
    
    ylim(AX(1),[0,320])
    ylim(AX(2),[0.5,1.6])
    ylabel(AX(1),'F_{10.7} (s.f.u.)','FontSize',BigFont); ylabel(AX(2),'log_{10}[\rho_{eq} (amu/cm^2)]','FontSize',BigFont);
    [legh,objh]=legend([H1;H2;H3;H4;H5],'F_{10.7}',...
        sprintf('GOES 6 - %2.2f',min(min(corrcoef(interptest(M6.DentonTime,M6.MassDensity,NewTime),F10727Day,'rows','pairwise')))), ...
        sprintf('GOES 2 - %2.2f',min(min(corrcoef(interptest(M2.DentonTime,M2.MassDensity,NewTime),F10727Day,'rows','pairwise')))), ...
        sprintf('GOES 5 - %2.2f',min(min(corrcoef(interptest(M5.DentonTime,M5.MassDensity,NewTime),F10727Day,'rows','pairwise')))), ...
        sprintf('GOES 7 - %2.2f',min(min(corrcoef(interptest(M7.DentonTime,M7.MassDensity,NewTime),F10727Day,'rows','pairwise')))), ...
        'Location','SouthEast');
    set(objh,'LineWidth',2)
    set(AX(1),'YTick',0:50:300);
    %set(AX(2),'YTick',.5:0.25:1.5);
    xlabel('Year','FontSize',BigFont);
    title('27-day medians','FontSize',BigFont)
    linkaxes(AX,'x')
    datetick('keeplimits');
    grid on
    print('-depsc2', '-r200', 'figures/F107MD27d-all.eps')
    print('-dpng', '-r200', 'figures/PNGs/F107MD27d-all.png')
    if(strcmp(visible,'off')),close(h);end;
    
    
    
    g=figure('Visible',visible);
    
    scale=@(sc) (300*(sc-0.5))/(1.6-0.5);
    
    h(1)=subplot('position',subplotstack(4,1));
    %[AX,H1,H2]=plotyy(NewTime,x,NewMTime2,log10(MD2),'plot','plot');
    plot(NewTime,x,'r');
    hold on; plot(NewMTime2,scale(log10(MD2)),'b');
    legend('F_{10.7}',sprintf('GOES 2 - %2.2f',min(min(corrcoef(interptest(M2.DentonTime,M2.MassDensity,NewTime),F10727Day,'rows','pairwise')))),'Location','SouthEast');
    ylim([0,300])
    
    
    h(2)=subplot('position',subplotstack(4,2));
    %[AX,H1,H2]=plotyy(NewTime,x,NewMTime5,log10(MD5),'plot','plot');
    plot(NewTime,x,'r');
    hold on; plot(NewMTime5,scale(log10(MD5)),'b');
    legend('F_{10.7}',sprintf('GOES 5 - %2.2f',min(min(corrcoef(interptest(M5.DentonTime,M5.MassDensity,NewTime),F10727Day,'rows','pairwise')))),'Location','SouthEast');
    ylim([0,300])
    
    h(3)=subplot('position',subplotstack(4,3));
    %[AX,H1,H2]=plotyy(NewTime,x,NewMTime6,log10(MD6),'plot','plot');
    plot(NewTime,x,'r');
    hold on; plot(NewMTime6,scale(log10(MD6)),'b');
    legend('F_{10.7}',sprintf('GOES 6 - %2.2f',min(min(corrcoef(interptest(M6.DentonTime,M6.MassDensity,NewTime),F10727Day,'rows','pairwise')))),'Location','SouthEast');
    ylim([0,300])
    
    set(findobj('type','axes'),'xticklabel',{[]})
    
    h(4)=subplot('position',subplotstack(4,4));
    %[AX,H1,H2]=plotyy(NewTime,x,NewMTime7,log10(MD7),'plot','plot');
    plot(NewTime,x,'r');
    hold on; plot(NewMTime7,scale(log10(MD7)),'b');
    legend('F_{10.7}',sprintf('GOES 7 - %2.2f',min(min(corrcoef(interptest(M7.DentonTime,M7.MassDensity,NewTime),F10727Day,'rows','pairwise')))),'Location','SouthEast');
    ylim([0,300])
    
    ylabel('F_{10.7} (s.f.u.)','FontSize',BigFont); %ylabel(AX(2),'log_{10}[\rho_{eq} (amu/cm^2)]','FontSize',BigFont);
    
    set(gca,'YTick',0:50:300);
    xlabel('Year','FontSize',BigFont);
    
    
    
    linkaxes(h,'x')
    datetick('keeplimits');
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',get(h(4),'xtick'))
    print('-depsc2', '-r200', 'figures/F107MD27d-all2.eps')
    print('-dpng', '-r200', 'figures/PNGs/F107MD27d-all2.png')
    if(strcmp(visible,'off')),close(g);end;
end