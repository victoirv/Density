function Density(stormcase,satnum,cutyears)
if nargin < 1
    stormcase=1;
    satnum=6;
    cutyears=[1900 2100]; %Just to use all data. The years only serve to constrict, if possible
elseif nargin < 2
    satnum=6;
    cutyears=[1900 2100];
elseif nargin < 3
    cutyears=[1900 2100];
end

if(any(stormcase==10:12))
    TakahashiCond=1;
else
    TakahashiCond=0;
end

MakePlots=0;
MakePaperPlots=1;
MakeBinPlots=0;
MakeDstThreshPlot=0;
MakeRandThreshPlot=0;
visible='off';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Script to read all Denton/OMNI data. In a separate script since multiple
%codes need it, and a script since it would be a lot of variables to return 
%from a function. Might turn it into three or four functions next
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DentonData


%-----------

yranges=zeros(4,4,2);
yranges(1,:,:)=[-2 2; 350 550; -60 0; 150 230];
yranges(2,:,:)=[-8 3; 350 600; -90 0; 160 210];
yranges(3,:,:)=[-8 3; 350 580; -90 0; 170 200];
yranges(4,:,:)=[-3 2; 400 550; -30 -10; 80 120];

LongTimeScale=1;%24*27; %How many hours to average over. Best stick to 1, 24, or 24*27
cutoffduration=1; %Minimum duration of events, in hours
cutconditions=0; %1=Look at only pre-noon
removef107=0;
DSTCut=0;
MDCut=0;
AECut=0;
nnanalysis=0;
ccanalysis=0;
figurename='paperfigures/stormavs-';
%Select kind of storm to look for

switch stormcase
    case 1
        storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
        DSTCut=-50;
        figurename=strcat(figurename,'dst.eps');
        yr=2;
        MakeBinPlots=1;
    case 2
        storms=diff([0 (MassDensityNanSpline>40)' 0]); %Mass Density Storm, started at 40
        MDCut=40;
        figurename=strcat(figurename,'mass.eps');
        yr=2;
        MakeBinPlots=1;
    case 3
        storms=[0 diff([0 (diff(MassDensityNanSpline)>10)' 0])];
        figurename=strcat(figurename,'diffden-10amu.eps');
        yr=3;
    case 4
        storms=[0 diff([0 (abs(MassDensityNanSpline(2:end)./MassDensityNanSpline(1:end-1))>1.3)' 0])];
        figurename=strcat(figurename,'diffden-30percent.eps');
        yr=3;
    case 5
        storms=diff([0 (FILLED(:,15)<-50)' 0]);
        DSTCut=-50;
        cutoffduration=12; %12 hour DST storm
        figurename=strcat(figurename,'dd12.eps');
        yr=2;
    case 6
        storms=diff([0 (MassDensityNanSpline>40)' 0]);
        MDCut=40;
        cutoffduration=12;
        figurename=strcat(figurename,'md12.eps');
        yr=2;
    case 7
         storms=diff([0 (FILLED(:,15)<-80)' 0]);
         DSTCut=-80;
         figurename=strcat(figurename,'d80.eps');
         yr=3;
    case 8
        storms=diff([0 (MassDensityNanSpline>70)' 0]);  
        MDCut=70;
        figurename=strcat(figurename,'m70.eps');
        yr=3;
    case 9
        storms=diff([0 (FILLED(:,15)<-50)' 0]);
        DSTCut=-50;
        removef107=1;
        figurename=strcat(figurename,'dst-nof107.eps');
        yr=3; %Dont know
    case 10 %Takahashi Fig 11
        storms=diff([0 (FILLED(:,15)<-50)' 0]);
        DSTCut=-50;
        cutconditions=1;
        LongTimeScale=24;
        figurename=strcat(figurename,'dst-50-tak.eps');
        yr=1;
    case 11 %Takahashi but Mass Storm
        storms=diff([0 (MassDensityNanSpline>40)' 0]); %Mass Density Storm, started at 40
        MDCut=40;
        cutconditions=1;
        LongTimeScale=24;
        figurename=strcat(figurename,'mass-tak.eps');
        yr=1;
    case 12 %Takahashi hourly dst
        storms=diff([0 (FILLED(:,15)<=-50)' 0]);
        DSTCut=-50;
        cutconditions=1;
        %cutoffduration=12;
        figurename=strcat(figurename,'dst-50-tak-hour.eps');
        yr=1;%Don't know
        %MakeBinPlots=1;
    case 13 %Full time range, daily medians
        storms=diff([0 (FILLED(:,15)<-50)' 0]);
        DSTCut=-50;
        LongTimeScale=24;
        figurename=strcat(figurename,'dst-day.eps');
        yr=1;
    case 14
        storms=diff([0 (MassDensityNanSpline>40)' 0]); 
        MDCut=40;
        LongTimeScale=24;
        figurename=strcat(figurename,'mass-day.eps');
        yr=1;
    case 15
        storms=diff([0 (FILLED(:,15)<-30)' 0]); %DST Storm
        DSTCut=-30; 
        figurename=strcat(figurename,'dst-30.eps');
        yr=2;
        MakeBinPlots=1;
    case 16 %AE Events
        Hrs=hour(FILLEDTime);
        for i=1:24
            avrhos(i)=nanmedian(AEFit(Hrs==(i-1)));
        end
        AEFit=AEFit-avrhos(Hrs+1)';
        AECut=300;
        storms=diff([0 (AEFit>AECut)' 0]);
        figurename=strcat(figurename,'AE.eps');
        yr=4;
        MakeBinPlots=1;
    case 17 %Random events
        starti = randsample(1:length(MassDensitySpline),700,false);
        storms=zeros(1,length(MassDensitySpline));
        storms(starti)=1; storms(starti+1)=-1;
        while(sum(storms)~=0)
            starti = randsample(1:length(MassDensitySpline),700,false);
            storms=zeros(1,length(MassDensitySpline));
            storms(starti)=1; storms(starti+1)=-1;
        end
        figurename=strcat(figurename,'random.eps');
        yr=4;
        MakeBinPlots=1;
    case 18    
    for i=1:24
        avrhos(i)=nanmedian(MassDensitySpline(FILLED(:,3)==(i-1)));
    end
    MassDensitySpline=MassDensitySpline-avrhos(FILLED(:,3)+1);
    storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
        DSTCut=-50;
        figurename=strcat(figurename,'dst-detrended.eps');
        yr=2;
        MakeBinPlots=1;
    case 19
        DSTCut=-80;
        storms=diff([0 (FILLED(:,15)<DSTCut)' 0]); %DST Storm
        figurename=strcat(figurename,'dst80.eps');
        yr=2;
        MakeBinPlots=1;
    case 20 %Specifically for a case later that loops over DST threshholds 
        DSTCut=-30; %This doesn't really matter
        storms=diff([0 (FILLED(:,15)<DSTCut)' 0]); %DST Storm
        figurename=strcat(figurename,'dst30.eps');
        yr=2;
        MakeDstThreshPlot=1;
    case 21 %Specifically for a case later that loops over DST threshholds 
        DSTCut=-30; %This doesn't really matter
        storms=diff([0 (FILLED(:,15)<DSTCut)' 0]); %DST Storm
        figurename=strcat(figurename,'dst30.eps');
        yr=2;
        MakeRandThreshPlot=1;
    case 22
        nnanalysis=1;
        storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
        DSTCut=-50;
        figurename=strcat(figurename,'dst.eps');
        yr=2;
    case 23
        ccanalysis=1;
        storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
        DSTCut=-50;
        figurename=strcat(figurename,'dst.eps');
        yr=2;
end
starti=find(storms>0);
endi=find(storms<0)-1;
duration=endi-starti+1;

%Shift event start points to next local DST minimum 
while(0) 
    ind=FILLED(starti+1,15)<FILLED(starti,15);
    if(sum(ind)==0)
        break
    end
    starti(ind)=starti(ind)+1;
end

durationcaveat=''; %empty string if no cutoff set 
if(cutoffduration>1)
    starti=starti(duration>cutoffduration);
    endi=endi(duration>cutoffduration);
    duration=duration(duration>cutoffduration);
    durationcaveat=sprintf('longer than %d hours',cutoffduration);
end

%Remove F10.7 influence
if(removef107)
    [~, ~, ~,xnew,~] = IR(MassDensitySpline,FILLED(:,30),0,1,0,0); %Remove F10.7 
    MassDensitySplineOriginal=MassDensitySpline;
    MassDensitySpline=MassDensitySpline-xnew;
end


timewidth=24; %How many hours before, and twice as many hours after, to select as the window for an event
if(LongTimeScale==24),timewidth=96; end
if(LongTimeScale==24*27),timewidth=96*27; end
maxwidth=timewidth*2;

%Remove storms where the window of interest extends past existing data
while(starti(1)-maxwidth<1)
    starti(1)=[]; endi(1)=[];
end
while(starti(end)+maxwidth>=length(MassDensitySpline))
    starti(end)=[]; endi(end)=[];
end

%Cut down found storms to only include pre-noon conditions
if(cutconditions) 
   endi(MLTFit(starti)>12 | MLTFit(starti)<6)=[];
   duration(MLTFit(starti)>12 | MLTFit(starti)<6)=[];
   starti(MLTFit(starti)>12 | MLTFit(starti)<6)=[];
end

%Build matrices storing all storms
stormi=1;
AVMat=[];
AVMDMat=[];
for i=1:length(starti)
    %datanum=sum(~isnan(MassDensitySpline(starti(i):endi(i))));
    AVMat(stormi,:,:)=FILLED((starti(i)-timewidth):starti(i)+timewidth*2,:); 
    AVMDMat(stormi,:)=MassDensitySpline((starti(i)-timewidth):starti(i)+timewidth*2);
    stormi=stormi+1;
end
AVs=nanmedian(AVMat,1);
AVMDs=nanmedian(AVMDMat);
AVnnans=sum(~isnan(AVMDMat));
AVMatBars=[AVs(1,:,:)-nanstd(AVMat(:,:,:))./sqrt(sum(~isnan(AVMat(:,:,:)))) ; AVs(1,:,:)+nanstd(AVMat(:,:,:))./sqrt(sum(~isnan(AVMat(:,:,:))))];
AVMDMatBars=[AVMDs(1,:)-nanstd(AVMDMat(:,:))./sqrt(sum(~isnan(AVMDMat(:,:)))) ; AVMDs(1,:)+nanstd(AVMDMat(:,:))./sqrt(sum(~isnan(AVMDMat(:,:))))];
AVs=squeeze(AVs);

if(LongTimeScale>1)
    AVs=interptest(-timewidth:1:timewidth*2,AVs,-timewidth:LongTimeScale:timewidth*2);
    AVMDs=interptest(-timewidth:1:timewidth*2,AVMDs',-timewidth:LongTimeScale:timewidth*2)';
    
    AVnnans=interptest(-timewidth:1:timewidth*2,AVnnans',-timewidth:LongTimeScale:timewidth*2)';
    NewAVMatBars(1,:,:)=interptest(-timewidth:1:timewidth*2,squeeze(AVMatBars(1,:,:)),-timewidth:LongTimeScale:timewidth*2); 
    NewAVMatBars(2,:,:)=interptest(-timewidth:1:timewidth*2,squeeze(AVMatBars(2,:,:)),-timewidth:LongTimeScale:timewidth*2);
    AVMatBars=NewAVMatBars;
    
    AVMDMatBars=[interptest(-timewidth:1:timewidth*2,AVMDMatBars(1,:)',-timewidth:LongTimeScale:timewidth*2)'; ...
        interptest(-timewidth:1:timewidth*2,AVMDMatBars(2,:)',-timewidth:LongTimeScale:timewidth*2)'];
    
    xa=(-timewidth:LongTimeScale:timewidth*2)./LongTimeScale;
    AVMDblock=nanmedian(reshape(AVMDMat(:,14:end),[],length(-timewidth:LongTimeScale:timewidth*2)-1));


    AVMD2=zeros(length(starti),13);
    AVMD2(:,1)=nanmedian(AVMDMat(:,1:12),2);
    for i=1:length(-timewidth:LongTimeScale:timewidth*2)-2
       AVMD2(:,i+1)=nanmedian(AVMDMat(:,i*LongTimeScale+1-12:i*LongTimeScale+12),2);
    end
    AVMD2(:,13)=nanmedian(AVMDMat(:,end-12:end),2);
    AVMD2=nanmedian(AVMD2);
end

%Shifting these to external functions for cleanliness and proof of independence/modularity
if(ccanalysis)
    CCAnalysis(AVMat,AVMDMat);
end

if(nnanalysis) 
    NNAnalysis(AVMat,AVMDMat);
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
end

%Median of medians vs block median
if(MakePaperPlots && LongTimeScale==24 && stormcase==10)
    h=figure('Visible',visible);
    plot(xa(2:end),AVMDblock)
    hold on; plot(xa(2:end),AVMDs(2:end),'r')
    plot(xa(2:end),AVMD2(2:end),'k')
    xlabel('Time from minimum of event (day)')
    ylabel('\rho_{eq} (amu/cm^3)')
    legend('Block median','MoM storm, day','MoM day, storm')
    print -depsc2 -r200 paperfigures/blockmedian.eps
end

%Tak2006 Fig 10 Kp vs M
if(MakePaperPlots && stormcase==1)
   figure('Visible',visible);
   subplot(311)
   x=starti;y=starti;
   for i=1:length(starti)
    x(i)=nanmedian(FILLED(starti(i)-3:starti(i),13) );
    y(i)=nanmean(FILLED(starti(i)-3:starti(i),32).*exp(-(3:-1:0)./3)');
   end
   plot(x,y,'.'); %KP vs M
   bins=unique(x); ymid=bins; ysd=bins;
   for i=1:length(bins)
       ymid(i)=nanmedian(y(x==bins(i)));
       ysd(i)=nanstd(y(x==bins(i)));
   end
   hold on; plot(bins, ymid,'r-.','MarkerSize',10); plot(bins, ymid-ysd,'k-.','MarkerSize',10); plot(bins, ymid+ysd,'k-.','MarkerSize',10);
   title(sprintf('%d events from GOES-%d, %d-%d',length(~isnan(y)),satnum,sy,ey))
   legend('3 hour means')
   
   subplot(312) %Go from 1 hour to 1 day
   for i=1:length(starti)
    x(i)=nanmedian(FILLED(starti(i)-24:starti(i),13));
    y(i)=nanmean(FILLED(starti(i)-24:starti(i),32).*exp(-(24:-1:0)./24)');
   end
   plot(x,y,'.');
   bins=unique(x); ymid=bins; ysd=bins;
   for i=1:length(bins)
       ymid(i)=nanmedian(y(x==bins(i)));
       ysd(i)=nanstd(y(x==bins(i)));
   end
   hold on; plot(bins, ymid,'r-.','MarkerSize',10); plot(bins, ymid-ysd,'k-.','MarkerSize',10); plot(bins, ymid+ysd,'k-.','MarkerSize',10);
   legend('1 day means')
   
   subplot(313) %1 hour to 3 days
   for i=1:length(starti)
    x(i)=nanmedian(FILLED(starti(i)-24*3:starti(i),13));
    y(i)=nanmean(FILLED(starti(i)-24*3:starti(i),32).*exp(-(24*3:-1:0)./(24*3))');
   end
   plot(x,y,'.');
   bins=unique(x); ymid=bins; ysd=bins;
   for i=1:length(bins)
       ymid(i)=nanmedian(y(x==bins(i)));
       ysd(i)=nanstd(y(x==bins(i)));
   end
   hold on; plot(bins, ymid,'r-.','MarkerSize',10); plot(bins, ymid-ysd,'k-.','MarkerSize',10); plot(bins, ymid+ysd,'k-.','MarkerSize',10);
   xlabel('Kp')
   ylabel('M (amu)')
   legend('3 day means')
   
   print -depsc2 -r200 paperfigures/KpvsM.eps
   
end


%Showing 'detrending' by removing F10.7 influencex`
if(MakePaperPlots && removef107)
    figure('Visible',visible);
    %plot(FILLEDTime,MassDensitySplineOriginal-(MassDensitySpline-nanmean(MassDensitySplineOriginal)),'-')
    plot(FILLEDTime,MassDensitySplineOriginal,'r.')
    hold on; plot(FILLEDTime,MassDensitySpline,'b.')
    legend('Original','F_{10.7} Removed');
    ylabel('\rho_{eq} (amu/cm^3)')
    xlabel('Year')
    datetick
    print -depsc2 -r200 paperfigures/f107removed.eps
end

if(MakePaperPlots && stormcase==1)
    figure('Visible',visible);
    xa=(-12:1:24);
    idx=arrayfun(@colon,starti(1:5)-12,starti(1:5)+24,'Uniform',false);
    plot(repmat(xa,5,1)',reshape(MassDensitySpline([idx{:}]),37,5),'r')
    hold on;
    idx2=arrayfun(@colon,starti(end-4:end)-12,starti(end-4:end)+24,'Uniform',false);
    plot(repmat(xa,5,1)',reshape(MassDensitySpline([idx2{:}]),37,5),'b')
    ylabel('\rho_{eq} (amu/cm^3)')
    xlabel('Hours from onset')
    title(sprintf('First 5 (red) and last 5 (blue) events from %d-%d',sy,ey));
    print -depsc2 -r200 paperfigures/TenEvents.eps
    print -dpng -r200 paperfigures/PNGs/TenEvents.png
end

if(MakePaperPlots && stormcase==1)
    
    stackplot(FILLEDTime, [FILLED(:,[5,6,15,30]) MassDensitySpline'],{'B_z (nT)','V_{SW} (km/s)','D_{st} (nT)','F_{10.7} (s.f.u.)','\rho_{eq} (amu/cm^3)'},satnum)
    
    if(satnum==6) %Make March 1989 Geomagnetic storm plot
        stackplot(FILLEDTime, [FILLED(:,[5,6,15,30]) MassDensitySpline'],{'B_z (nT)','V_{SW} (km/s)','D_{st} (nT)','F_{10.7} (s.f.u.)','\rho_{eq} (amu/cm^3)'},satnum,0,[datenum('Mar-10-1989') datenum('Mar-18-1989')])
    end
    %{
    h=figure('Visible',visible);
    orient tall;
    h1=subplot('position',subplotstack(5,1));plot(FILLEDTime,FILLED(:,5),'.'); text(0.01,0.9,'B_z (nT)','Units','normalized','FontSize',14); %Bz
    h2=subplot('position',subplotstack(5,2));plot(FILLEDTime,FILLED(:,6),'.'); text(0.01,0.9,'V_{SW} (km/s)','Units','normalized','FontSize',14); ylim([300,900]);%V_sw
    h3=subplot('position',subplotstack(5,3));plot(FILLEDTime,FILLED(:,15),'.');text(0.01,0.9,'D_{st} (nT)','Units','normalized','FontSize',14); ylim([-300,100]);
    hold on; plot([FILLEDTime(1) FILLEDTime(end)],[-40 -40],'r-.','LineWidth',4); hold off;
    h4=subplot('position',subplotstack(5,4));plot(FILLEDTime,FILLED(:,30),'.');text(0.01,0.9,'F_{10.7} (s.f.u.)','Units','normalized','FontSize',14); %f107
    h5=subplot('position',subplotstack(5,5));plot(FILLEDTime,MassDensitySpline,'r.');text(0.01,0.85,'\rho_{eq} (amu/cm^3)','Units','normalized','FontSize',14);%f107
    hold on; plot([FILLEDTime(1) FILLEDTime(end)],[40 40],'b-.','LineWidth',4);
    set(findobj('type','axes'),'xticklabel',{[]});
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on')
    axis tight;
    datetick('x','keeplimits')
    set(findobj('type','axes'),'xtick',get(h5,'xtick'))
    linkaxes([h5 h1 h2 h3 h4],'x')
    xlabel('Year')
    %set(gcf,'NextPlot','add'); axes; h = title(sprintf('All data',length(duration)));set(gca,'Visible','off');set(h,'Visible','on');
    print('-depsc2','-r200',sprintf('paperfigures/alldata-GOES%d-%d-%d.eps',satnum,sy,ey));
    print('-dpng','-r200',sprintf('paperfigures/PNGs/alldata-GOES%d-%d-%d.png',satnum,sy,ey));
    %}
    
    h=figure('Visible',visible);
    plot(-timewidth:1:timewidth*2,AVMDMat,'.')
    hold on;
    xa=(-timewidth:LongTimeScale:timewidth*2);
    plot(xa,AVMDs(1,:),'r','LineWidth',3);
    hold on; plot(xa,AVMDMatBars(:,:),'r-.','LineWidth',2);
    print -depsc2 -r200 paperfigures/allstorms.eps
    print -dpng -r200 paperfigures/PNGs/allstorms.png
    
    avrhos=zeros(1,24);
    for i=1:24
        avrhos(i)=nanmedian(MassDensity(round(MLT)==(i-1)));
    end
    figure('Visible',visible);
    plot(0:23,avrhos,'+-')
    ylabel('Median \rho_{eq} (amu/cm^3)')
    xlabel('MLT (hour)')
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[0:4:23])
    print -depsc2 -r200 paperfigures/rhoMLT.eps
    print -dpng -r200 paperfigures/PNGs/rhoMLT.png
    
    DHr=str2num(datestr(DentonTime,'HH'));
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
    if(DSTCut<0), plot([xa(1) xa(end)],[DSTCut DSTCut],'r-.','LineWidth',2); hold off; end
    h4=subplot('position',subplotstack(5,4));plot(xa,AVs(:,30),'+-','LineWidth',2); text(0.01,0.9,'F_{10.7} (s.f.u)','Units','normalized','FontSize',14); %ylabel('F10.7 (s.f.u.)');%f107
    hold on; plot(xa,AVMatBars(:,:,30),'r-.'); ylim(yranges(yr,4,:))
    set(findobj('type','axes'),'xticklabel',{[]})
    xv=[xa(1) xa(end)];
    subplot('position',subplotstack(5,5)); [AX,H5,H6]=plotyy(xa,AVMDs(1,:),xa,AVnnans,'plot','bar');
    hold on; plot(xa,AVMDMatBars(:,:),'r-.'); %ylim(AX(1),[0 90])
    if(MDCut>0), plot([xa(1) xa(end)],[MDCut MDCut],'r-.','LineWidth',2); end
    set(H5,'LineWidth',2);set(AX(2),'Xlim',xv); set(AX(1),'Xlim',xv);  set(H5,'marker','+','color','red'); set(AX(1),'YColor','r'); set(AX(2),'YColor',[0 0.5 0.5]); set(get(H6,'child'),'FaceColor',[0 0.5 0.5]); uistack(AX(1)); set(AX(1),'Color','none'); set(AX(2),'Color','w');
    text(0.01,0.85,'\rho_{eq} (amu/cm^3)','Units','normalized','FontSize',14); %ylabel(AX(1),'\rho_{eq} (amu/cm^3)'); 
    ylabel(AX(2),'# valid hourly values');
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[-timewidth:timewidth/2:timewidth*2]./LongTimeScale)
    linkaxes([AX h1 h2 h3 h4],'x')
    if(DSTCut<0), title(h1,sprintf('%d events of D_{st} < %d nT for %d-%d',length(starti),DSTCut,sy,ey)); end
    if(AECut<0), title(h1,sprintf('%d events of AE > %d for %d-%d',length(starti),AECut,sy,ey)); end
    if(MDCut>0), title(h1,sprintf('%d events of \\rho_{eq} > %damu/cm^3 for %d-%d',length(starti),MDCut,sy,ey)); end
    if(LongTimeScale>1),xlabel('Time from start of event (day)');
    else xlabel('Time from start of event (hour)'); end
    fprintf('Average of %d storms %s (%d to %d)\n',length(duration),durationcaveat, sy,ey);
    print('-depsc2','-r200',figurename);
    
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
end

%Nans per hour
if(MakePaperPlots && stormcase==1) 
    h=figure('Visible',visible);
    hist(FILLED(isnan(MassDensitySpline),3),0:23)
    axis([-1 24 0 3000])
    xlabel('UT Hour of event start')
    ylabel('Frequency')
    print -depsc2 -r200 paperfigures/nansbyhour.eps
    print -dpng -r200 paperfigures/PNGs/nansbyhour.png

    h=figure('Visible',visible);
    subplot(2,1,1)
    plot(xa,AVs(:,3))
    axis tight;
    ylabel('Hour average')
    grid on
    subplot(2,1,2)
    plot(xa,AVnnans,'r')
    axis tight;
    ylabel('Data available')
    xlabel('Time from event start (hr)')
    grid on
    print -depsc2 -r200 paperfigures/nansbyhour_storm.eps
    print -dpng -r200 paperfigures/PNGs/nansbyhour_storm.png
end

%DST vs rho_eq for 1 hour and 1 day
if(MakePaperPlots && stormcase==1) 
    F107Day=interptest(FILLEDTime,FILLED(:,30),FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end));
    MDDay=interptest(FILLEDTime,MassDensitySpline',FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end));
    F10727Day=interptest(FILLEDTime,FILLED(:,30),FILLEDTime(1):24*27*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end));
    MD27Day=interptest(FILLEDTime,MassDensitySpline',FILLEDTime(1):24*27*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end));
    
    h=figure('Visible',visible);
    
    plot(FILLED(~isnan(MassDensitySpline),30),log10(MassDensitySpline(~isnan(MassDensitySpline))),'.','MarkerSize',3);
    hold on;
    plot(F107Day,log10(MDDay),'r.','MarkerSize',8);    
    plot(F10727Day,log10(MD27Day),'go','LineWidth',1.5,'MarkerSize',12) ;
    xlabel('F_{10.7} (s.f.u.)','FontSize',16)
    ylabel('log(\rho_{eq}) (amu/cm^3)','FontSize',16)
    
    cc1=corrcoef(FILLED(~isnan(MassDensitySpline),30),log10(MassDensitySpline(~isnan(MassDensitySpline))),'rows','pairwise');
    cc1=cc1(1,2);
    cc2=corrcoef(F107Day,log10(MDDay),'rows','pairwise');
    cc2=cc2(1,2);
    cc3=corrcoef(F10727Day,log10(MD27Day),'rows','pairwise');
    cc3=cc3(1,2);
    
    legend(sprintf('One Hour cc=%2.2f',cc1),sprintf('One Day Median cc=%2.2f',cc2),sprintf('27-Day Median cc=%2.2f',cc3),'Location','SouthEast');
    legend boxoff;
    axis tight;
    print -depsc2 -r200 paperfigures/ccplot.eps 
    print -dpng -r200 paperfigures/PNGs/ccplot.png 
end


if(MakePaperPlots && stormcase==1)
    h=figure('Visible',visible);
    NewTime=FILLEDTime(1):24*27*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end);
    x=interptest(FILLEDTime,FILLED(:,30),NewTime);
    y=log10(interptest(FILLEDTime,MassDensitySpline',NewTime));
    cc=corrcoef(x,y,'rows','pairwise');
    [AX,H1,H2]=plotyy(NewTime,x,NewTime,y,'plot','plot');
    set(H1,'marker','.','color','red'); set(AX(1),'YColor','r'); set(AX(2),'XTick',[]);
    ylim(AX(2),[0.5,1.5])
    ylabel(AX(1),'F_{10.7\_27d}'); ylabel(AX(2),'log(\rho_{eq\_27d})');
    xlabel('Year','FontSize',16);
    title(sprintf('Data comparison for GOES %d - CC: %2.2f',satnum,cc(1,2)))
    linkaxes(AX,'x')
    datetick('keeplimits');
    grid on
    print('-depsc2', '-r200', sprintf('paperfigures/F107MD27d-GOES%d.eps',satnum))
    print('-dpng', '-r200', sprintf('paperfigures/PNGs/F107MD27d-GOES%d.png',satnum))
    
    h=figure('Visible',visible);
    NewTime=FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end);
    x=interptest(FILLEDTime,FILLED(:,30),NewTime);
    y=log10(interptest(FILLEDTime,MassDensitySpline',NewTime));
    cc=corrcoef(x,y,'rows','pairwise');
    [AX,H1,H2]=plotyy(NewTime,x,NewTime,y,'plot','plot');
    set(H1,'marker','.','color','red'); set(AX(1),'YColor','r'); set(AX(2),'XTick',[]);
    ylim(AX(2),[0.5,1.5])
    ylabel(AX(1),'F_{10.7\_1d}'); ylabel(AX(2),'log(\rho_{eq\_1d})');
    xlabel('Year','FontSize',16);
    title(sprintf('Data comparison for GOES %d - CC: %2.2f',satnum,cc(1,2)))
    linkaxes(AX,'x')
    datetick('keeplimits');
    grid on
    print('-depsc2', '-r200', sprintf('paperfigures/F107MD1d-GOES%d.eps',satnum))
    print('-dpng', '-r200', sprintf('paperfigures/PNGs/F107MD1d-GOES%d.png',satnum))
end

%Statistical comparisons
if(stormcase==1)
    %Show DST is significantly different during missing data times
    [p,h]=ranksum(FILLED(~isnan(MassDensitySpline),15),FILLED(isnan(MassDensitySpline),15),'Alpha',0.01)
    [h,p]=ttest2(FILLED(~isnan(MassDensitySpline),15),FILLED(isnan(MassDensitySpline),15),'Alpha',0.01)
    %Or when mass density gets above/below 40
    [p,h]=ranksum(FILLED((MassDensitySpline<40),15),FILLED((MassDensitySpline>40),15),'Alpha',0.01)
    [h,p]=ttest2(FILLED((MassDensitySpline<40),15),FILLED((MassDensitySpline>40),15),'Alpha',0.01)
    %Or pre-noon vs post-noon
    [p,h]=ranksum(FILLED(FILLED(:,3)>12,15),FILLED(FILLED(:,3)<12,15),'Alpha',0.01)
    [h,p]=ttest2(FILLED(FILLED(:,3)>12,15),FILLED(FILLED(:,3)<12,15),'Alpha',0.01)
end




%%%%%%%%%%%%%%%%%%%%%%
%Older plots, possibly for thesis, just batch comparing variables and
%making correlation tables

%{

VBS=1/2*FILLED(:,6).*(abs(FILLED(:,5))-FILLED(:,5));
VBz=FILLED(:,6).*FILLED(:,5);
BS=1/2*(abs(FILLED(:,5))-FILLED(:,5));
OMNIDensity=FILLED(:,7);

%Get re-scaled values for the sake of plotting overlays
OverlayFilled=FILLED;
OHr=FILLED(:,3);
DHr=hour(DentonTime);
for i=1:length(headers)
    OverlayFilled(:,i)=((250-150).*(OverlayFilled(:,i)-min(OverlayFilled(:,i))))./(max(OverlayFilled(:,i))-min(OverlayFilled(:,i))) + 150;
end
OverlayBS=((250-150).*(BS-min(BS)))./(max(BS)-min(BS)) + 150;
OverlayVBS=((250-150).*(VBS-min(VBS)))./(max(VBS)-min(VBS)) + 150;
OverlayVBz=((250-150).*(VBz-min(VBz)))./(max(VBz)-min(VBz)) + 150;
OverlayF107=((250-150).*(F107-min(F107)))./(max(F107)-min(F107)) + 150;
OverlayDBz=((250-150).*(DBz-min(DBz)))./(max(DBz)-min(DBz)) + 150;
OverlayDBS=((250-150).*(DBS-min(DBS)))./(max(DBS)-min(DBS)) + 150;


if(MakePlots) %Stack plots
    for i=1:length(headers)
        subplot('position',subplotstack(5,mod(i,4)+1));plot(AVs(:,i),'+-');
        ylabel(headers{i})
        if(mod(i,4)==0 || i==length(headers))
            subplot('position',subplotstack(5,5)); [AX,H1,H2]=plotyy(1:length(AVMDs),AVMDs(1,:),1:length(AVMDs),AVnnans,'plot','bar');
            set(AX(2),'Xlim',[1 length(AVMDs)]); set(AX(1),'Xlim',[1 length(AVMDs)]); set(get(H2,'child'),'facea',.3); set(H1,'marker','+','color','red'); set(AX(1),'YColor','r');
            %drawnow; hold on; xlim manual; h=bar(AVnnans); ch = get(h,'child'); set(ch,'facea',.3)
            ylabel(AX(1),'\rho_{eq}'); ylabel(AX(2),'# of values');
            %xlabel('Time before storm onset (hr)')
            print('-dpng',sprintf('figures/stormavs-%d',i));
            close all
        end
    end
end

if(MakePlots)
    NanPH=zeros(1,24);
    for i=0:23
        NanPH(i+1)=sum(isnan(MassDensitySpline(hour(FILLEDTime)==i)));
    end
    plot(0:23,NanPH);
    ylabel('Nans per hour')
    xlabel('Hour (UT)')
    print -dpng figures/nanph.png
end

if(MakePlots) %Make plot showing that most of F10.7 and Mass Density correlation is from long term effects
    [~, cb, ~,xnew,corr] = IR(log(MassDensity),log(F107),0,12,0,0);
    h=figure('Visible',visible);
    plot(DentonTime,log(MassDensity));
    hold on;
    plot(DentonTime,xnew,'g')
    xlabel('Time')
    ylabel('log(Mass Density)')
    legend('Data','Predicted')
    title(sprintf('Predicting Mass Density from F10.7 - cc(12)=%2.2f',corr))
    print -dpng figures/longtermf107.png
end

if(MakePlots)
    plot(OverlayBS)
    hold on; plot(OverlayVBz,'r')
    plot(MassDensitySpline,'g')
    plot(OverlayFilled(:,15),'k')
    legend('BS','VBz','MassDensity','DST','F10.7','Location','SouthWest')
    timestr=datestr(startt)
    title(timestr)
    print -dpng figures/onestorm-4.png
end


Na=0;
lag=0;
advance=0;
slice=0;
slicemin=10;
slicemax=24;

%x=detrend(MassDensitySpline);
x=MassDensitySpline;

%{
%Figure out what's going on with the year
[ca, cb, cc,xnew,corr] = IR(x,FILLED(:,1),Na,Nb,lag,advance);
figure; plot(xnew,MassDensitySpline,'+')

xm=mean(x);
ym=mean(FILLED(:,1));
sx=std(x);
sy=std(FILLED(:,1));
zx=(x-xm)./sx;
zy=(FILLED(:,1)-ym)./sy;
c=sum((zx.*zy))/(length(zx)-1);
title(sprintf('Denton density vs Year predicted density, my c: %2.2f - corrcoef c: %2.2f',c,corr));
print -dpng figures/yearmystery.png

lx=length(x);
lx5=floor(length(x)/5);
test=ones(1,lx);
test(lx5:2*lx5)=2;
test(2*lx5:3*lx5)=3;
test(3*lx5:4*lx5)=4;
test(4*lx5:5*lx5)=5;

test2=1:lx;

corrcoef(test,x)
%}

%Calculate correlation coefficients
cc=corrcoef(F107,DVsw);
fprintf('F10.7 Vsw: %2.2f \n',cc(2))
cc=corrcoef(F107,DBz);
fprintf('F10.7 Bz: %2.2f \n',cc(2))
cc=corrcoef(x,FILLED(:,15),'rows','complete');
fprintf('DST MassDensity: %2.2f \n',cc(2))

if(MakePlots)
    for i=0:23
        avf107(i+1)=mean(F107(DHr==i));
        avdens(i+1)=mean(MassDensity(DHr==i));
        avvsw(i+1)=mean(DVsw(DHr==i));
        avbz(i+1)=mean(DBz(DHr==i));
    end
    avf107=(avf107-min(avf107))/(max(avf107)-min(avf107));
    avdens=(avdens-min(avdens))/(max(avdens)-min(avdens));
    avvsw=(avvsw-min(avvsw))/(max(avvsw)-min(avvsw));
    avbz=(avbz-min(avbz))/(max(avbz)-min(avbz));
    h=figure('Visible',visible);
    plot(0:23,avf107,'r+-')
    hold on
    plot(0:23,avdens,'b+-')
    plot(0:23,avvsw,'k+-')
    plot(0:23,avbz,'m+-')
    legend('F107','Density','Vsw','Bz')
    xlabel('UT Hour')
    ylabel('Normalized Average')
    print '-dpng' 'figures/avf107.png'
end

%Add extra variables
headers=[headers;{'VBS';'BS'; 'VBz'}];
FILLED=[FILLED, VBS, BS, VBz];
OverlayFilled=[OverlayFilled, OverlayVBS, OverlayBS, OverlayVBz];

dheaders={'DBS';'F107';'lnF107';'DBz'};
DFILLED=[DBS, F107, log(F107), DBz];
DOverlayFilled=[OverlayDBS, OverlayF107, log(OverlayF107), OverlayDBz];

if(MakePlots && ~exist(sprintf('figures/OMNI_%s.png',headers{1}),'file'))
    for i=1:length(headers)
        close all;
        h=figure('Visible',visible);plot(FILLEDTime,FILLED(:,i))
        datetick
        ylabel(headers{i})
        xlabel('Time')
        title(sprintf('OMNI %s',headers{i}))
        filestring=sprintf('figures/OMNI_%s.png',headers{i});
        print('-dpng',filestring)
    end
end
if(MakePlots && ~exist(sprintf('figures/Denton_%s.png',dheaders{1}),'file'))
    for i=1:length(dheaders)
        close all;
        h=figure('Visible',visible);plot(DentonTime,DFILLED(:,i))
        datetick
        ylabel(dheaders{i})
        xlabel('Time')
        title(sprintf('Denton %s',dheaders{i}))
        filestring=sprintf('figures/Denton_%s.png',dheaders{i});
        print('-dpng',filestring)
    end
end
if(MakePlots && ~exist(sprintf('figures/Density_vs_%s.png',dheaders{1}),'file'))
    for i=1:length(headers)
        close all;
        h=figure('Visible',visible);plot(MassDensitySpline,FILLED(:,i),'+')
        ylabel(headers{i})
        xlabel('Mass Density')
        cc=corrcoef(MassDensitySpline,FILLED(:,i),'rows','complete');
        title(sprintf('%s - cc:%2.2f',headers{i},cc(2)))
        filestring=sprintf('figures/Density_vs_%s.png',headers{i});
        print('-dpng',filestring)
    end
end

%x=detrend(MassDensitySpline);

%Open the table for writing
if(slice)
    table=fopen('tableslice.txt','w');
else
    table=fopen('table.txt','w');
end
fprintf(table,'<pre>\n');
Nb1=1; %Just in case?
Nb=12;
x=log(MassDensitySpline);
fprintf(table,'Input \t CC1 \t CC12 \t PE(1) \t PE(12)\n');
BigTable={};
for i=1:length(headers)
    f=FILLED(:,i);
    if(slice)
        [~,~,~,xnew1,corr1] = IRslice(x,f,Na,Nb1,lag,advance,OHr,slicemin,slicemax);
        [~,cb,~,xnew,corr] = IRslice(x,f,Na,Nb,lag,advance,OHr,slicemin,slicemax);
    else
        [~,~,~,xnew1,corr1] = IR(x,f,Na,Nb1,lag,advance);
        [~, cb, ~,xnew,corr] = IR(x,f,Na,Nb,lag,advance);
    end
    pe1=pe_nonflag(x,xnew1);
    if(MakePlots), densityplot(FILLEDTime,[x,xnew1,OverlayFilled(:,i)],headers{i},Na,Nb1,corr1,pe1,visible), end
    pe=pe_nonflag(x,xnew);
    BigTable=[BigTable;{headers{i},corr1,corr,pe1,pe}];
    if(MakePlots), densityplot(FILLEDTime,[x,xnew,OverlayFilled(:,i)],headers{i},Na,Nb,corr,pe,visible), end
    if(MakePlots), densitycoefplot(-advance:Nb-advance-1,flipud(cb),headers{i},Na,Nb,corr,pe,visible), end
end

x=log(MassDensity);
for i=1:length(dheaders)
    f=DFILLED(:,i);
    if(slice)
        [~,~,~,xnew1,corr1] = IRslice(x,f,Na,Nb1,lag,advance,DHr,slicemin,slicemax);
        [~,cb,~,xnew,corr] = IRslice(x,f,Na,Nb,lag,advance,DHr,slicemin,slicemax);
    else
        [~,~,~,xnew1,corr1] = IR(x,f,Na,Nb1,lag,advance);
        [~, cb, ~,xnew,corr] = IR(x,f,Na,Nb,lag,advance);
    end
    pe1=pe_nonflag(x,xnew1);
    if(MakePlots), densityplot(DentonTime,[x,xnew1,DOverlayFilled(:,i)],dheaders{i},Na,Nb1,corr1,pe1,visible), end
    pe=pe_nonflag(x,xnew);
    BigTable=[BigTable;{dheaders{i},corr1,corr,pe1,pe}];
    if(MakePlots), densityplot(DentonTime,[x,xnew,DOverlayFilled(:,i)],headers{i},Na,Nb,corr,pe,visible), end
    if(MakePlots), densitycoefplot(-advance:Nb-advance-1,flipud(cb),headers{i},Na,Nb,corr,pe,visible), end
end

%Now for doubles
x=log(MassDensitySpline);
f=[FILLED(:,3) FILLED(:,5)]; %Hr and Bz
if(slice)
    [~,~,~,xnew1,corr1] = IRslice(x,f,Na,Nb1,lag,advance,OHr,slicemin,slicemax);
    [~,~,~,xnew,corr] = IRslice(x,f,Na,Nb,lag,advance,OHr,slicemin,slicemax);
else
    [~,~,~,xnew1,corr1] = IR(x,f,Na,Nb1,lag,advance);
    [~, ~, ~,xnew,corr] = IR(x,f,Na,Nb,lag,advance);
end
pe1=pe_nonflag(x,xnew1);
pe=pe_nonflag(x,xnew);
BigTable=[BigTable;{'Hr+Bz',corr1,corr,pe1,pe}];

f=[FILLED(:,5) FILLED(:,6)]; %Bz and V
if(slice)
    [~,~,~,xnew1,corr1] = IRslice(x,f,Na,Nb1,lag,advance,OHr,slicemin,slicemax);
    [~,~,~,xnew,corr] = IRslice(x,f,Na,Nb,lag,advance,OHr,slicemin,slicemax);
else
    [~,~,~,xnew1,corr1] = IR(x,f,Na,Nb1,lag,advance);
    [~, ~, ~,xnew,corr] = IR(x,f,Na,Nb,lag,advance);
end
pe1=pe_nonflag(x,xnew1);
pe=pe_nonflag(x,xnew);
BigTable=[BigTable;{'Bz+V',corr1,corr,pe1,pe}];

x=log(MassDensity);
f=[F107 DBz]; %F107 and Bz
if(slice)
    [~,~,~,xnew1,corr1] = IRslice(x,f,Na,Nb1,lag,advance,DHr,slicemin,slicemax);
    [~,~,~,xnew,corr] = IRslice(x,f,Na,Nb,lag,advance,DHr,slicemin,slicemax);
else
    [~,~,~,xnew1,corr1] = IR(x,f,Na,Nb1,lag,advance);
    [~, ~, ~,xnew,corr] = IR(x,f,Na,Nb,lag,advance);
end
pe1=pe_nonflag(x,xnew1);
pe=pe_nonflag(x,xnew);
BigTable=[BigTable;{'F107+Bz',corr1,corr,pe1,pe}];
f=[F107 DVsw]; %Bz and V
if(slice)
    [~,~,~,xnew1,corr1] = IRslice(x,f,Na,Nb1,lag,advance,DHr,slicemin,slicemax);
    [~,~,~,xnew,corr] = IRslice(x,f,Na,Nb,lag,advance,DHr,slicemin,slicemax);
else
    [~,~,~,xnew1,corr1] = IR(x,f,Na,Nb1,lag,advance);
    [~, ~, ~,xnew,corr] = IR(x,f,Na,Nb,lag,advance);
end
pe1=pe_nonflag(x,xnew1);
pe=pe_nonflag(x,xnew);
BigTable=[BigTable;{'F107+V',corr1,corr,pe1,pe}];


BigTable=sortrows(BigTable,2);
for i=1:size(BigTable,1)
    if(BigTable{i,5}>-10)
        fprintf(table,'%s\t %2.2f \t %2.2f \t %2.2f \t %2.2f\n',BigTable{i,1},BigTable{i,2},BigTable{i,3},BigTable{i,4},BigTable{i,5});
    else
        fprintf(table,'%s\t %2.2f \t %2.2f \t %2.2f \t %2.0e\n',BigTable{i,1},BigTable{i,2},BigTable{i,3},BigTable{i,4},BigTable{i,5});
    end
end
fprintf(table,'\n</pre>');

fclose(table);
system('cat README.txt table.txt > README.md');
%}

end

function densityplot(x,ys,string,Na,Nb,corr,eff,visible)
close all;
test=size(ys);
if test(1)>test(2)
    ys=ys';
end
h=figure('Visible',visible);plot(x,ys)
legend('Data','Prediction',string,'Location','NorthEast')
ylabel('Density')
xlabel('Time')
datetick
title(sprintf('Denton Density from OMNI %s: Nx:%d Nf:%d, corr: %2.3f, eff: %2.3f',string,Na,Nb,corr,eff))
filestring=sprintf('figures/density_%s_%d_%d.png',string,Na,Nb);
print(h,'-dpng',filestring)
end

function densitycoefplot(x,ys,string,Na,Nb,corr,eff,visible)
close all;
h=figure('Visible',visible);
plot(x,ys)
ylabel('Impulse Coefficient')
xlabel('Time Lags')
grid
title(sprintf('Denton Density from OMNI %s: Nx:%d Nf:%d, corr: %2.3f, eff: %2.3f',string,Na,Nb,corr,eff))
filestring=sprintf('figures/density_%s_%d_%d_Cb.png',string,Na,Nb);
print(h,'-dpng',filestring)
end
