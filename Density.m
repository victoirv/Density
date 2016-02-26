function Density(stormcase,satnum,cuttime)
if nargin < 1
    stormcase=1;
    satnum=6;
    cuttime=[0 datenum('Jan-01-2100')]; %Just to use all data. The years only serve to constrict, if possible
elseif nargin < 2
    satnum=6;
    cuttime=[0 datenum('Jan-01-2100')];
elseif nargin < 3
    cuttime=[0 datenum('Jan-01-2100')];
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

profile on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Script to read all Denton/OMNI data. In a separate script since multiple
%codes need it, and a script since it would be a lot of variables to return 
%from a function. Might turn it into three or four functions next
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DentonData


%-----------

yranges=zeros(5,4,2);
yranges(1,:,:)=[-2 2; 350 550; -60 0; 150 230];
yranges(2,:,:)=[-8 3; 350 600; -90 0; 160 210];
yranges(3,:,:)=[-8 3; 350 580; -90 0; 170 200];
yranges(4,:,:)=[-3 2; 400 550; -30 -10; 80 120];
yranges(5,:,:)=[-10 10; 0 1000; -100 0; 00 200]; %Made for overwriting with specific cases

LongTimeScale=1;%24*27; %How many hours to average over. Best stick to 1, 24, or 24*27
cutoffduration=1; %Minimum duration of events, in hours
cutconditions=0; %1=Look at only pre-noon. Change in FindStorms.m for other conditions
removef107=0; %Remove an IR predicted f10.7 dependence from mass density
DSTCut=0; %Set to 0 here, then if set by a case the code knows what variable to draw a threshold on
MDCut=0;
AECut=0;
nnanalysis=0; %Also set 0 here, changed on per-case basis.
ccanalysis=0;
figurename='paperfigures/stormavs-';
BigFont=16; %Set size for larger, readable fonts

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


%Remove F10.7 influence if requested
if(removef107)
    [~, ~, ~,xnew,~] = IR(MassDensitySpline,FILLED(:,30),0,1,0,0); %Remove F10.7 
    MassDensitySplineOriginal=MassDensitySpline;
    MassDensitySpline=MassDensitySpline-xnew;
end

timewidth=24; %How many hours before, and twice as many hours after, to select as the window for an event
if(LongTimeScale==24),timewidth=96; end %Make window larger for daily or monthly(ish) averages
if(LongTimeScale==24*27),timewidth=96*27; end
maxwidth=timewidth*2;

durationcaveat=''; %empty string if no cutoff set 
if(cutoffduration>1)
    durationcaveat=sprintf('longer than %d hours',cutoffduration);
end

%%%%%%
%Find storm indices and apply conditions to prune storms based on duration
%or overlapping edges of data
[starti,endi,duration]=FindStorms(storms,cutoffduration,cutconditions,maxwidth,MLTFit);

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%    PLOT SECTION    %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Default x-axis for plots
xa=(-timewidth:LongTimeScale:timewidth*2);

%Plots just made to prove data is accurate, or make some simple comparison
DataProofPlots

%Plots that are a simple scatter plot comparison
BasicPlots

%Plots meant to reproduce a previously published figure
PlotReproductions


%Every plot left here is more advanced analysis

profile viewer

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
    if(strcmp(visible,'off')),close(h);end;
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
if(strcmp(visible,'off')),close(h);end;
    
 
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



%Older plots, possibly for thesis, just batch comparing variables and
%making correlation tables
if(0)
    OldPlots
end