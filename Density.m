function Density(stormcase,satnum,cuttime,skiptoplots)
if nargin < 1
    stormcase=1;
end
if nargin < 2
    satnum=6;
end
if nargin < 3
    cuttime=[0 datenum('Jan-01-2100')]; %Just to use all data. The years only serve to constrict, if possible
end
if nargin < 4
    skiptoplots=0;
end

savefilename=sprintf('data/DensitySaveState_%d_%d_%2.2f_%2.2f.mat',stormcase,satnum,cuttime(1),cuttime(2));
if(~(exist(savefilename,'file') && skiptoplots))

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

%profile on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Script to read all Denton/OMNI data. In a separate script since multiple
%codes need it, and a script since it would be a lot of variables to return 
%from a function. Might turn it into three or four functions next
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DentonData


%-----------

yranges=zeros(5,5,2);
yranges(1,:,:)=[-2 2; 350 550; -60 0; 150 230; 8 16];
yranges(2,:,:)=[-8 3; 400 600; -80 0; 140 180; 8 18];
yranges(3,:,:)=[-8 3; 350 580; -90 0; 170 200; 12 22];
yranges(4,:,:)=[-3 2; 400 550; -30 -10; 80 120; 12 22];
yranges(5,:,:)=[-10 10; 0 1000; -100 0; 00 200; 12 22]; %Made for overwriting with specific cases

LongTimeScale=1;%24*27; %How many hours to average over. Best stick to 1, 24, or 24*27
cutoffduration=1; %Minimum duration of events, in hours
cutconditions=0; %1=Look at only pre-noon. Change in FindStorms.m for other conditions
removef107=0; %Remove an IR predicted f10.7 dependence from mass density
DSTCut=0; %Set to 0 here, then if set by a case the code knows what variable to draw a threshold on
MDCut=0;
AECut=0;
nnanalysis=0; %Also set 0 here, changed on per-case basis.
ccanalysis=0;
arxanalysis=0;
figurename='paper/figures/stormavs-';
BigFont=16; %Set size for larger, readable fonts

%Select kind of storm to look for
switch stormcase
    case 1
        storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
        DSTCut=-50;
        figurename=strcat(figurename,'dst',sprintf('-GOES%d.eps',satnum));
        yr=2;
        cutconditions=1;
        MakeBinPlots=1;
    case 2
        storms=diff([0 (MassDensityNanSpline>40)' 0]); %Mass Density Storm, started at 40
        MDCut=40;
        figurename=strcat(figurename,'mass',sprintf('-GOES%d.eps',satnum));
        yr=5;
        MakeBinPlots=1;
        yranges(5,:,:)=[-1 3; 350 450; -20 0; 180 210; 10 60];
    case 3
        storms=[0 diff([0 (diff(MassDensityNanSpline)>5)' 0])];
        figurename=strcat(figurename,'diffden-5amu',sprintf('-GOES%d.eps',satnum));
        yr=5;
        yranges(5,:,:)=[-1 3; 350 450; -20 0; 180 210; 10 60];
    case 4
        storms=[0 diff([0 (abs(MassDensityNanSpline(2:end)./MassDensityNanSpline(1:end-1))>1.3)' 0])];
        figurename=strcat(figurename,'diffden-30percent',sprintf('-GOES%d.eps',satnum));
        yr=3;
    case 5
        storms=diff([0 (FILLED(:,15)<-50)' 0]);
        DSTCut=-50;
        cutoffduration=12; %12 hour DST storm
        figurename=strcat(figurename,'dd12',sprintf('-GOES%d.eps',satnum));
        yr=5;
        cutconditions=1;
        yranges(5,:,:)=[-10 5; 400 700; -100 0; 150 200; 0 70];
    case 6
        storms=diff([0 (MassDensityNanSpline>40)' 0]);
        MDCut=40;
        cutoffduration=12;
        figurename=strcat(figurename,'md12',sprintf('-GOES%d.eps',satnum));
        yr=2;
    case 7
         storms=diff([0 (FILLED(:,15)<-80)' 0]);
         DSTCut=-80;
         figurename=strcat(figurename,'d80',sprintf('-GOES%d.eps',satnum));
         yr=3;
    case 8
        storms=diff([0 (MassDensityNanSpline>70)' 0]);  
        MDCut=70;
        figurename=strcat(figurename,'m70.',sprintf('-GOES%d.eps',satnum));
        yr=3;
    case 9
        storms=diff([0 (FILLED(:,15)<-50)' 0]);
        DSTCut=-50;
        removef107=1;
        figurename=strcat(figurename,'dst-nof107',sprintf('-GOES%d.eps',satnum));
        yr=3; %Dont know
    case 10 %Takahashi Fig 11
        storms=diff([0 (FILLED(:,15)<-50)' 0]);
        DSTCut=-50;
        cutconditions=1;
        cutoffduration=12;
        LongTimeScale=24;
        figurename=strcat(figurename,'dst-50-tak',sprintf('-GOES%d.eps',satnum));
        yr=5;
        yranges(5,:,:)=[-2 2; 400 600; -80 0; 150 230; 10 30];
        
    case 11 %Takahashi but Mass Storm
        storms=diff([0 (MassDensityNanSpline>40)' 0]); %Mass Density Storm, started at 40
        MDCut=40;
        cutconditions=1;
        LongTimeScale=24;
        figurename=strcat(figurename,'mass-tak',sprintf('-GOES%d.eps',satnum));
        yr=1;
    case 12 %Takahashi hourly dst
        storms=diff([0 (FILLED(:,15)<=-50)' 0]);
        DSTCut=-50;
        cutconditions=1;
        %cutoffduration=12;
        figurename=strcat(figurename,'dst-50-tak-hour',sprintf('-GOES%d.eps',satnum));
        yr=1;%Don't know
        %MakeBinPlots=1;
    case 13 %Full time range, daily medians
        storms=diff([0 (FILLED(:,15)<-50)' 0]);
        DSTCut=-50;
        LongTimeScale=24;
        cutoffduration=12;
        cutconditions=1;
        figurename=strcat(figurename,'dst-day',sprintf('-GOES%d.eps',satnum));
        yr=5;
        yranges(5,:,:)=[-2 2; 400 600; -80 0; 150 230; 5 25];
    case 14
        storms=diff([0 (MassDensityNanSpline>40)' 0]); 
        MDCut=40;
        LongTimeScale=24;
        figurename=strcat(figurename,'mass-day',sprintf('-GOES%d.eps',satnum));
        yr=1;
    case 15
        storms=diff([0 (FILLED(:,15)<-30)' 0]); %DST Storm
        DSTCut=-30; 
        figurename=strcat(figurename,'dst-30',sprintf('-GOES%d.eps',satnum));
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
        figurename=strcat(figurename,'AE',sprintf('-GOES%d.eps',satnum));
        yr=4;
        MakeBinPlots=1;
    case 17 %Random events
        starti = randsample(1:floor(length(MassDensitySpline)/3),700,false).*3; %So each storm is at least 4 hours long
        storms=zeros(1,length(MassDensitySpline));
        storms(starti)=1; storms(starti+1)=-1;
        while(sum(storms)~=0)
            starti = randsample(1:length(MassDensitySpline),700,false);
            storms=zeros(1,length(MassDensitySpline));
            storms(starti)=1; storms(starti+1)=-1;
        end
        figurename=strcat(figurename,'random',sprintf('-GOES%d.eps',satnum));
        yr=5;
        yranges(5,:,:)=[-1 2; 400 500; -20 0; 100 120; 5 30];
        MakeBinPlots=1;
    case 18    
    for i=1:24
        avrhos(i)=nanmedian(MassDensitySpline(FILLED(:,3)==(i-1)));
    end
    MassDensitySpline=MassDensitySpline-avrhos(FILLED(:,3)+1);
    storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
        DSTCut=-50;
        figurename=strcat(figurename,'dst-detrended',sprintf('-GOES%d.eps',satnum));
        yr=2;
        MakeBinPlots=1;
    case 19
        DSTCut=-80;
        storms=diff([0 (FILLED(:,15)<DSTCut)' 0]); %DST Storm
        figurename=strcat(figurename,'dst80',sprintf('-GOES%d.eps',satnum));
        yr=2;
        MakeBinPlots=1;
    case 20 %Specifically for a case later that loops over DST threshholds 
        DSTCut=-30; %This doesn't really matter
        storms=diff([0 (FILLED(:,15)<DSTCut)' 0]); %DST Storm
        figurename=strcat(figurename,'dst30',sprintf('-GOES%d.eps',satnum));
        yr=2;
        MakeDstThreshPlot=1;
    case 21 %Specifically for a case later that loops over DST threshholds 
        DSTCut=-30; %This doesn't really matter
        storms=diff([0 (FILLED(:,15)<DSTCut)' 0]); %DST Storm
        figurename=strcat(figurename,'dst30',sprintf('-GOES%d.eps',satnum));
        yr=2;
        MakeRandThreshPlot=1;
    case 22
        nnanalysis=1;
        storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
        DSTCut=-50;
        figurename=strcat(figurename,'dst',sprintf('-GOES%d.eps',satnum));
        yr=2;
    case 23
        ccanalysis=1;
        storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
        DSTCut=-50;
        figurename=strcat(figurename,'dst',sprintf('-GOES%d.eps',satnum));
        yr=2;
    case 24
                storms=diff([0 (MassDensityNanSpline>20)' 0]); %Mass Density Storm, started at 20
        MDCut=20;
        figurename=strcat(figurename,'mass-gt20',sprintf('-GOES%d.eps',satnum));
        yr=5;
        yranges(5,:,:)=[-2 2; 350 550; -60 0; 150 230; 10 30];
    case 25 %Specifically for 27 day F10.7 vs mass density plot
        storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
        DSTCut=-50;
        figurename=strcat(figurename,'dst',sprintf('-GOES%d.eps',satnum));
        yr=2;
    case 26
        storms=diff([0 (FILLED(:,13)>6)' 0]); %KP Storm
        figurename=strcat(figurename,'Kp',sprintf('-GOES%d.eps',satnum));
        yr=5;
        yranges(5,:,:)=[-2 2; 350 550; 1 10; 150 230; 10 30];
    case 27 %For doing 3day/27day F107 plots
        storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
        DSTCut=-50;
        figurename=strcat(figurename,'ignore',sprintf('-GOES%d.eps',satnum));
        yr=2;
        MakeBinPlots=1;
    case 28 %For making specific stack plots of events (i.e. alldata-* plots with specific cuttimes)
        storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
        DSTCut=-50;
        figurename=strcat(figurename,'ignore',sprintf('-GOES%d.eps',satnum));
        yr=2;
    case 29
        storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
        DSTCut=-50;
        figurename=strcat(figurename,'ignore',sprintf('-GOES%d.eps',satnum));
        yr=2;
    case 30
        storms=diff([0 (MassDensityNanSpline>20)' 0]); %Density Storm
        MDCut=20;
        figurename=strcat(figurename,'ignore',sprintf('-GOES%d.eps',satnum));
        yr=2;
    case 31 %Random daily events
        starti = randsample(1:floor(length(MassDensitySpline)/3),700,false).*3; %So each storm is at least 4 hours long
        storms=zeros(1,length(MassDensitySpline));
        storms(starti)=1; storms(starti+1)=-1;
        while(sum(storms)~=0)
            starti = randsample(1:length(MassDensitySpline),700,false);
            storms=zeros(1,length(MassDensitySpline));
            storms(starti)=1; storms(starti+1)=-1;
        end
        LongTimeScale=24;
        figurename=strcat(figurename,'randomdaily',sprintf('-GOES%d.eps',satnum));
        yr=5;
        yranges(5,:,:)=[-1 2; 400 500; -20 0; 100 120; 5 30];
    case 32
        storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
        DSTCut=-50;
        figurename=strcat(figurename,'ignore',sprintf('-GOES%d.eps',satnum));
        yr=2;
        arxanalysis=1;
    case 33
        storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
        DSTCut=-50;
        figurename=strcat(figurename,'dst-normalized',sprintf('-GOES%d.eps',satnum));
        yr=5;
        yranges(5,:,:)=[-6 2; 350 550; -60 0; 150 230; -4 4];
    case 34
        storms=diff([0 (FILLED(:,15)<-50)' 0]);
        DSTCut=-50;
        cutoffduration=6; %12 hour DST storm
        figurename=strcat(figurename,'dd6',sprintf('-GOES%d.eps',satnum));
        yr=5;
        yranges(5,:,:)=[-10 5; 400 700; -100 0; 150 200; 0 30];
     case 35
        storms=diff([0 (FILLED(:,15)<-50)' 0]);
        DSTCut=-50;
        LongTimeScale=24;
        cutoffduration=6; %12 hour DST storm
        figurename=strcat(figurename,'dd6-day',sprintf('-GOES%d.eps',satnum));
        yr=5;
        yranges(5,:,:)=[-10 5; 400 700; -100 0; 150 200; 0 30];
end

eventtype='D_{st}';
if(MDCut>0)
    eventtype='\rho_{eq}';
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
[starti,endi,duration]=FindStorms(storms,cutoffduration,cutconditions,maxwidth,MLTFit,FILLED);

%Build matrices storing all storms
stormi=1;
AVMat=[];
AVMDMat=[];
for i=1:length(starti)
    %datanum=sum(~isnan(MassDensitySpline(starti(i):endi(i))));
    AVMat(stormi,:,:)=FILLED((starti(i)-timewidth):starti(i)+timewidth*2-1,:); 
    AVMDMat(stormi,:)=MassDensitySpline((starti(i)-timewidth):starti(i)+timewidth*2-1);
    stormi=stormi+1;
end
AVs=nanmedian(AVMat,1);
AVMDs=nanmedian(AVMDMat);
if(stormcase==33)
    NormAVMDMat=AVMDMat;
    for i=1:length(AVMat)
        NormAVMDMat(i,:)=(AVMDMat(i,:)-min(AVMDMat(i,:)))./(max(AVMDMat(i,:))-min(AVMDMat(i,:)));
    end
    AVMDs=nanmedian(NormAVMDMat);
end
AVnnans=sum(~isnan(AVMDMat));
AVMatBars=[AVs(1,:,:)-nanstd(AVMat(:,:,:))./sqrt(sum(~isnan(AVMat(:,:,:)))) ; AVs(1,:,:)+nanstd(AVMat(:,:,:))./sqrt(sum(~isnan(AVMat(:,:,:))))];
AVMDMatBars=[AVMDs(1,:)-nanstd(AVMDMat(:,:))./sqrt(sum(~isnan(AVMDMat(:,:)))) ; AVMDs(1,:)+nanstd(AVMDMat(:,:))./sqrt(sum(~isnan(AVMDMat(:,:))))];
AVs=squeeze(AVs);

if(LongTimeScale>1)
    TimeGridOld=[-timewidth:1:(timewidth*2-1)];
    TimeGridNew=[-timewidth:LongTimeScale:(timewidth*2-1)];
    AVs=interptest(TimeGridOld,AVs,TimeGridNew);
    AVMDs=interptest(TimeGridOld,AVMDs',TimeGridNew)';
    
    AVnnans=interptest(TimeGridOld,AVnnans',TimeGridNew)';
    NewAVMatBars(1,:,:)=interptest(TimeGridOld,squeeze(AVMatBars(1,:,:)),TimeGridNew); 
    NewAVMatBars(2,:,:)=interptest(TimeGridOld,squeeze(AVMatBars(2,:,:)),TimeGridNew);
    AVMatBars=NewAVMatBars;
    
    AVMDMatBars=[interptest(TimeGridOld,AVMDMatBars(1,:)',TimeGridNew)'; ...
        interptest(TimeGridOld,AVMDMatBars(2,:)',TimeGridNew)'];
    
    xa=(-timewidth:LongTimeScale:timewidth*2)./LongTimeScale;
    AVMDblock=nanmedian(reshape(AVMDMat(:,14:end),[],length(TimeGridNew)-1));


    AVMD2=AVMDs;
    %{
    AVMD2=zeros(length(starti),13);
    AVMD2(:,1)=nanmedian(AVMDMat(:,1:12),2);
    for i=1:length(TimeGridNew)-2
       AVMD2(:,i+1)=nanmedian(AVMDMat(:,i*LongTimeScale+1-12:i*LongTimeScale+12),2);
    end
    AVMD2(:,13)=nanmedian(AVMDMat(:,end-12:end),2);
    AVMD2=nanmedian(AVMD2);
    %}
end

%Shifting these to external functions for cleanliness and proof of independence/modularity
if(ccanalysis)
    CCAnalysis(AVMat,AVMDMat,satnum);
end

if(nnanalysis) 
    NNAnalysis(AVMat,AVMDMat,satnum);
end

if(arxanalysis)
    ARXAnalysis(FILLED,MassDensitySpline,headers,satnum);
end


save(savefilename)

else
    load(savefilename)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%    PLOT SECTION    %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Default x-axis for plots
xa=(-timewidth:LongTimeScale:timewidth*2-1);

%Plots just made to prove data is accurate, or make some simple comparison
DataProofPlots

%Plots that are a simple scatter plot comparison
BasicPlots

%Plots meant to reproduce a previously published figure
PlotReproductions

%Every plot left here is more advanced analysis
ComplexPlots

%profile viewer



%Statistical comparisons
if(stormcase==1)
    %Show DST is significantly different during missing data times
    [p,h]=ranksum(FILLED(~isnan(MassDensitySpline),15),FILLED(isnan(MassDensitySpline),15),'Alpha',0.01);
    [h,p]=ttest2(FILLED(~isnan(MassDensitySpline),15),FILLED(isnan(MassDensitySpline),15),'Alpha',0.01);
    %Or when mass density gets above/below 40
    [p,h]=ranksum(FILLED((MassDensitySpline<40),15),FILLED((MassDensitySpline>40),15),'Alpha',0.01);
    [h,p]=ttest2(FILLED((MassDensitySpline<40),15),FILLED((MassDensitySpline>40),15),'Alpha',0.01);
    %Or pre-noon vs post-noon
    [p,h]=ranksum(FILLED(FILLED(:,3)>12,15),FILLED(FILLED(:,3)<12,15),'Alpha',0.01);
    [h,p]=ttest2(FILLED(FILLED(:,3)>12,15),FILLED(FILLED(:,3)<12,15),'Alpha',0.01);
end



%Older plots, possibly for thesis, just batch comparing variables and
%making correlation tables
if(0)
    OldPlots
end