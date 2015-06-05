function Density

close all;clear all;
%load('../OMNI_OMNI2_merged')
%{
getyears=Year>=2000;
getyears=getyears+(Year<=2011);
getyears=(getyears==2);
VBS=VBS(getyears);
DST=Dst_index(getyears);
ION=Ion_density(getyears);
HOUR=Hour(getyears);
KP=Kp_index(getyears);
%}
FILLED=dlmread('WGhourFS_72_13.txt',',',1,0);
%FILLED=FILLED((FILLED(:,15)>-30),:); %DST less than -30

%Find a good storm

%FILLED=FILLED(FILLED(:,1)>1980,:); %Just to get into Denton time

startt=727547.333333; 
%DST Storm: 726479.750000;
%Mass Storm: 725858.250000 - 20
selectduration=22;
OMNITime=datenum(FILLED(:,1),1,0)+FILLED(:,2)+FILLED(:,3)/24;
starti=find(floor(OMNITime*100000)==floor(startt*100000)); %Rounding error?
%starti=find(OMNITime==startt);
%FILLED=FILLED(starti:starti+selectduration,:);
%{
%}

%FILLED=FILLED((FILLED(:,1)==2000),:); %Year 2000
%FILLED=FILLED((FILLED(:,2)==43),:); %February 12 (Day 31+12)

headers=textread('WGhourFS_72_13.txt','%s',28,'delimiter',',');
VBS=1/2*FILLED(:,6).*(abs(FILLED(:,5))-FILLED(:,5));
VBz=FILLED(:,6).*FILLED(:,5);
BS=1/2*(abs(FILLED(:,5))-FILLED(:,5));
OMNIDensity=FILLED(:,7);
OMNITime=datenum(FILLED(:,1),1,0)+FILLED(:,2)+FILLED(:,3)/24;
%VBS=1/2*Plasma_bulk_speed.*(abs(Bz_GSM)-Bz_GSM);
%OMNITime=Time;

if(exist('DentonDensityAndTime.mat','file'))
    load('DentonDensityAndTime')
else
    IN=dlmread('massdensity.txt');
    IN(IN(:,1)~=7,:)=[]; %Only use satellite 7 (6 should also work)
    %IN=sortrows(IN,2);
    IN(IN==9999)=NaN;
    
    [DentonTime,uniquerows]=unique(datenum(IN(:,2),1,0) + IN(:,3)+IN(:,4)/24+IN(:,5)/(24*60),'rows');
    
    %88 is density, 10 is f10.7, 31 is BZ_sw, 32 is Vsw
    F107=IN(uniquerows,10);
    DBz=IN(uniquerows,31);
    DVsw=IN(uniquerows,32);
    DBS=1/2*IN(uniquerows,32).*(abs(DBz)-DBz);
    MLT=IN(uniquerows,8);
    MassDensity=IN(uniquerows,88);
    
    save('DentonDensityAndTime','DentonTime','MassDensity','F107','DBz','MLT','DBS','DVsw');
end

if(~exist('figures','dir'))
    mkdir('figures');
end

%Test linear fit
%{
x=F107;
y=MassDensity;

p = polyfit(x,y,1)
yfit =  p(1) * x + p(2);
yresid = y - yfit;
rsq = 1 -sum(yresid.^2)/((length(y)-1) * var(y))

%Test without removing expanded plasma
t27=DentonTime(1):27:DentonTime(end); %Make 27 day time vector
y27=interptest(DentonTime,log(y),t27);
x27=interptest(DentonTime,x,t27);
p = polyfit(x27,y27,1)
yfit =  p(1) * x27 + p(2);
yresid = y27 - yfit;
rsq = 1 -sum(yresid.^2)/((length(y27)-1) * var(y27))

%Remove expanded plasma
x2=x(y<20);
y2=y(y<20);
DentonTime2=DentonTime(y<20);
y27=interptest(DentonTime2,log(y2),t27);
x27=interptest(DentonTime2,x2,t27);
p = polyfit(x27,y27,1)
yfit =  p(1) * x27 + p(2);
yresid = y27 - yfit;
rsq = 1 -sum(yresid.^2)/((length(y27)-1) * var(y27))

%Remove afternoon
prenoon=(MLT<12 & MLT>=7);
x2=x(y<25 & prenoon);
y2=y(y<25 & prenoon);
DentonTime2=DentonTime(y<25 & prenoon);
y27=interptest(DentonTime2,log(y2),t27);
x27=interptest(DentonTime2,x2,t27);
p = polyfit(x27,y27,1)
yfit =  p(1) * x27 + p(2);
yresid = y27 - yfit;
rsq = 1 -sum(yresid.^2)/((length(y27)-1) * var(y27))

ignore=1;%For breakpoint
%}

%Test subtracting F10.7 trend from mass density
%[~, cb, ~,xnew,corr] = IR(log(MassDensity),log(F107),0,12,0,0);
%MassDensity=exp(log(MassDensity)-xnew);


gettime=OMNITime>min(DentonTime);
gettime=gettime+(OMNITime<max(DentonTime));
gettime=(gettime==2);

OMNITime=OMNITime(gettime);
VBS=VBS(gettime);
VBz=VBz(gettime);
BS=BS(gettime);

OMNIDensity=OMNIDensity(gettime);
FILLED=FILLED(gettime,:);
OverlayFilled=FILLED;
OHr=FILLED(:,3);
DHr=hour(DentonTime);
%Rescale to overlay on plots of density (where there's space from 150-250)
for i=1:length(headers)
    OverlayFilled(:,i)=((250-150).*(OverlayFilled(:,i)-min(OverlayFilled(:,i))))./(max(OverlayFilled(:,i))-min(OverlayFilled(:,i))) + 150;
end
OverlayBS=((250-150).*(BS-min(BS)))./(max(BS)-min(BS)) + 150;
OverlayVBS=((250-150).*(VBS-min(VBS)))./(max(VBS)-min(VBS)) + 150;
OverlayVBz=((250-150).*(VBz-min(VBz)))./(max(VBz)-min(VBz)) + 150;
OverlayF107=((250-150).*(F107-min(F107)))./(max(F107)-min(F107)) + 150;
OverlayDBz=((250-150).*(DBz-min(DBz)))./(max(DBz)-min(DBz)) + 150;
OverlayDBS=((250-150).*(DBS-min(DBS)))./(max(DBS)-min(DBS)) + 150;

%MassDensitySpline=interp1(DentonTime,MassDensity,OMNITime,'linear');
if(exist('InterpedVals.mat','file'))
    load('InterpedVals')
else
    MassDensitySpline=interptest(DentonTime,MassDensity,OMNITime,(OMNITime(2)-OMNITime(1))/2);
    F107Spline=interptest(DentonTime,F107,OMNITime,(OMNITime(2)-OMNITime(1))/2);
    MassDensitySpline=MassDensitySpline';
    F107Spline=F107Spline';
    save('InterpedVals','MassDensitySpline','F107Spline');
end

FILLED=[FILLED F107Spline];
headers{end+1}='F107';

LongTimeScale=0;%24*27; %How many hours to average over
LongTimeInc=LongTimeScale*(OMNITime(2)-OMNITime(1));
if(LongTimeScale>0) 
   MassDensitySpline=interp1(OMNITime,MassDensitySpline,OMNITime(1):LongTimeInc:OMNITime(end));
   F107Spline=interp1(OMNITime,F107Spline,OMNITime(1):LongTimeInc:OMNITime(end));
   FILLED=interp1(OMNITime,FILLED,OMNITime(1):LongTimeInc:OMNITime(end));
   OMNITime=OMNITime(1):LongTimeInc:OMNITime(end);
end


%Find storm with enough data to analyze
MassDensityNanSpline=interp1(OMNITime(~isnan(MassDensitySpline)),MassDensitySpline(~isnan(MassDensitySpline)),OMNITime,'linear');
%storms=diff([0 (FILLED(:,15)<-50)' 0]); %DST Storm
storms=diff([0 (MassDensityNanSpline>40)' 0]); %Mass Density Storm, started at 40
%storms=[0 diff([0 (diff(MassDensityNanSpline)>10)' 0])];
%storms=[0 diff([0 (abs(MassDensityNanSpline(2:end)./MassDensityNanSpline(1:end-1))>1.3)' 0])];
starti=find(storms>0);
endi=find(storms<0)-1;
duration=endi-starti+1;

durationcaveat=''; %empty string if no cutoff set
cutoffduration=1; %If you want to get rid of short (spurious?) storms
if(cutoffduration>1)
    starti=starti(duration>cutoffduration);
    endi=endi(duration>cutoffduration);
    duration=duration(duration>cutoffduration);
    durationcaveat=sprintf('longer than %d hours',cutoffduration);
end

%Remove F10.7 influence
removef107=0;
if(removef107)
    [~, ~, ~,xnew,~] = IR(MassDensitySpline,F107Spline,0,1,0,0); %1 hour IR model
    MassDensitySplineOriginal=MassDensitySpline;
    MassDensitySpline=MassDensitySpline-xnew; 
end

stormi=1;
AVMat=[];
AVMDMat=[];
timewidth=24;
if(LongTimeScale>0),timewidth=4; end
maxwidth=timewidth*2;

while(starti(1)-maxwidth<1)
    starti(1)=[]; endi(1)=[]; 
end
while(endi(end)+maxwidth>length(MassDensitySpline))
    starti(end)=[]; endi(end)=[];
end
for i=1:length(starti)
   datanum=sum(~isnan(MassDensitySpline(starti(i):endi(i))));
   %if(datanum>(endi(i)-starti(i))/2 && datanum>18)
       AVMat(stormi,:,:)=FILLED((starti(i)-timewidth):starti(i)+timewidth*2,:);
       AVMDMat(stormi,:)=MassDensitySpline((starti(i)-timewidth):starti(i)+timewidth*2);
       stormi=stormi+1;
       %fprintf('%6f - %3d \n',OMNITime(starti(i)),endi(i)-starti(i))
  % end
end
AVs=nanmean(AVMat,1);
AVMDs=nanmean(AVMDMat);
AVnnans=sum(~isnan(AVMDMat));

%%%%%Make Plots?
MakePlots=0;
MakePaperPlots=1;
visible='off';

%Compare the two densities
if(MakePlots)
    h=figure('Visible',visible); plot(OMNITime,MassDensitySpline);hold on; plot(OMNITime,OMNIDensity,'r')
    legend('Denton','OMNI','Location','NorthEast')
    title('OMNI Density vs Denton Density')
    ylabel('Density')
    xlabel('Time')
    datetick
    print -dpng figures/densitycomp.png
end

if(MakePaperPlots && removef107)
   plot(OMNITime,MassDensitySplineOriginal-(MassDensitySpline-nanmean(MassDensitySplineOriginal)),'-') 
   title('F10.7 trend')
   ylabel('Difference (amu/cm^3)')
   xlabel('Time')
   datetick
   print -depsc2 -r200 paperfigures/f107removed.eps
end

if(MakePaperPlots)
h=figure('Visible',visible); 
orient tall;
h1=subplot('position',subplotstack(5,1));plot(OMNITime,FILLED(:,5),'.'); ylabel('Bz (nT)'); %Bz
h2=subplot('position',subplotstack(5,2));plot(OMNITime,FILLED(:,6),'.'); ylabel('V_{SW} (km/s)'); ylim([300,900]);%V_sw
h3=subplot('position',subplotstack(5,3));plot(OMNITime,FILLED(:,15),'.');yl=ylabel('D_{ST} (nT)'); set(yl,'Units','Normalized','Position',[-.03 0.5 0]); ylim([-300,100]);
h4=subplot('position',subplotstack(5,4));plot(OMNITime,FILLED(:,29),'.');ylabel('F10.7 (s.f.u.)');%f107
h5=subplot('position',subplotstack(5,5));plot(OMNITime,MassDensitySpline,'r.');ylabel('\rho_{eq} (amu/cm^3)');%f107
set(findobj('type','axes'),'xticklabel',{[]}); 
set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on')
datetick('x')
set(findobj('type','axes'),'xtick',get(h5,'xtick'))
axis tight;
        linkaxes([h5 h1 h2 h3 h4],'x')
        xlabel('Time')
        set(gcf,'NextPlot','add'); axes; h = title(sprintf('All data',length(duration)));set(gca,'Visible','off');set(h,'Visible','on'); 
        print -depsc2 -r200 paperfigures/alldata.eps
        
        
        xa=(1:length(AVMDs))-timewidth-1;
        h=figure('Visible',visible); 
        orient tall;
h2=subplot('position',subplotstack(5,2));plot(xa,AVs(1,:,6),'+-'); ylabel('V_{SW} (km/s)');%V_sw
hold on; plot(xa,AVs(1,:,6)+std(AVMat(:,:,6)),'r-.'); plot(xa,AVs(1,:,6)-std(AVMat(:,:,6)),'r-.');
h1=subplot('position',subplotstack(5,1));plot(xa,AVs(1,:,5),'+-'); ylabel('Bz (nT)'); %Bz
hold on; plot(xa,AVs(1,:,5)+std(AVMat(:,:,5)),'r-.'); plot(xa,AVs(1,:,5)-std(AVMat(:,:,5)),'r-.');
h3=subplot('position',subplotstack(5,3));plot(xa,AVs(1,:,15),'+-');ylabel('D_{ST} (nT)'); %dst
hold on; plot(xa,AVs(1,:,15)+std(AVMat(:,:,15)),'r-.'); plot(xa,AVs(1,:,15)-std(AVMat(:,:,15)),'r-.');
h4=subplot('position',subplotstack(5,4));plot(xa,AVs(1,:,29),'+-');ylabel('F10.7 (s.f.u.)');%f107
hold on; plot(xa,AVs(1,:,29)+nanstd(AVMat(:,:,29)),'r-.'); plot(xa,AVs(1,:,29)-nanstd(AVMat(:,:,29)),'r-.');
set(findobj('type','axes'),'xticklabel',{[]})
        xv=[xa(1) xa(end)];
        subplot('position',subplotstack(5,5)); [AX,H5,H6]=plotyy(xa,AVMDs(1,:),xa,AVnnans,'plot','bar'); 
        hold on; plot(xa,AVMDs(1,:)+nanstd(AVMDMat(:,:)),'r-.'); plot(xa,AVMDs(1,:)-nanstd(AVMDMat(:,:)),'r-.');
        set(AX(2),'Xlim',xv); set(AX(1),'Xlim',xv);  set(H5,'marker','+','color','red'); set(AX(1),'YColor','r'); set(AX(2),'YColor',[0 0.5 0.5]); set(get(H6,'child'),'FaceColor',[0 0.5 0.5]); uistack(AX(1)); set(AX(1),'Color','none'); set(AX(2),'Color','w');
        ylabel(AX(1),'\rho_{eq} (amu/cm^3)'); ylabel(AX(2),'non-nan datapoints');
        set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[-timewidth:timewidth/2:timewidth*2])
        linkaxes([AX h1 h2 h3 h4],'x')
        xlabel('Time from storm onset (27 day)')
        set(gcf,'NextPlot','add'); axes; h = title(sprintf('Average of %d storms %s (%d to %d)',length(duration),durationcaveat, year(OMNITime(1)),year(OMNITime(end))));set(gca,'Visible','off');set(h,'Visible','on'); 
        print -depsc2 -r200 paperfigures/stormavs-1.eps
        
end

if(MakePaperPlots)
    h=figure('Visible',visible); 
    plot(OMNITime,MassDensitySpline); hold on;
    plot([OMNITime(1) OMNITime(end)],[40 40],'r-.','LineWidth',6);
    ylabel('Mass Density')
    xlabel('Time')
    datetick('x')
    title('Mass Density Storms')
    print -depsc2 -r200 paperfigures/massdensitystorms.eps
    
    h=figure('Visible',visible); 
    plot(OMNITime,FILLED(:,15)); hold on;
    plot([OMNITime(1) OMNITime(end)],[-50 -50],'r-.','LineWidth',6);
    ylabel('D_{st}')
    xlabel('Time')
    datetick('x')
    title('D_{st} Storms')
    print -depsc2 -r200 paperfigures/dststorms.eps
end

if(MakePlots) %Stack plots
    for i=1:length(headers)
    subplot('position',subplotstack(5,mod(i,4)+1));plot(AVs(1,:,i),'+-');
    ylabel(headers{i})
    if(mod(i,4)==0 || i==length(headers))
        subplot('position',subplotstack(5,5)); [AX,H1,H2]=plotyy(1:length(AVMDs),AVMDs(1,:),1:length(AVMDs),AVnnans,'plot','bar'); 
        set(AX(2),'Xlim',[1 length(AVMDs)]); set(AX(1),'Xlim',[1 length(AVMDs)]); set(get(H2,'child'),'facea',.3); set(H1,'marker','+','color','red'); set(AX(1),'YColor','r');
        %drawnow; hold on; xlim manual; h=bar(AVnnans); ch = get(h,'child'); set(ch,'facea',.3)
        ylabel(AX(1),'Mass Density'); ylabel(AX(2),'non-nan datapoints');
        %xlabel('Time before storm onset (hr)')
        print('-dpng',sprintf('figures/stormavs-%d',i));
        close all
    end
    end
end

if(MakePlots)
    NanPH=zeros(1,24);
    for i=0:23
   NanPH(i+1)=sum(isnan(MassDensitySpline(hour(OMNITime)==i)));
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
        h=figure('Visible',visible);plot(OMNITime,FILLED(:,i))
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
    if(MakePlots), densityplot(OMNITime,[x,xnew1,OverlayFilled(:,i)],headers{i},Na,Nb1,corr1,pe1,visible), end
    pe=pe_nonflag(x,xnew);
    BigTable=[BigTable;{headers{i},corr1,corr,pe1,pe}];
    if(MakePlots), densityplot(OMNITime,[x,xnew,OverlayFilled(:,i)],headers{i},Na,Nb,corr,pe,visible), end
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
