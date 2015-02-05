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
FILLED=FILLED((FILLED(:,3)>10),:); %Hour greater than 10

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
    IN(IN(:,1)~=6,:)=[]; %Only use satellite 7 (6 should also work)
    %IN=sortrows(IN,2);
    IN(IN==9999)=NaN;
    
    [DentonTime,uniquerows]=unique(datenum(IN(:,2),1,0) + IN(:,3)+IN(:,4)/24+IN(:,5)/(24*60),'rows');
    
    %88 is density, 10 is f10.7, 31 is BZ_sw, 32 is Vsw
    F107=IN(uniquerows,10);
    DBz=IN(uniquerows,31);
    DBS=1/2*IN(uniquerows,32).*(abs(DBz)-DBz);
    MLT=IN(uniquerows,8);
    MassDensity=IN(uniquerows,88);
    
    save('DentonDensityAndTime','DentonTime','MassDensity','F107','DBz','MLT','DBS');
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
MassDensitySpline=interptest(DentonTime,MassDensity,OMNITime,(OMNITime(2)-OMNITime(1))/2);
MassDensitySpline=MassDensitySpline';


%Compare the two densities
figure; plot(OMNITime,MassDensitySpline);hold on; plot(OMNITime,OMNIDensity,'r')
legend('Denton','OMNI','Location','NorthEast')
title('OMNI Density vs Denton Density')
ylabel('Density')
xlabel('Time')
datetick
print -dpng figures/densitycomp.png

Na=0;
lag=0;
advance=0;

%x=detrend(MassDensitySpline);
x=MassDensitySpline;


%{ 
%Figure out what's going on with the year
[ca, cb, cc,xnew,corr] = IRboot(x,FILLED(:,1),Na,Nb,lag,advance);
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

MakePlots=0;
visible='off';

if(MakePlots)
    Hr=hour(DentonTime);
    for i=0:23
        avf107(i+1)=mean(F107(Hr==i));
        avdens(i+1)=mean(MassDensity(Hr==i));
    end
    avf107=(avf107-min(avf107))/(max(avf107)-min(avf107));
    avdens=(avdens-min(avdens))/(max(avdens)-min(avdens));
    h=figure('Visible',visible);
    plot(0:23,avf107,'r+:')
    hold on
    plot(0:23,avdens,'b+:')
    legend('F107','Density')
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

%x=detrend(MassDensitySpline);

%Open the table for writing
table=fopen('table.txt','w');
fprintf(table,'<pre>\n');
Nb1=1;
Nb=120;
x=MassDensitySpline;
fprintf(table,'Input \t CC(1) \t CC(120) \t PE(1) \t PE(120)\n');
BigTable={};
for i=1:length(headers)
    f=FILLED(:,i);
    [~,~,~,xnew,corr1] = IRboot(x,f,Na,Nb1,lag,advance); 
    pe1=pe_nonflag(x,xnew);
    if(MakePlots), densityplot(OMNITime,[x,xnew,OverlayFilled(:,i)],headers{i},Na,Nb1,corr1,pe1), end
    [~, cb, ~,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
    pe=pe_nonflag(x,xnew);
    BigTable=[BigTable;{headers{i},corr1,corr,pe1,pe}];
    if(MakePlots), densityplot(OMNITime,[x,xnew,OverlayFilled(:,i)],headers{i},Na,Nb,corr,pe), end
    if(MakePlots), densitycoefplot(-advance:Nb-advance-1,flipud(cb),headers{i},Na,Nb,corr,pe), end
end

x=MassDensity;
for i=1:length(dheaders)
    f=DFILLED(:,i);
    [~,~,~,xnew,corr1] = IRboot(x,f,Na,Nb1,lag,advance);
    pe1=pe_nonflag(x,xnew);
    if(MakePlots), densityplot(OMNITime,[x,xnew,DOverlayFilled(:,i)],dheaders{i},Na,Nb1,corr1,pe1), end
    [~, cb, ~,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
    pe=pe_nonflag(x,xnew);
    BigTable=[BigTable;{dheaders{i},corr1,corr,pe1,pe}];
    if(MakePlots), densityplot(OMNITime,[x,xnew,OverlayFilled(:,i)],headers{i},Na,Nb,corr,pe), end
    if(MakePlots), densitycoefplot(-advance:Nb-advance-1,flipud(cb),headers{i},Na,Nb,corr,pe), end
end

%Now for doubles
x=MassDensitySpline;
f=[FILLED(:,3) FILLED(:,5)]; %Hr and Bz
[~,~,~,xnew,corr1] = IRboot(x,f,Na,Nb1,lag,advance); 
    pe1=pe_nonflag(x,xnew);
    [~, ~, ~,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
    pe=pe_nonflag(x,xnew);
    BigTable=[BigTable;{'Hr+Bz',corr1,corr,pe1,pe}];
f=[FILLED(:,5) FILLED(:,6)]; %Bz and V
[~,~,~,xnew,corr1] = IRboot(x,f,Na,Nb1,lag,advance); 
    pe1=pe_nonflag(x,xnew);
    [~, ~, ~,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
    pe=pe_nonflag(x,xnew);
    BigTable=[BigTable;{'Bz+V',corr1,corr,pe1,pe}];




%{
[ca, ca2, cb, cb2, cc, cc2, xnew, xnew2, corr, corr2] = IRsortboot(MassDensity,FILLED(:,5),FILLED(:,5),100,Na,1,lag,advance);
pe=pe_nonflag(MassDensity,xnew);
pe2=pe_nonflag(MassDensity,xnew2);
densityplot(DentonTime,[MassDensity,xnew,OverlayFilled(:,5)],'BzHigh',Na,Nb,corr,pe)
densityplot(DentonTime,[MassDensity,xnew2,OverlayFilled(:,5)],'BzLow',Na,Nb,corr2,pe2)
densitycoefplot(-advance:Nb-advance-1,flipud(cb),'BzHigh',Na,Nb,corr,pe)
densitycoefplot(-advance:Nb-advance-1,flipud(cb2),'BzLow',Na,Nb,corr2,pe2)
%fprintf(table,'BZHigh \t -- \t %2.5f \t -- \t %2.5f\n',corr,pe);
%fprintf(table,'BZLow \t -- \t %2.5f \t -- \t %2.5f\n',corr2,pe2);
BigTable=[BigTable;{'BZHigh',NaN,corr,NaN,pe}];
BigTable=[BigTable;{'BZHigh',NaN,corr2,NaN,pe2}];
%}
BigTable=sortrows(BigTable,2);
for i=1:size(BigTable,1)
   if(BigTable{i,5}>-10)
       fprintf(table,'%s \t %2.2f \t %2.2f \t %2.2f \t %2.2f\n',BigTable{i,1},BigTable{i,2},BigTable{i,3},BigTable{i,4},BigTable{i,5}); 
   else
       fprintf(table,'%s \t %2.2f \t %2.2f \t %2.2f \t %2.0e\n',BigTable{i,1},BigTable{i,2},BigTable{i,3},BigTable{i,4},BigTable{i,5}); 
   end
end
fprintf(table,'\n</pre>');

fclose(table);
system('cat README.txt table.txt > README.md');


end

function densityplot(x,ys,string,Na,Nb,corr,eff)
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

function densitycoefplot(x,ys,string,Na,Nb,corr,eff)
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
