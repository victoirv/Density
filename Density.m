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
headers=textread('WGhourFS_72_13.txt','%s',28,'delimiter',',');
VBS=1/2*FILLED(:,6).*(abs(FILLED(:,5))-FILLED(:,5));
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



gettime=OMNITime>min(DentonTime);
gettime=gettime+(OMNITime<max(DentonTime));
gettime=(gettime==2);

OMNITime=OMNITime(gettime);
VBS=VBS(gettime);
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
Nb=1;
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

if(~exist(sprintf('figures/OMNI_%s.png',headers{1}),'file'))
    
    for i=1:length(headers)
        close all;
        figure;plot(OMNITime,FILLED(:,i))
        datetick
        ylabel(headers{i})
        xlabel('Time')
        title(sprintf('OMNI %s',headers{i}))
        filestring=sprintf('figures/OMNI_%s.png',headers{i});
        print('-dpng',filestring)
    end
    
end

corrs_1=1:(length(headers)+2);
pes_1=1:(length(headers)+2);
Nb=1;
for i=1:length(headers)
    f=FILLED(:,i);
    [ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
    corrs_1(i)=corr;  
    pe=pe_nonflag(x,xnew);
    pes_1(i)=pe;
    densityplot(OMNITime,[x,xnew,OverlayFilled(:,i)],headers{i},Na,Nb,corr,pe)
    
end

%Add VBS_1 and BS_1
%x=detrend(MassDensitySpline);

f=VBS;
[ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
corrs_1(end-1)=corr;
pes_1(end-1)=pe_nonflag(x,xnew);
densityplot(OMNITime,[x,xnew,OverlayVBS],'VBS',Na,Nb,corr,pe_nonflag(x,xnew))

f=BS;
[ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
corrs_1(end)=corr;
pes_1(end)=pe_nonflag(x,xnew);
densityplot(OMNITime,[x,xnew,OverlayBS],'BS',Na,Nb,corr,pe_nonflag(x,xnew))

%Open the table for writing
table=fopen('table.txt','w');
fprintf(table,'<pre>');
Nb=120;
fprintf(table,'Input \t CC(1) \t CC(120) \t PE(1) \t PE(120)\n');
BigTable={};
for i=1:length(headers)
    f=FILLED(:,i);
    [ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
    pe=pe_nonflag(x,xnew);
    %fprintf(table,'%s \t %2.5f \t %2.5f \t %2.5f \t %2.5f\n',headers{i},corrs_1(i),corr,pes_1(i),pe);
    BigTable=[BigTable;{headers{i},corrs_1(i),corr,pes_1(i),pe}];
    densityplot(OMNITime,[x,xnew,OverlayFilled(:,i)],headers{i},Na,Nb,corr,pe)
    densitycoefplot(-advance:Nb-advance-1,flipud(cb),headers{i},Na,Nb,corr,pe)
end

%Add VBS, BS, and F107 for 120 coef
f=VBS;

[ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
pe=pe_nonflag(x,xnew);
%fprintf(table,'VBS \t %2.5f \t %2.5f \t %2.5f \t %2.5f\n',corrs_1(end-1),corr,pes_1(end-1),pe);
BigTable=[BigTable;{'VBS',corrs_1(end-1),corr,pes_1(end-1),pe}];
densityplot(OMNITime,[x,xnew,OverlayVBS],'VBS',Na,Nb,corr,pe)
densitycoefplot(-advance:Nb-advance-1,flipud(cb),'VBS',Na,Nb,corr,pe)

f=BS;
[ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
pe=pe_nonflag(x,xnew);
%fprintf(table,'BS \t %2.5f \t %2.5f \t %2.5f \t %2.5f\n',corrs_1(end),corr,pes_1(end),pe);
BigTable=[BigTable;{'BS',corrs_1(end),corr,pes_1(end),pe}];
densityplot(OMNITime,[x,xnew,OverlayBS],'BS',Na,Nb,corr,pe)
densitycoefplot(-advance:Nb-advance-1,flipud(cb),'BS',Na,Nb,corr,pe)


[ca, cb, cc,xnew,corr1] = IRboot(MassDensity,DBS,Na,1,0,0);
pe1=pe_nonflag(MassDensity,xnew);
densityplot(DentonTime,[MassDensity,xnew,OverlayDBS],'DBS',Na,1,corr,pe1)
[ca, cb, cc,xnew,corr] = IRboot(MassDensity,DBS,Na,Nb,lag,advance);
pe=pe_nonflag(MassDensity,xnew);
densityplot(DentonTime,[MassDensity,xnew,OverlayDBS],'DBS',Na,Nb,corr,pe_nonflag(MassDensity,xnew))
densitycoefplot(-advance:Nb-advance-1,flipud(cb),'DBS',Na,Nb,corr,pe_nonflag(MassDensity,xnew))
%fprintf(table,'DBS \t %2.5f \t %2.5f \t %2.5f \t %2.5f\n',corr1,corr,pe1,pe);
BigTable=[BigTable;{'DBS',corr1,corr,pe1,pe}];


[ca, cb, cc,xnew,corr1] = IRboot(MassDensity,F107,Na,1,0,0);
pe1=pe_nonflag(MassDensity,xnew);
densityplot(DentonTime,[MassDensity,xnew,OverlayF107],'F107',Na,1,corr,pe1)
[ca, cb, cc,xnew,corr] = IRboot(MassDensity,F107,Na,Nb,lag,advance);
pe=pe_nonflag(MassDensity,xnew);
densityplot(DentonTime,[MassDensity,xnew,OverlayF107],'F107',Na,Nb,corr,pe_nonflag(MassDensity,xnew))
densitycoefplot(-advance:Nb-advance-1,flipud(cb),'F107',Na,Nb,corr,pe_nonflag(MassDensity,xnew))
%fprintf(table,'F107 \t %2.5f \t %2.5f \t %2.5f \t %2.5f\n',corr1,corr,pe1,pe);
BigTable=[BigTable;{'F107',corr1,corr,pe1,pe}];

[ca, cb, cc,xnew,corr1] = IRboot(MassDensity,log(F107),Na,1,0,0);
pe1=pe_nonflag(MassDensity,xnew);
densityplot(DentonTime,[MassDensity,xnew,OverlayF107],'lnF107',Na,1,corr,pe1)
[ca, cb, cc,xnew,corr] = IRboot(MassDensity,log(F107),Na,Nb,lag,advance);
pe=pe_nonflag(MassDensity,xnew);
densityplot(DentonTime,[MassDensity,xnew,OverlayF107],'lnF107',Na,Nb,corr,pe_nonflag(MassDensity,xnew))
densitycoefplot(-advance:Nb-advance-1,flipud(cb),'lnF107',Na,Nb,corr,pe_nonflag(MassDensity,xnew))
%fprintf(table,'logF107 \t %2.5f \t %2.5f \t %2.5f \t %2.5f\n',corr1,corr,pe1,pe);
BigTable=[BigTable;{'lnF107',corr1,corr,pe1,pe}];

[ca, cb, cc,xnew,corr1] = IRboot(MassDensity,DBz,Na,1,0,0);
pe1=pe_nonflag(MassDensity,xnew);
densityplot(DentonTime,[MassDensity,xnew,OverlayDBz],'DBz',Na,1,corr,pe1)
[ca, cb, cc,xnew,corr] = IRboot(MassDensity,DBz,Na,Nb,lag,advance);
pe=pe_nonflag(MassDensity,xnew);
densityplot(DentonTime,[MassDensity,xnew,OverlayDBz],'DBz',Na,Nb,corr,pe_nonflag(MassDensity,xnew))
densitycoefplot(-advance:Nb-advance-1,flipud(cb),'DBz',Na,Nb,corr,pe_nonflag(MassDensity,xnew))
%fprintf(table,'DBz \t %2.5f \t %2.5f \t %2.5f \t %2.5f\n',corr1,corr,pe1,pe);
BigTable=[BigTable;{'DBz',corr1,corr,pe1,pe}];

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
   fprintf(table,'%s \t %2.2f \t %2.2f \t %2.2f \t %2.2f\n',BigTable{i,1},BigTable{i,2},BigTable{i,3},BigTable{i,4},BigTable{i,5}); 
end
fprintf(table,'</pre>');

fclose(table);
system('cat README.txt table.txt > README.md');


end

function densityplot(x,ys,string,Na,Nb,corr,eff)
close all;
test=size(ys);
if test(1)>test(2) 
   ys=ys'; 
end
figure;plot(x,ys)
    legend('Data','Prediction',string,'Location','NorthEast')
    ylabel('Density')
    xlabel('Time')
    datetick
    title(sprintf('Denton Density from OMNI %s: Nx:%d Nf:%d, corr: %2.3f, eff: %2.3f',string,Na,Nb,corr,eff))
    filestring=sprintf('figures/density_%s_%d_%d.png',string,Na,Nb);
print('-dpng',filestring)
end

function densitycoefplot(x,ys,string,Na,Nb,corr,eff)
close all;
figure;plot(x,ys)
ylabel('Impulse Coefficient')
xlabel('Time Lags')
grid    
    title(sprintf('Denton Density from OMNI %s: Nx:%d Nf:%d, corr: %2.3f, eff: %2.3f',string,Na,Nb,corr,eff))
    filestring=sprintf('figures/density_%s_%d_%d_Cb.png',string,Na,Nb);
    print('-dpng',filestring)
end
